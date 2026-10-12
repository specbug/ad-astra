# Methods

This document describes how the search works, what each threshold is and where it came from, and
how false alarms were measured. Numbers quoted here come from the pilot runs described below; the
final survey statistics are in `results/`.

## 1. Data access
Everything is read from the public MAST bucket on AWS (`s3://stpubdata`, read over HTTPS without
credentials). The other astronomy archives (MAST API, ExoFOP, NASA Exoplanet Archive, SIMBAD, VizieR,
the Gaia archive, arXiv) were not reachable from the compute environment, which shaped several choices.

- **Light curves.** QLP full-frame light curves (sectors 1–98; MIT Quick-Look Pipeline), TESS-SPOC
  full-frame light curves (1–81) and SPOC 2-minute light curves (1–107). For QLP we use the
  systematics-corrected flux `SYS_RM_FLUX` and our own detrending, not QLP's detrended `DET_FLUX`:
  QLP's spline detrending erases dips longer than about half a day (pilot example: TIC 229434519's
  6%-deep, 0.8-day dip is present in the raw flux but flattened away in `DET_FLUX`). Only cadences
  with QLP `QUALITY == 0` are used; for SPOC products the standard lightkurve quality bitmask.
- **Stellar parameters.** TESS Input Catalog v8.2 (HATS parquet copy in the bucket): 23.3 million
  stars with Tmag < 13.5 or distance < 200 pc were extracted by fetching only the needed parquet column
  chunks with HTTP range requests (`pipeline/s3parquet.py`, `build_tic.py`).
- **Pixels.** Full-frame-image cutouts are cut straight out of the per-sector image cubes in the
  bucket with `astrocut` (about 10 s per 15×15-pixel cutout).
- **Gaia DR3** neighbours come from the HATS copy of `gaia_source` in the same bucket.
- **Known objects.** ExoFOP TOI table export of February 2026 (TOIs up to 7689; from the public
  `shporer/autottv` repository), the NASA Exoplanet Archive TOI table snapshot of 5 October 2026
  (positions and periods; from a public Hugging Face mirror), and SPOC threshold-crossing events read
  from the data-validation (DV) XML reports stored next to the 2-minute light curves.
  Not available offline: community TOIs (CTOIs), the AAVSO variable star index, the TESS-SPOC
  full-frame TCE lists, and the literature beyond web search.

## 2. Survey (Tracks B, C, D)
Sample: 1.95 million stars, all dwarfs or subgiants (TIC radius ≤ 2.5 R☉) with Tmag ≤ 12, plus every M
dwarf (Teff < 3900 K, R < 0.65 R☉) with Tmag ≤ 13.5; 10.2 million QLP light curves in sectors 1–98, and
all SPOC 2-minute light curves of sectors 99–107 (which have no full-frame light curves in the bucket yet).

Per light curve (`survey.py`, `dips.py`):
1. Bin to 10 minutes; normalise each contiguous segment; drop segments shorter than 0.75 d.
2. Flatten twice with a Tukey biweight (c = 5, as in `wotan`) evaluated on a knot grid: a 0.75-day
   window for dips up to 4.5 h and a 2.5-day window for longer dips. Upward outliers (> 4σ) are clipped;
   downward points never are.
3. Box filter at 10 durations from 0.75 h to 18 h. The noise for each duration is the robust spread of
   the box-filtered series itself, so correlated (red) noise on that timescale is included. White-noise
   error bars overstate the significance of hours-long dips by factors of 2–3 in TESS data.
4. Every local minimum with depth/noise ≥ 6 is an event; events at several durations are merged.
   Each event carries: depth in the raw SAP flux and (QLP) in QLP's detrended small- and
   large-aperture fluxes, background and centroid excursions, local noise within ±3 d, distance to
   data gaps, coverage, and simple shape numbers.
5. 10% of light curves (25% in pilots) are also searched upside-down (flux flipped about its median),
   through the identical pipeline. Anything found there is a false alarm by construction.

Per sector (`survey_reduce.py`):
- **Systematics.** For each camera/CCD, the number of stars with a dip in each 0.1-day bin is compared
  with the CCD's median; bins at ≥ 5× the median (and Poisson z > 4) are flagged, plus the bins either
  side. The same at 0.05 d across all cameras catches spacecraft events. In the sector-80 pilot 83%
  of all raw dips fell in flagged bins.
- **Blends.** Dips within half a duration of each other in stars < 2.5′ apart: only the deepest is kept
  as a possible source.
- **Clean dips** must also: be ≥ max(0.25 d, one duration) from a data edge, have ≥ 80% coverage, not
  coincide with a background rise (z_bkg ≤ 4), appear in the raw SAP flux at ≥ 40% of the flattened depth,
  and have local SNR ≥ 5.

## 3. Track B: single and duo transits
1. Candidate dips: clean, SNR ≥ 8 (local ≥ 6), duration ≥ 1 h, ≤ 2 clean dips in that light curve,
   implied radius < 2.2 R_Jup (TIC radius, contamination-corrected), dwarf host.
2. **Stage 2** (`vet_events.py`, `stage2.py`) re-fetches the light curve and measures:
   - dip-vs-bump: the strongest *upward* excursion of the same duration elsewhere in the light curve
     (noise and variability make bumps as well as dips; a transit only makes a dip);
   - a local fit of quadratic baseline × trapezoid, against baseline-only, an asymmetric trapezoid and
     a Gaussian dip: depth, T14, ingress fraction, asymmetry, Δχ²;
   - the same fit on raw SAP flux.
   *Planet-like* requires all of: strongest bump < 0.5× the dip SNR; no other dip > 0.5× unless it repeats;
   ingress fraction ≤ 0.5 (flat-bottomed); asymmetry Δχ² < max(10, 2% of the transit Δχ²); SAP depth
   within 0.6–1.6× of the flattened depth; transit Δχ² > 100; reduced χ² < 2; radius < 2.2 R_Jup.
   Validation: in the sector-80 pilot the only planet-like survivor out of 69 single dips was
   TOI-2583.01, a confirmed planet whose two transits the survey found blind.
3. **Star-level check** (`starcheck.py`, `finalize_b.py`): the dip finder is run on every sector of the
   star. Deeper dips anywhere mean an eclipsing binary (pilot example: TIC 259128745's "transit" is the
   6% secondary eclipse of an 80.5-day binary with 39% primary eclipses). Consistent repeats give a
   period via the longest ephemeris fitting all events within a quarter of a duration.
4. **Known?** TOI by TIC ID and by position (150″, period-matched, to catch signals attributed to a
   neighbour), SPOC TCEs from DV reports, web search of identifiers.
5. **Pixels** (`pixels.py`): difference image (out-of-transit minus in-transit) and its flux-weighted
   centroid, compared with the direct image's centroid measured identically (this cancels the
   pixel-phase bias; for WASP-126 b the offset is 0.04 px), plus Gaia DR3 neighbours overlaid.

## 4. Track C: eclipsing triples and 2+2 quadruples
For light curves with ≥ 3 clean dips of SNR ≥ 10 in a sector (`ebsearch.py`): BLS period on a
geometric frequency grid, P vs 2P from odd/even depths, a phase-binned median template is subtracted,
eclipse phases are masked, and the residuals are searched twice more for periodic signals and once
for isolated extra dips. A second signal counts if SNR ≥ 15, ≥ 2 eclipses, its period is not within 2%
of any ratio a/b (a, b ≤ 6) of the main period, and it is eclipse-shaped: box depth / peak-to-trough of
a sinusoid at the same period > 3 (pulsations and rotation give ≈ 1). Validation: the known sextuple
TIC 168789840 (Powell et al. 2021) is recovered blind with all three of its eclipse periods
(1.570, 8.22 and 1.306 d) in every sector tested. Candidates recurring in several sectors with the
same second period are ranked first; difference images are made separately for each eclipse set.

## 5. Track A: small planets around under-searched K/M dwarfs
Targets (`track_a_targets.py`): K/M dwarfs (2700–5400 K, R < 1 R☉, log g > 4.2) within 150 pc,
contamination < 0.2, ≥ 3 sectors spanning ≥ 300 d, at most two of them at 2-minute cadence, not TOI hosts:
100,146 stars; the 1,494 with the highest expected SNR for a 1.5 R⊕ planet on a 10-day orbit (≥ 7, using
the measured QLP noise) are searched. All sectors are stitched (SPOC 2-min > TESS-SPOC > QLP per sector),
flattened with a 0.75-day biweight, binned to 30 min, and searched with an iterative BLS (0.5–30 d; grid
fine enough that a grazing transit stays within half a duration across the full multi-year baseline).
Half of the stars are also searched flipped. Red flags follow the reference project
(paraschopra/tess-planet-candidates): odd/even, secondary eclipse, duration vs stellar density,
single-event dominance, edges, per-sector consistency, broad depression, depth, variability repeat,
common period across stars; plus empirical (red) noise from 40 off-transit box placements.

## 6. Deep vetting of individual candidates (`vet_deep.py`)
Run on every Track B shortlist star brighter than T = 10.5 (the shortlist: not a TOI, Rp ≤ 16 R⊕, single or
repeating on periods ≥ 15 d, flagged in fewer than five SPOC searches):
- **Star:** TIC parameters; Gaia DR3 RUWE, non-single-star flag and radial-velocity scatter of the target.
- **Known:** TOIs by TIC ID or within 10′ (ExoFOP February 2026 and NASA 5 October 2026 tables), SPOC TCEs.
- **Events:** the depth of every listed dip in its sector's best light curve (SPOC 2-min > TESS-SPOC > QLP).
- **Period:** for two or more events, the longest period consistent with all of them and every alias P/n
  (n ≤ 12) tested against all sectors: a predicted transit inside the data that is not there at > 5σ excludes
  that alias. For one event, the shortest period the other sectors allow (`duo.single_period_constraints`).
- **Pixels (source test):** for each transit, the difference image (out-of-transit minus in-transit, with a
  linear baseline per pixel) is fitted as a scaled copy of the out-of-transit image shifted to each candidate
  source position: the target and every Gaia star within 150″ bright enough to produce the dip with an eclipse of
  at most 50%. The χ² difference (alternative − target), summed over transits, is reported; positive favours the
  target. 2-minute target pixel files are used where they exist, otherwise 15 × 15 full-frame cutouts.
  Full-frame cutouts give much larger |Δχ²| than 2-minute pixels (more cadences, and residual pointing structure
  in the difference image); only the sign and the 2-minute values are read as evidence.

Verdicts: *off-target* if any source beats the target by Δχ² > 9; *binary host* if RUWE > 1.4 or Gaia has a
non-single-star solution; *no consistent period* if no alias fits all events; *published* if the literature
check (below) finds the star in this role; otherwise *candidate*.

## 7. Transit fits (`transit_fit.py`)
SPOC 2-minute PDCSAP windows (±0.6–1 d) are normalised by a quadratic fitted outside the transit and fitted
jointly with `batman` (quadratic limb darkening fixed to TESS-band values for the star's Teff) and `emcee`
(40 walkers, 3000 steps, first third discarded). Periodic candidates are parameterised by reference mid-time,
period, Rp/R★, impact parameter and log stellar density, so the fitted density can be compared with the TIC
value; duos by Rp/R★, b, T14 and one mid-time per transit.

## 8. Literature checks
TOI tables (above), SPOC DV reports, and web searches of the TIC, HIP and Gaia DR3 identifiers (one identifier
per query; a result counts only if the page names the same star). ExoFOP community TOIs, SIMBAD, VizieR and
arXiv listings could not be queried from the compute environment, so "no mention found" is weaker than "not
catalogued".

On 11 October the web searches were extended to every star the report names (the headline candidates had
only been checked against the TOI tables and SPOC). That found HIP 71639 (TIC 128703021) in the Planet Hunters
TESS candidate list: a single sector-11 transit, published without a period by Eisner et al. (2021, MNRAS 501,
4669) and sent to ExoFOP as a community TOI. The report now credits that detection; the new contribution is the
23.884-day period from four transits. Each entry in `results/candidates/literature_checks.json` has one of three
statuses: "known in this role" (named as the same kind of object), "mentioned elsewhere" (named in another
role, e.g. as a radial-velocity planet host) and "no mentions found".

Of the 150 stars the report names, 10 are known in their role: HIP 71639 (above); three Track B transit
candidates that are confirmed planets or brown dwarfs found from TESS transits by ground-based follow-up and
never made TOIs (NGTS-38 b = TIC 65910228 b, 180.5 d; TIC 279401253 b, 76.8 d; TIC 4672985 b, 69.0 d); one
published machine-learning candidate (TIC 221567884); and five Track C stars (the sextuple TIC 168789840, the
pulsating binary TIC 309658221, the catalogued quadruple TIC 321471064, the doubly eclipsing TIC 389836747, and
TIC 293954660, a published doubly eclipsing pair that turned out to be a blend). Gaia DR3 identifiers never
produced a hit; TIC numbers did. Snippets of large candidate tables (TESS multiply eclipsing lists, the T16
hunt's 10,000 Cycle 1 candidates) do not show every row, so those lists remain unchecked row by row.

The search budget is 200 web searches per turn, shared by all subagents. The first parallel sweep hit it
part-way, and the remainder was searched in a later turn.

Confirmed planets: `known.flag_confirmed` matches candidates within 10″, or by a "TIC <id>" host name, against
the NASA Exoplanet Archive confirmed-planet table (ps, default_flag = 1) of 5 October 2026. That table is copied
from the `juliensimon/nasa-exoplanets` Hugging Face mirror into `catalogs/hf_ps/` page by page through the
Hugging Face connector, the only route to it from the compute environment. `match_confirmed.py` adds the
match to every Track B, C and D table and writes `results/candidates/confirmed_matches.csv`.

This found what the web searches could not, because the literature uses other names. Nine confirmed systems
without TOIs were among Track B's candidates, each found blind, with the published period where the survey had
one:
- NGTS-38 b, 180.5 d, one transit
- TIC 279401253 b, 76.8 d, one transit
- TIC 4672985 b, 69.0 d, one transit
- NGTS-39 b, 58.2 d
- K2-290 c, 48.4 d, one transit (TIC 70803960)
- HD 148797 b and c, 42.1 and 68.2 d
- HD 3167 c, 29.85 d (HIP 2736, one transit; the 0.96-day SPOC detection is HD 3167 b)
- HD 278555 b, 28.84 d
- TIC 65612701 b, 10.22 d

Three of them, HD 3167, HD 278555 and K2-290, had been presented as new bright candidates; the report now marks
them "published". Two more Track B matches are circumbinary planets around eclipsing binaries the survey had
already rejected as binaries (Kepler-16, BEBOP-3). HD 11506 matches its three known radial-velocity planets. One confirmed planet was lost: all nine transits of TIC 245076932 b (21.6 d) were found, then
vetoed as an eclipsing binary because deeper dips elsewhere in its data had not been checked.

## 8b. Direct catalogue checks (12 October)
Once the environment's network allowed them, every star the report names (136) was checked directly:
- **ExoFOP** (`exofop_check.py`): the TOI and CTOI tables of 12 October 2026, and each star's target page
  (TOIs, CTOIs, every user-uploaded set of planet parameters with its period, author and date, imaging,
  spectra). The CTOI table was also matched against every Track B, C and D star
  (`exofop_ctoi_matches.csv`: 325 stars).
- **SIMBAD** (TAP, 10″ cone, nearest object): main identifier, object types, identifiers and the number of
  references (`simbad_named.csv`). Binary flags count only within 3″; farther ones are reported as a
  neighbouring binary.
- **VizieR** (10″): the TESS eclipsing-binary catalogue (Prša et al. 2022), Gaia DR3 eclipsing binaries, VSX,
  Kostov et al. 2022 quadruple candidates, Zasche et al. 2022 multiply eclipsing and doubly eclipsing systems,
  Zasche et al. 2025, Borkovits et al. 2025 (`vizier_named.csv`).
- **NASA Exoplanet Archive**: confirmed transiting planets within 2′ of each Track B star whose predicted
  transits coincide with the survey's dips (`neighbour_planet_transits.csv`). This catches blends from a
  neighbour's planet.
- **ADS** (`ads_check.py`, run in a worker session holding the API key): full-text search by every identifier
  of each named star (`ads_mentions.json`).

`prior.py` gathers all of this into one label per star (blend, confirmed, community candidate, published
multiple, catalogued binary, or nothing on record), and the report uses that label everywhere.

What it changed:
- **HIP 71639.** Another ExoFOP user posted its 23.8843518-day period on 15 May 2026, from the same four
  transits. The survey's period agrees within 0.1σ, so it is an independent confirmation, not a new period.
- **HIP 56966.** A single-transit CTOI since March 2025 (Salinas et al.), with no period. The other three
  transits and the 19.826-day period appear nowhere else.
- **HD 11506.** The sector-30 transit has been a CTOI since August 2024 (TESS Single Transit Planet Candidate
  working group), listed at 7,210 ± 3,331 ppm and 26.9 ± 7.4 h; the 2-minute data give 668 ppm and 3.2 h. The
  second transit (August 2026) appears nowhere else.
- **Bright deep-vetted list.** Of the 45 stars that pass, 17 are on record and one dip (TIC 457134361 =
  HD 80607) is a transit of HD 80606 b on the neighbouring star, 21″ away. 27 have nothing on record.
- **Track B shortlist.** 86 of 743 are on record, mostly single-transit CTOIs (TESS Single Transit Planet
  Candidate working group; Salinas et al. 2025).
- **Track C table.** 7 of 30 are published multiples or a known blend. 9 are catalogued binaries on the star
  itself and 3 have one 5–6″ away. 11 have nothing on record.
- **ADS.** 22 of the 136 stars are named in papers, all of them already known in other ways (confirmed planets,
  the HD 80606 system, HD 11506's radial-velocity planets, published quadruple candidates and eclipsing binaries)
  or named only incidentally (a 1933 position catalogue, a 1988 paper on another star). No paper names HIP 56966
  under any of its eight identifiers. None of HD 11506's 38 papers reports a transit; a targeted search for
  HD 11506 with "transit" or "TESS" returns only radial-velocity and statistics papers. HIP 71639 appears only in
  Planet Hunters TESS II.

## 9. Track C eclipse-shape test, revised (10 October)
The first version compared the box depth with the peak-to-trough of a single sine wave at the second period.
A pulsation at half that period has no power at the period itself, so it passed (the published pulsating
binary TIC 309658221 did). The test now fits a three-harmonic Fourier series at the period
(`ebsearch.box_over_fourier`): smooth signals, including pulsations at P/2 and P/3, give a ratio of about 0.5;
box-shaped eclipses with duty cycles of 3–10% give 1.2–3.3. Second signals need a ratio above 1.0.

## 10. Calibration results (interim)
- **Flipped light curves (false alarms by construction):** planet-shaped dips longer than 2 h on stars hotter
  than 4000 K: none in 11,843 flipped two-minute light curves; 0.97 per 1000 flipped full-frame light curves,
  concentrated in QLP sectors 20–22, 90 and 91 (up to 3.3 per 1000), against 0.33 per 1000 real full-frame light
  curves. Because only 10% of light curves are flipped, the CCD-wide systematics veto is ten times weaker for
  them, so these rates are upper bounds.
- **Blind recovery of known TOIs** (period > 15 d, depth > 1.5 ppt, predicted transits in surveyed full-frame
  sectors): 14% of 1.5–3 ppt transits, about 40% of deeper ones, listed as single-transit candidates. Losses:
  data gaps, ephemeris drift since discovery, and deliberately strict cuts.
- **Independent recoveries:** seven Track B "new" stars matched TOIs created in 2026 (TOI-7701 to 7850), after
  the February TOI list the survey was first compared with.

## 11. Extension (12 October): the newest TESS images, other telescopes, a period hunt

### 11.1 Our own light curves of sectors 99–108 (`ffi_new.py`)
Public full-frame light curves end at sector 98 (QLP) and near sector 80 (TESS-SPOC), but the calibrated
full-frame image cubes of sectors 99–108 are on MAST's AWS bucket. For a target on silicon in those sectors
(`tess-point`) without SPOC 2-minute data there, a 21 × 21-pixel cutout (astrocut) gives a light curve:
background per cadence from the faintest 40% of pixels away from the star; an aperture of the pixels within
r = 1–4 px of the star, choosing the r with the lowest 1-hour scatter after a 1-day running median; SPOC quality
flags and scattered-light glints (background far above its running median) removed. 261 light curves were made
(period-hunt targets, the 28 bright single transits, the headline stars). Typical 1-hour scatter: 200–600 ppm at
T = 9–11.

### 11.2 Period hunt (`period_hunt.py`, `period_review.py`, `period_vet.py`, `period_deep.py`)
Targets: the 1003 TOIs and CTOIs (ExoFOP, 12 October 2026) with no period, minus false positives. Every sector is
loaded (SPOC 2-min > TESS-SPOC > QLP, plus §11.1 light curves in a second pass), flattened and binned to 30
minutes; transits of other known planets on the star are masked; the catalogued transit is re-measured (duration
and time from the box SNR, because catalogued durations are sometimes far off); the rest of the data is searched
with a box of that duration for dips of 0.5–2× the depth at ≥ 7σ with data on both sides. The alias test
(`duo.alias_table`, vectorised on 12 October: identical results, 100× faster) lists every period that puts all
events on one ephemeris without predicting a transit the data do not show (depth tested at 0.6× the measured one).
Review: each new dip is re-measured on its own sector's light curve against a local quadratic baseline, with a
red-noise term, and kept only if SNR ≥ 7, ≥ 70% covered with data on both sides, ≥ 0.5 d from sector edges and
gaps, no background jump (> 5 MAD), at least half as deep in SAP flux, and 0.5–2× the catalogued depth. On
record: NASA archive planets within 10″ and every period uploaded to the star's ExoFOP page. New periods (on no
list, three or fewer aliases) then get the pixel test (§11.3) and a transit fit at each alias with the implied
stellar density compared with the TIC value.

### 11.3 Pixel test, revised: calibrated difference-image centroids (`eb_pixels.py`, `centroid_vet.py`)
The source test of §6 (difference image fitted as a shifted copy of the whole out-of-transit image) proved
unreliable in crowded fields: on Kostov et al.'s quadruples and against ZTF it called on-target systems blends and
blends on-target, because the shifted image moves every star, not just one. It is replaced, for every verdict in
the report, by the standard centroid test:
- per transit or eclipse, the difference image (linear fit to the frames just outside, minus the mean in transit)
  is fitted with a Gaussian of free position, width, amplitude and offset;
- the position is measured from where the target really is on that frame: the out-of-transit image is fitted with
  Gaussians at the Gaia DR3 positions moved to the observing date with their proper motions (fluxes in Gaia's
  proportions, one common shift, width, scale and background), which removes errors of the frame's coordinate
  solution (often several arcseconds) and the target's own motion;
- eclipses are combined within a sector (weighted mean, 3.5σ clipping, errors scaled by √χ²ᵣ), then sectors are
  combined with a 2″ systematic floor each (errors of the model are shared by every eclipse of a sector);
- verdict: *on target* if the target lies within 3σ and every Gaia star bright enough to fake the dip (a 100%
  eclipse at most) lies beyond 3σ; *on* a named neighbour if only it is within 3σ; *unresolved* if both are.
Saturated stars (T ≲ 6.8) bleed charge along CCD columns and their centroids move along the column; the test
does not apply to them (HIP 71639).

Calibration (`centroid_calib.py`, `eb_pixels.py --kostov`):
- 40 confirmed TESS planets (TOIs dispositioned CP or KP; T = 8–11, depth 0.07–1%, P = 2–40 d; up to four
  transits from two 2-minute sectors and two from a full-frame-only sector): 30 on their star, 5 unresolved (a
  star a few arcseconds away), 3 offset without a catalogued source, 2 untestable, none on a neighbour.
- 68 quadruple candidates of Kostov et al. (2022) with both periods below 10 d, whose authors traced both eclipse
  sets to one star: 36 of the 37 with usable eclipses have both sets on the target, 1 is unresolved, none is
  called a blend (31 had no sector where both eclipse sets were re-found).
So a verdict that a dip belongs to a named neighbour has not once been wrong on known on-target signals; an
"unresolved" or "offset" verdict is inconclusive, not evidence of a blend.

### 11.4 ZTF (`ztf_check.py`)
ZTF DR23 light curves are read from IRSA's public HATS/Parquet copy on AWS (IRSA's own services were down): the
position columns of the HEALPix partition, then only the needed rows. For each Track C system north of −28°, the
target and every Gaia star within 120″ able to make either eclipse depth: ZTF objects within 1.5″, points with
catflags = 0 and fainter than 12.5 mag (saturation), each object normalised to its median. For each eclipse set, a
box search over P(1 ± 5×10⁻⁴) with free phase (the TESS ephemeris drifts over ZTF's 2018–2024 span) and the same
search in 300 random period windows for the false-alarm probability. A star *eclipses* if SNR ≥ 7, FAP < 1% and
depth > 0.2× what it would need to be the source; it is *ruled out* if the depth + 3σ is below 0.3× what it would
need.

### 11.5 Gaia DR3 (`gaia_check.py`)
Each named star is matched to the Gaia DR3 source within 15″ whose position, moved back to J2000, is closest to
the TIC position and whose G fits T. From the ESA archive (TAP): RUWE, ipd_frac_multi_peak (a second star
resolved in Gaia's windows), radial-velocity χ² p-value and renormalised goodness of fit, non-single-star orbits
(`nss_two_body_orbit`), the eclipsing-binary model (`vari_eclipsing_binary`) and the variability class. For Track C
neighbours, parallax and proper-motion consistency with the target (co-moving if |Δϖ| < 3σ and Δv_t < 5 km/s).

### 11.6 KELT (`kelt_check.py`)
For the 28 bright single transits, KELT light curves within 15″ (TFA versions; archive times are JD for KELT-North
and JD + 2,450,000 for KELT-South) are searched at every period the TESS data allow, anchored on the TESS transit:
night-by-night means of the points within ±0.4 T14 of t0 − kP, weighted with a red-noise floor, computed with
difference arrays over ~10⁵–10⁶ trial periods. Calibration per star: the same search anchored on 40 false transit
times, and the TESS transit injected at three periods.

### 11.7 Vetting of the new candidates (`fit_new.py`)
Beyond §6 and §11.3: every transit the period predicts inside the data is measured; a joint fit (batman + emcee,
mixed cadences integrated over their exposure) with free stellar density, compared with the TIC density by a
second fit with a density prior; and, for any neighbour the centroids cannot exclude, a shape test: the transits
are fitted as a diluted eclipse with a large true radius ratio (k = 0.3–0.9) and compared with the planet fit.

## 12. Catalogue corrections and the literature (12 October)

### 12.1 Posted periods the hunt does not reproduce (`posted_check.py`, `posted_cases.py`)
Where a period has been posted to ExoFOP for the same candidate (same TIC and planet number, from the ExoFOP
planet-parameter uploads) and none of the hunt's allowed periods is within 0.2% of it, every transit that the
posted period predicts is measured (`period_review.measure`: box depth against a local quadratic baseline, with a
red-noise error), with the ephemeris through the candidate's catalogued transit and the period refined within
0.3% (the posted values are rounded, and the error grows over hundreds of cycles). The hunt's own best periods and
half and double the posted one are measured the same way. A predicted transit is *missing* when the data cover it
(coverage ≥ 50%) and show less than a third of the catalogued depth at more than 3σ below it; it is *seen* when its
SNR ≥ 3 and its depth is more than 0.4× the catalogued depth. Odd and even transits are compared, and phase 0.5 is
measured for an eclipsing-binary secondary. For the cases where the posted value fails, the competing ephemerides
(the TOI's listed one, the community TOI's own transits) are measured transit by transit, drawn, and tested at the
pixel level (`centroid_vet.vet`, up to four sectors).

### 12.2 Promoted community TOIs (`promoted_check.py`)
For every community TOI that ExoFOP marks as promoted to a TOI, the offset of its transit from that TOI's listed
ephemeris is computed from the two tables (12 October 2026), in hours and in units of the ephemeris error at that
epoch (period error × cycles, in quadrature with the epoch error). Those more than one transit duration and 5σ off
are checked against the light curves: is there a dip at the community TOI's epoch; does it match another TOI on
the same star (a label swap); are the transits that the TOI's ephemeris predicts present; and, if not, which period
within 0.5% of the TOI's, through the community TOI's transit, fits.

### 12.3 Fits of the corrected ephemerides (`corrections_fit.py`)
Joint transit fits (§7) of the corrected ephemerides. For a community TOI known only from transits separated by
a gap Δ, the period is Δ/n; for a fixed transit shape a/R★ grows in proportion to the assumed period, so the
stellar density a fit implies does too (ρ = 3π(a/R★)³/(GP²) ∝ P). One fit at a reference alias therefore gives the
implied density at every alias, and the aliases whose density matches the TIC value (within 2σ) are the likeliest
periods for a circular orbit.

### 12.4 Literature check of the new periods (`ads_periods.py`)
For each star behind a new period: all identifiers of the nearest SIMBAD object within 3″ of the TIC position
(HD, HIP, TYC, 2MASS, Gaia DR2/DR3, BD/CD/CPD, UCAC4, WISE, TOI, and others), plus the TIC and TOI/CTOI names, in
one NASA ADS full-text query; for stars with papers, a second query adds the period itself (the alias to one and
two decimals). Matches are read by eye: identifiers like "CD-73 14" also match unrelated text.
