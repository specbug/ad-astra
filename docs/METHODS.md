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
