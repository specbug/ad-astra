# TESS QLP exploratory pilot — 11 October 2026

**No convincing new planet, previously unknown star, or novel stellar variability source has been established.** The completed pilot produced one strong repeated aperture modulation matching a known nearby RR Lyrae, plus a weaker unresolved signal. Contamination is the leading explanation for the stronger signal, with pixel attribution still inconclusive. These are useful vetting results, not discoveries.

## Completed observations and controls

The reconstructed run contains 24 stars, 48 QLP light curves from sectors 1 and 2, and 48 matching sign-reversed controls. Every requested light curve and control completed; no failed targets were dropped or replaced after flux inspection. Source hashes stayed unchanged throughout the run. Later science seasons remain unopened.

Selection began with a complete 1,155-row TIC field query. Fixed stellar cuts and dated NASA TOI/confirmed-planet and ExoFOP CTOI exclusions left a seed-ordered pool of 1,138. Three preserved metadata revisions added SIMBAD and publication checks; the final sample used the first 24 eligible targets among 47 seed positions. The [protocol](research-protocol.md) records the cuts and screening choices. The selection intentionally excludes some previously studied systems and therefore has catalog-dependent selection bias.

Known-signal controls recovered π Mensae c in sector 1 but missed it in sector 27, even among the three searched peaks. All 18 deeper injections were recovered in period and phase: periods 1/3/7 days, depths 1,000/3,000/10,000 ppm, using two products of the same bright star. This does not measure completeness on the pilot population, at shallow depths, or through upstream QLP processing. Strong peaks in inverted controls demonstrate that formal conditional signal-to-noise alone is inadequate.

The box search produced three direct cross-sector priorities and one harmonic priority. Independent visual and fixed-ephemeris checks promoted **none** as planet candidates. The [transit vetting](qlp-transit-vetting.md) records broad duty cycles, timing failures, weak fixed-ephemeris depths, and matching inverted signals. Two sign-reversed targets also produced direct cross-sector priorities; this is a diagnostic count, not a calibrated false-positive rate.

## TIC 149392671: likely known RR Lyrae contamination

The two strongest Lomb–Scargle periods are 0.476711 and 0.478290 days, with fitted semiamplitudes approximately 283 and 336 ppm. At the fixed sector-1 period, the signal is also present in SAP and alternative QLP apertures; removing any one-day interval leaves it present. The target changes from camera 4/CCD 4 to camera 4/CCD 1. It is the strongest of these 24 stars at that fixed frequency in both sectors. These checks support a repeated signal in the aperture, but do not locate it on the named target.

Wider neighbor queries identify **OGLE LMC-RRLYR-34351 / Gaia DR3 4757558810204649856**, 29.25 arcseconds from the target, or about 1.39 TESS pixels. Gaia's RR Lyrae table gives a fundamental period **0.477113596 days**, and VSX also lists 0.477114 days. VSX may incorporate Gaia information; this agreement is a catalog cross-reference, not an independent photometric experiment. Its approximately 1.20-mag peak-to-peak G variation can plausibly supply the much smaller blended signal. This is a strong, specific contamination explanation.

At the published neighbor period, fixed-frequency fits give 283 and 330 ppm fundamental semiamplitudes. The fitted phases differ by only about 0.021 cycles across sectors; substantial second and third harmonics are also present. No optimization of the published period was performed. Phase agreement and harmonics strengthen the known-RR explanation; they do not by themselves constitute pixel localization.

Additional 7×7 TESSCut full-frame cutouts from sectors 1 and 2 passed original-byte FITS checksum validation. Fixed published-period pixel maps do **not** reliably localize the signal: the strongest pixels move with sector and background treatment. The pixel nearest the RR Lyrae shows repeated modulation (0.715 ± 0.136 and 0.337 ± 0.041 electrons/second in one background model), but these errors omit temporal and interpixel covariance. Strong target-centered power in sector 1 also depends on background treatment. Attribution to the RR Lyrae remains a likely explanation, not a demonstrated pixel-level identification; the [pixel report](qlp-rrlyrae-pixels.md) retains both models. No later-sector pixels were opened.

The target is also the bright member of a **known common-proper-motion white-dwarf pair**, reported by [Raddi et al. 2022](https://cdsarc.cds.unistra.fr/viz-bin/cat/J/A+A/658/A22). Its faint white-dwarf companion is a different source, Gaia DR3 4757558814501499392, about 6.20 arcseconds away. Do not conflate this companion with the RR Lyrae at 29.25 arcseconds, or describe the primary as the white dwarf.

## TIC 150104580: weak unresolved structure, not a discovery candidate

The initial sector periods are 1.193842 and 1.205756 days, with semiamplitudes approximately 107 and 101 ppm. Holding the sector-1 period fixed gives approximately 105 and 92 ppm in KSPSAP, with SAP counterparts around 92 and 82 ppm. Robust regression and leave-one-day-out checks preserve the low-amplitude structure. However, its fitted phase changes by approximately **0.414 cycles**, and other pilot stars have comparable or larger amplitudes at that frequency. The evidence is insufficient for a stable, localized stellar signal.

A known RR Lyrae, **OGLE LMC-RRLYR-37553 / Gaia DR3 5476797641361944064**, lies about 46 arcseconds away and has a published period of 0.544290167 days. Its frequency minus one cycle per day corresponds to **1.194378807 days**, close to the apparent longer period. This arithmetic coincidence is a warning, not proof of an alias: with the actual TESS sampling and weights, a pure sinusoid at the RR period leaks only about 4% of its amplitude into that daily alias under our linear model. A joint fit retains both components. Upstream processing, additional contaminants and correlated noise remain untested explanations.

The nearby RR period itself fits approximately 96 and 82 ppm in the aperture. A common offset-and-slope plus RR-harmonics and daily-alias model gives alias amplitudes around 103 and 95 ppm, but the alias phase differs by roughly 0.408 cycles. No significance calibrated for selection or red noise is assigned.

Subsequent 12 October prior-art checks found a published TARS sector-2 header peak at 1.1989171 days. This approximately 1.2-day periodicity therefore has an existing public counterpart; its physical source and interpretation remain unresolved. See the [rotation review](rotation-prior-review.md).

GALAH's low rotational-broadening estimates add context: DR4 reports roughly 4.57 ± 1.59 km/s for the exact target identity. Combining the TIC radius with an assumed 1.2-day rotation implies an equatorial speed near 50 km/s, requiring a small inclination if broadening traces projected rotation. This does not rule rotation out, and the broadening model also includes nonrotational contributions.

## Scope of prior-art checks

The final 24 targets have complete ASAS-SN, VSX and GCVS responses with no returned rows within 10 arcseconds, and Gaia DR3 variability/rotation queries within 5 arcseconds with no reported classifications or periods. Those small cones missed the crucial variables farther out in TESS's large pixels. Missing classifications mean **unchecked or unclassified**, not constant stars.

Additional exact-identity publication checks excluded a previously published CHESPA transit candidate, known/potential binaries, a published rotation source and earlier survey artifacts before science flux was opened. The complete CHESPA response has 222 rows, including likely binaries; 116 in the paper title is not the returned row count. Avallone study-membership exclusions preserve each original class and do not assert every excluded star has a measured period or confirmed binary.

An initial GCVS query used an incorrect table name and therefore did not check that catalog. Its source and results are preserved; the corrected query explicitly verifies that all three requested catalog tables are returned. Harmless, understood legacy VOTable ID/unit warnings are recorded; unexpected warnings fail the parse. See `data/qlp-prior-review/review-corrections.json` and the final review outputs. ADS was omitted during this pilot's original review. Following the user's later authorization, scoped ADS checks completed on 12 October; see the [ADS review](ads-literature-review.md). No exhaustive literature or novelty clearance is claimed.

## Reproducibility and limitations

Core artifacts:

- `data/qlp-frozen/manifest.json` and `config.json`: exact reviewed targets, product URIs, protected later products, input/source hashes and limitations.
- `runs/qlp-pilot/results.json` and `source-snapshot/`: all 24 target results, 48 original/control analyses and original screened source.
- `runs/qlp-pilot-diagnostics/summary.csv`, target dossiers, PNGs and PDF contact sheets: all-target inspection, not only successes.
- `runs/qlp-variability-followup-v3/results.json`: fixed sector-1 frequencies, SAP/aperture/centroid/background checks, robust fits, daily omissions and all-target comparison amplitudes.
- `runs/qlp-catalog-period-tests/`: published neighbor periods, harmonic/alias fits, sampling leakage and figures.
- `data/qlp-stellar-metadata/` and `data/qlp-prior-review-v2/`: public catalog responses and query provenance.

The first variability follow-up left only its declaration after its execution session became unavailable. A recorded repeat then exposed a rank-deficient nuisance model for isolated one-sample segments. A corrected version retains each such sample with an offset and omits its unidentifiable slope; regression tests verify signal recovery and robust fitting. The failed runs, exact previous source versions and recovery decisions remain preserved. The already completed pilot screen was not changed or rerun.

All subsequent diagnostics use data already examined during screening. They are exploratory, with no independent held-out confirmation, calibrated survey false-alarm probability or transit completeness estimate. QLP correction and our added detrending can change shapes; alternative apertures are not independent instruments. Photometric ppm values refer to the blended aperture and fitted semiamplitude unless otherwise specified.

The original pre-reset 88-star work cannot currently be reproduced from this machine. Its statements remain historical conversation evidence only; see [recovery notes](recovery-notes.md). New science retrieval is restricted literally to sectors 1/2 to protect every potential old holdout. A saved environment configuration or local code commit is not a backup of ignored raw data; source and science artifacts require a separate verified checkpoint.
