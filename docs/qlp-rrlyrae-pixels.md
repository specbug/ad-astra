# TIC 149392671: known RR Lyrae neighbour and pixel check

The approximately 0.477-day aperture signal has a plausible known contaminant:
Gaia DR3 4757558810204649856, also **OGLE LMC-RRLYR-34351**, is an RRab variable
29.2515 arcseconds from TIC 149392671. Its published Gaia fundamental period is
**0.47711359612012993 days**. This makes the aperture signal unsuitable for a
claim of newly discovered stellar variability without resolving its source.
The TESS pixel analysis supports keeping this contamination explanation open,
but **does not establish a unique source location**.

## What was measured

On 2026-10-11 UTC, two 7×7 TESSCut products were retrieved for sectors 1 and 2
only. Sector 1 was observed on camera 4, CCD 4; sector 2 on camera 4, CCD 1.
There were 1,267 and 1,196 quality-zero cadences respectively. No later sector
or protected holdout was downloaded.

The fixed-period diagnostic was declared before the pixel data were read.
It fits three harmonics of the published Gaia period, with an intercept and
linear baseline for each continuous observing segment (two per sector).
Two models are retained: those baselines alone, and the same model with the
median outer-border flux as an additional background nuisance regressor.
There was no period optimization, clipping, optimized aperture selection, or
formal centroid localization. Positions use Gaia coordinates propagated from
epoch 2016.0 to each sector with the catalogued proper motions.

For the model including the border background, the fundamental semiamplitudes
at the pixels nearest the two sources were:

| Sector | Pixel nearest bright TIC target (e−/s) | Pixel nearest known RRab (e−/s) |
| --- | ---: | ---: |
| 1 | 6.603 ± 2.278 | 0.715 ± 0.136 |
| 2 | 0.149 ± 0.356 | 0.337 ± 0.041 |

These uncertainties assume independent residuals and are only formal
diagnostics. Temporal and interpixel correlations, the background estimate,
and the point-spread function are not incorporated in them.

The sector 1 largest-amplitude pixel is close to the bright target and its
0.30-pixel companion; those two sources are unresolved. The sector 2 maximum
is between the target and the RRab neighbour. The nearest-RRab pixel has
repeat modulation at the fixed period, but changing the background treatment
changes the amplitudes considerably. The signed phase-projection plots use
a fixed 3×3 aperture around the RRab's position; that aperture also receives
flux from the bright star and is not an independent RRab phase template.

The repeated period near an already catalogued RRab, coupled with the
unresolved pixel evidence, is sufficient to withhold a novelty claim. It is
not sufficient to claim the bright target is constant, to identify a new
variable companion, or to establish that every part of the aperture signal
comes from the RRab. No PRF deblending, injection localization calibration,
or new held-out observations were used.

## Provenance and reproduction

- `data/qlp-stellar-metadata/neighbor-rrlyrae.json` holds the Gaia variable-star
  period record; `149392671/simbad-cone.json` gives the exact OGLE/Gaia alias
  match.
- `data/qlp-rrlyrae-pixels/s0001` and `s0002` preserve each original ZIP, original
  extracted FITS, request parameters, UTC timestamps, and SHA-256 hashes.
- All three HDUs in both FITS products supplied CHECKSUM and DATASUM cards.
  They passed the repository's original-byte verification before analysis.
- The first attempts failed because the generic S3 delivery hostname was
  blocked. Those failures remain in
  `runs/qlp-rrlyrae-pixels/first-download-failures.json`. The public bucket's
  already allowed official regional endpoint delivered the same generated
  object keys successfully. No authentication or unsigned-query rewriting
  was involved.
- `runs/qlp-rrlyrae-pixels/declaration.json`, `results.json`, and four PNG maps
  preserve the declared plan, final source hash, complete pixel coefficients,
  positions, and both background treatments. The declaration's first-fetch
  script hash intentionally predates completing the analysis implementation;
  this was a post-selection exploratory diagnostic throughout.

Run from the repository root:

```sh
.venv/bin/python scripts/qlp_rrlyrae_pixels.py
.venv/bin/pytest -q tests/test_rrlyrae_pixels.py
```

The tests verify recovery of known harmonic coefficients with segment offsets,
slopes, and a varying background; the literal sector-1/2 network guard; and
rejection of pixels without positive uncertainties.
