# TIC 402710280.01: data products and analysis code

Supporting record for the preprint *TIC 402710280.01: A Sub-Neptune Candidate on a 36.34-Day Orbit around a Nearby
K Dwarf, Identified in TESS Full-Frame Images* (R. Vora, 2026).

Ephemeris: P = 36.33990 d, T0 = 2460838.2019 BJD_TDB. All times in this record are BJD_TDB.

## Contents

| Path | Description |
|---|---|
| `data/transit_windows/` | The four transits used (Sectors 39, 93, 100, 102) and the unused Sector 66 window, cadence by cadence, flux normalised by the out-of-transit quadratic baseline used in the shape fit. |
| `data/ffi_photometry/` | Our light curves of Sectors 100-103 from the calibrated TESS full-frame images, with per-sector extraction metadata (camera, CCD, aperture radius, scatter). No public full-frame light curves of these sectors existed when the analysis was done. |
| `data/results/transit_times.csv` | Individual mid-times, uncertainties (red-noise scaled), O-C and depths (Table 2 of the paper). |
| `data/results/parameters.json` | Ephemeris and transit-shape posterior summaries, convergence diagnostics and derived quantities (Table 3). |
| `data/results/predicted_transits.csv` | Predicted transits with ephemeris-only and total (including the observed 21.5-min timing scatter) uncertainties, and visibility from La Silla (Table 4). |
| `data/results/alias_test.csv` | Every transit predicted by P/2-P/5 inside the data, with the measured flux and the exclusion significance. |
| `data/results/secondary_eclipse.csv` | Flux at orbital phase 0.5. |
| `data/results/blend_test.json` | Diluted-eclipse fits for Gaia DR3 5828441289812582656 (Figure 4). |
| `data/results/triceratops*.json` | TRICERATOPS runs with the photometric apertures and with 5x5-pixel apertures. |
| `data/results/centroids.json` | Difference-image centroids, per transit and combined. |
| `data/results/sector66_diagnostics.json` | Quality-flag and coverage summary of the Sector 66 window. |
| `code/pipeline/` | The analysis pipeline (Python), git commit `a5e7bbbe2a631f771b6556cb966e42c39f6c1f73`. |
| `manuscript/` | LaTeX source, bibliography, figures and PDF of the preprint. |

## Reproducing the paper

With the packages in `code/requirements.txt` (and TRICERATOPS 1.1 for the false-positive calculation), from
`code/pipeline/`:

1. `ffi_new.py` builds the Sector 100-103 light curves from the calibrated full-frame images on the MAST AWS bucket.
2. `audit_402710280.py` measures the individual transit times and depths and the Sector 66 diagnostics.
3. `fit_new.py` fits the transit shape (each window centred on its measured mid-time) and the ephemeris.
4. `blend_shape_402710280.py` runs the diluted-eclipse test; `fpp_triceratops.py` the TRICERATOPS calculation.
5. `paper_figures.py` makes the figures; `build_release.py` builds this record.

Light curves are read from the public MAST archive (TESS-SPOC and QLP high-level science products, SPOC target
pixel files and calibrated full-frame images).

## Licences

Data, manuscript and figures: CC BY 4.0 (`LICENSE-DATA`). Code: MIT (`LICENSE-CODE`).

## Integrity

`SHA256SUMS.txt` lists the SHA-256 checksum of every file.
