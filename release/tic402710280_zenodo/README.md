# TIC 402710280.01: data products and analysis code

Supporting record for the preprint *A 36.34-day Sub-Neptune Candidate Transiting the Nearby K Dwarf TIC 402710280,
Identified in TESS Full-Frame Images* (R. Vora, 2026).

Ephemeris: P = 36.33986 d, T0 = 2460838.2008 BJD_TDB. All times in this record are BJD_TDB.

## Contents

| Path | Description |
|---|---|
| `data/transit_windows/` | The four transits used (Sectors 39, 93, 100, 102) and the unused Sector 66 window, cadence by cadence, flux normalised by the out-of-transit quadratic baseline used in the shape fit. |
| `data/ffi_photometry/` | Our light curves of Sectors 100-103 from the calibrated TESS full-frame images (quality-flagged cadences removed; times barycentric for the target), with per-sector extraction metadata (camera, CCD, circular aperture radius in pixels around the target's position, scatter). No public full-frame light curves of these sectors existed when the analysis was done. |
| `data/results/transit_times.csv` | Individual mid-times, uncertainties (red-noise scaled), O-C and depths (Table 2 of the paper). |
| `data/results/parameters.json` | Ephemeris and transit-shape posterior summaries, convergence diagnostics and derived quantities (Table 3). |
| `data/results/predicted_transits.csv` | The next 30 predicted transits, with ephemeris-only and total (including the observed timing scatter) uncertainties, and topocentric UTC and visibility at La Silla (Table 4). |
| `data/results/alias_test.csv` | Every transit predicted by P/2-P/5 inside the data, with the measured flux and the exclusion significance. |
| `data/results/secondary_eclipse.csv` | Flux at orbital phase 0.5. |
| `data/results/blend_test.json` | Diluted-eclipse fits for Gaia DR3 5828441289812582656 (Figure 4). |
| `data/results/triceratops*.json` | TRICERATOPS runs with the photometric apertures and with 5x5-pixel apertures. |
| `data/results/centroids.json` | Difference-image centroids, per transit and combined (Figure 5). |
| `data/results/shape_chain.npy` | Posterior samples of the shape fit (columns Rp/R*, b, log10 rho*/rho_sun; float32). |
| `data/results/shape_free_limb_darkening.json` | The shape fit with free limb darkening (Kipping 2013 q1, q2). |
| `data/results/timing_diagnostics.json` | Mid-times under varied baselines and windows, the Sector 93 flux-series test, variability periodograms, Sector 103 coverage. |
| `data/results/secondary_scan.json` | Secondary-eclipse search over orbital phases 0.3-0.7. |
| `data/results/signal_3.78d.json`, `triceratops_3.78d.json` | The 3.78-day signal: searches, per-sector depths, transit fit, and TRICERATOPS (Figure 6). |
| `data/results/paper_numbers.json` | Every number filled into the manuscript template by `paper_numbers.py`. |
| `data/results/sector66_diagnostics.json` | Quality-flag and coverage summary of the Sector 66 window. |
| `code/pipeline/` | The analysis pipeline (Python), git commit `39e7e4e98d253ecf08f9dcd5c39ec57c55bfef34`. |
| `manuscript/` | LaTeX source, bibliography, figures and PDF of the preprint. |

## Reproducing the paper

With the packages in `code/requirements.txt` (and TRICERATOPS 1.1 for the false-positive calculation), from
`code/pipeline/`:

1. `ffi_new.py` builds the Sector 100-103 light curves from the calibrated full-frame images on the MAST AWS bucket.
2. `audit_402710280.py` measures preliminary transit times, the Sector 66 diagnostics and the combined centroid.
3. `fit_new.py --only-402710280` fits the transit shape (each window centred on its measured mid-time), re-times the
   transits with that shape, and fits the ephemeris; `shape_chains.py free` repeats the shape fit with free limb
   darkening.
4. `timing_diagnostics_402710280.py`, `secondary_scan_402710280.py`, `blend_shape_402710280.py` and
   `fpp_triceratops.py` run the timing, secondary-eclipse, blend and TRICERATOPS tests;
   `second_signal_402710280.py` and `fpp_second_signal.py` characterise the 3.78-day signal.
5. `paper_numbers.py` fills the manuscript from the result files, `paper_figures.py` makes the figures, and
   `build_release.py` builds this record.

Light curves are read from the public MAST archive (TESS-SPOC and QLP high-level science products, SPOC target
pixel files and calibrated full-frame images).

## Licences

Data, manuscript and figures: CC BY 4.0 (`LICENSE-DATA`). Code: MIT (`LICENSE-CODE`).

## Integrity

`SHA256SUMS.txt` lists the SHA-256 checksum of every file.
