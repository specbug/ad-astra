# pipeline

Every script is standalone (`python <script>.py`, run from this directory) and starts with a docstring giving its
inputs, method and outputs. `env.py` sets the paths: committed results go to `../results/`, bulk data to
`$ADASTRA_WORK`.

## Data access and catalogues
| Script | Role |
|---|---|
| `env.py` | Paths and process settings shared by every script |
| `s3io.py`, `s3parquet.py` | Light curves and parquet catalogues read directly from the MAST S3 bucket |
| `build_tic.py`, `build_manifest.py`, `build_spoc_manifest.py`, `ticdb.py` | TIC extract and per-sector light-curve manifests |
| `known.py`, `match_confirmed.py`, `prior.py` | Cross-matching with TOIs, community TOIs, confirmed planets and earlier checks |
| `ffi_new.py` | Photometry of the newest full-frame images (Sectors 99–108), with target-barycentric time stamps |

## Survey (Tracks B, C, D)
| Script | Role |
|---|---|
| `make_sample.py` | Survey sample and light-curve counts per sector |
| `lcprep.py`, `dips.py`, `survey.py` | Detrending and matched-filter dip search, on real and flipped light curves |
| `run_sector.py`, `worker.sh`, `survey_reduce.py`, `aggregate.py` | Per-sector driver, systematics reduction, cross-sector aggregation |
| `classify.py`, `stage2.py`, `vet_events.py`, `starcheck.py` | Candidate lists per track and stage-2 vetting of single dips |
| `finalize_b.py`, `finalize_d.py`, `vet_singles.py`, `duo.py`, `final_batch.sh` | Final Track B/D vetting and single/duo-transit period constraints |
| `ebsearch.py`, `eb_pixels.py`, `ztf_check.py` | Track C: eclipsing-binary subtraction, pixel tests, ZTF light curves |
| `completeness.py` | Blind recovery of known TOIs |

## Track A (small planets around nearby K/M dwarfs)
| Script | Role |
|---|---|
| `track_a_targets.py`, `periodic.py`, `run_track_a.py`, `worker_a.sh` | Target selection and multi-sector BLS search (real and flipped) |
| `triage_a.py`, `tracka_v2.py`, `tracka_vet.py` | Red-flag triage, comparison with known objects, pixel vetting |

## Vetting and follow-up of candidates
| Script | Role |
|---|---|
| `pixels.py`, `centroid_vet.py`, `centroid_calib.py` | Difference-image centroids and their calibration on confirmed planets |
| `vet_deep.py`, `dossier.py`, `plots.py` | Deep vetting and per-candidate dossiers |
| `gaia_check.py`, `kelt_check.py`, `ads_check.py`, `ads_periods.py`, `exofop_check.py` | Gaia DR3, KELT, NASA ADS and ExoFOP checks |
| `transit_fit.py`, `fit_headline.py`, `replot_fits.py` | Transit models (batman + emcee) of the headline candidates |

## Periods for period-less TOIs and the TOI ephemeris audit
| Script | Role |
|---|---|
| `period_hunt.py`, `period_review.py`, `period_recheck.py`, `period_signal_review.py`, `period_vet.py`, `period_deep.py` | New periods for single-transit TOIs and community TOIs |
| `posted_check.py`, `posted_cases.py`, `promoted_check.py`, `merged_ctoi_test.py`, `corrections_fit.py`, `s108_alias.py` | Posted periods and merged community TOIs the data disagree with |
| `toi_audit.py`, `toi_audit_review.py`, `toi_audit_timing.py` | Ephemeris audit of the TOI catalogue |
| `exofop_upload.py`, `upload_recheck.py` | Draft ExoFOP parameter upload and its live duplicate check |

## TIC 402710280.01 paper
| Script | Role |
|---|---|
| `audit_402710280.py` | Preliminary transit times, Sector 66 diagnostics, combined centroid |
| `fit_new.py` | Transit shape, re-timing, ephemeris and predictions (`--only-402710280`) |
| `shape_chains.py` | Shape fit with free limb darkening |
| `timing_diagnostics_402710280.py`, `secondary_scan_402710280.py`, `blend_shape_402710280.py`, `fpp_triceratops.py` | Timing, secondary-eclipse, blend and false-positive tests |
| `second_signal_402710280.py`, `fpp_second_signal.py` | The 3.78-day signal |
| `report_extras.py` | Derived quantities and transit predictions with site times |
| `paper_figures.py`, `paper_numbers.py` | Figures, and the manuscript filled from the result files |
| `build_release.py`, `build_prereg.py`, `prereg_protocol.template.tex` | Zenodo data record and pre-registration |

## Report
| Script | Role |
|---|---|
| `report.py` | The survey report page, `../results/report/index.html` |
