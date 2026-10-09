# ad-astra: an amateur search of NASA TESS data for things nobody has reported yet

*Ad astra per aspera.* A from-scratch survey of public TESS light curves for new transiting planets,
eclipsing multiple-star systems and unusual dimming events. Status: **in progress**. This page is
updated as results arrive. Nothing here is a confirmed discovery.

## Search tracks
| Track | Looking for | Why it might still hide in the data |
|---|---|---|
| A | Small planets (≈1–2 Earth radii) around nearby K/M dwarfs | Stars with only full-frame data, or with at most two 2-minute sectors, never had all their sectors searched together |
| B | Long-period giant planets seen as single or "duo" transits | Standard searches need ≥ 2–3 transits per search window; QLP's own detrending flattens dips longer than about half a day |
| C | Eclipsing triples and 2+2 quadruples | Extra eclipse sets are found by subtracting the main binary's eclipses, which few pipelines do systematically |
| D | Anomalies: deep aperiodic dips, exocomet-like asymmetric dips | No standard pipeline looks for them |

## Data (all public, read straight from the MAST open-data bucket on AWS, `s3://stpubdata`)
- TESS SPOC 2-minute light curves, sectors 1–107
- QLP (MIT Quick-Look Pipeline) full-frame light curves, sectors 1–98; TESS-SPOC full-frame light curves, sectors 1–81
- Full-frame image cubes, for pixel-level checks (via `astrocut`)
- TESS Input Catalog v8.2 and Gaia DR3, from the HATS parquet copies in the same bucket
- Known-object lists: ExoFOP TOI catalogue (Feb 2026 export) and a NASA TOI snapshot (Oct 2026)

## Pipeline (`pipeline/`)
| Module | Role |
|---|---|
| `s3io.py`, `s3parquet.py` | Light-curve and catalogue access on S3 (HTTP range reads of only the needed parquet columns) |
| `lcprep.py` | Binning and a fast Tukey-biweight detrender (numba) |
| `dips.py`, `survey.py` | Matched-filter dip finder with empirically measured (red-noise-aware) noise; survey driver |
| `survey_reduce.py` | Spacecraft and CCD-level systematics, blends between neighbours, quality flags |
| `classify.py` | Turns dips into Track B/C/D candidate lists |
| `ebsearch.py` | Eclipsing-binary template subtraction and residual search (Track C) |
| `periodic.py`, `run_track_a.py`, `triage_a.py` | Multi-sector BLS search and red-flag triage (Track A) |
| `duo.py` | Period aliases of duo-transits |
| `run_sector.py`, `worker.sh` | Per-sector driver used by parallel worker sessions |

Every search is also run on light curves flipped upside-down (dips become bumps). Anything the
pipeline "finds" in flipped data is a false alarm by construction, which gives a measured false-alarm
rate instead of an assumed one.

## Reference and prior art
Inspired by [paraschopra/tess-planet-candidates](https://github.com/paraschopra/tess-planet-candidates),
which searched the 526 stars within 60 pc that never got 2-minute data. This project deliberately
covers different stars and different kinds of signals.
