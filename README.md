# ad-astra

A search of public NASA TESS data for transiting planets, eclipsing multiple stars and unusual dimming events that
had not been reported, run end to end from the MAST open-data archive. Its principal result is the planet candidate
described below. Nothing in this repository is a confirmed planet.

## TIC 402710280.01: a 36.34-day sub-Neptune candidate

> R. Vora (2026), *A 36.34-day Sub-Neptune Candidate Transiting the Nearby K Dwarf TIC 402710280*. Preprint: [doi:10.5281/zenodo.23390317](https://doi.org/10.5281/zenodo.23390317) ([`paper/tic402710280.pdf`](paper/tic402710280.pdf)).
> Data and code: [doi:10.5281/zenodo.23390274](https://doi.org/10.5281/zenodo.23390274). Pre-registered test of the 2027 transits: [doi:10.5281/zenodo.23390325](https://doi.org/10.5281/zenodo.23390325).

| | |
|---|---|
| Host | TIC 402710280 (TYC 9041-03677-1), K dwarf, T = 9.24, 66.8 pc |
| Period | 36.33986 ± 0.00051 d, from four transits (TESS Sectors 39, 93, 100, 102) |
| Radius | 3.75 ± 0.22 R⊕ (R<sub>p</sub>/R<sub>★</sub> = 0.0421) |
| Equilibrium temperature | 498 ± 21 K |
| False-positive probability | 0.003 ± 0.004 (TRICERATOPS); no high-resolution imaging yet |
| Next transits | 2027 Feb 19 (TESS Sector 114); 2027 May 3, 03:41 UTC (full night at La Silla) |

The preprint also reports a 407 ± 26 ppm signal with a period of 3.78 days whose source (the target or a background
star 3.5″ away) the present data cannot identify.

| Path | Content |
|---|---|
| [`paper/`](paper) | Manuscript (AASTeX): `tic402710280.template.tex` is the source; `pipeline/paper_numbers.py` fills every number from the result files into `tic402710280.tex`. Figures in `paper/figures/`. |
| [`release/tic402710280_zenodo.zip`](release) | The data-and-code record for Zenodo: light curves of every transit window, our photometry of Sectors 100–103, transit times, posterior summaries and chains, all vetting results, the pipeline snapshot and the manuscript. |
| [`release/prereg_tic402710280/`](release/prereg_tic402710280) | Pre-registered test of the 2027 transits: protocol, frozen test code, predictions and validation. |
| [`release/ZENODO_UPLOAD.md`](release/ZENODO_UPLOAD.md) | Upload order and metadata for the three Zenodo records. |

Reproducing the paper from `pipeline/` (after the Sector 100–103 light curves exist, `ffi_new.py`):

```
python audit_402710280.py                 # preliminary times, Sector 66 diagnostics, combined centroid
python fit_new.py                         # transit shape, re-timing, ephemeris, predictions
python shape_chains.py free               # free limb darkening
python timing_diagnostics_402710280.py    # baseline/window/flux-series timing tests, variability
python secondary_scan_402710280.py        # secondary eclipse over phases 0.3-0.7
python blend_shape_402710280.py           # diluted-eclipse test of the 3.5" neighbour
python fpp_triceratops.py 10              # TRICERATOPS
python second_signal_402710280.py         # the 3.78-day signal; then fpp_second_signal.py
python paper_figures.py && python paper_numbers.py && (cd ../paper && latexmk -pdf tic402710280.tex)
python build_release.py && python build_prereg.py
```

## The survey

| Track | Looking for | Why it might still be hidden |
|---|---|---|
| A | Small planets (≈1–2 R⊕) around nearby K and M dwarfs | Stars with only full-frame data, or few two-minute sectors, were never searched across all their sectors |
| B | Long-period planets seen as single or "duo" transits | Standard searches need several transits per search window |
| C | Eclipsing triples and 2+2 quadruples | Extra eclipse sets appear only after the main binary is subtracted |
| D | Deep aperiodic or asymmetric dips | No standard pipeline looks for them |

The search covered 10,324,827 public light curves (QLP full-frame light curves of Sectors 1–98 and SPOC two-minute
light curves of nine sectors), plus our own photometry of the 2026 full-frame images (Sectors 99–108). Every search
is repeated on light curves flipped upside down, so that the false-alarm rate is measured rather than assumed. The
methods are described in [`docs/METHODS.md`](docs/METHODS.md); the survey report is
[`results/report/index.html`](results/report/index.html).

### Data

All inputs are public and are read directly from the MAST open-data bucket on AWS (`s3://stpubdata`): SPOC
two-minute, TESS-SPOC and QLP light curves, the calibrated full-frame image cubes (via `astrocut`), the TESS Input
Catalog v8.2 and Gaia DR3. Known objects come from ExoFOP (TOI and community-TOI tables), the NASA Exoplanet Archive,
SIMBAD, VizieR and NASA ADS.

### Layout

| Path | Content |
|---|---|
| [`pipeline/`](pipeline) | All analysis code; [`pipeline/README.md`](pipeline/README.md) groups the scripts by role. |
| [`results/`](results) | Compact, committed outputs: per-sector survey tables, candidate lists, vetting results, dossiers and the report. |
| [`docs/`](docs) | Methods, and the brief used to run the survey in parallel worker sessions. |

Bulk data (light curves, catalogue extracts, caches) are kept outside the repository, in `$ADASTRA_WORK`
(default: `../work` next to the repository).

## Setup

```
python3.13 -m venv .venv && .venv/bin/pip install -r requirements.txt
export ADASTRA_WORK=/path/to/scratch   # optional
```

## Licences and citation

Code: MIT ([`LICENSE`](LICENSE)). Data products, figures and text: CC BY 4.0 ([`LICENSE-DATA`](LICENSE-DATA)).
To cite this work, see [`CITATION.cff`](CITATION.cff).

The analysis code and the manuscript were developed with the assistance of an AI system (Claude, Anthropic); the
author directed the work and is responsible for its content.
