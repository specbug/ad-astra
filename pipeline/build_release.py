"""Build the Zenodo data-and-code record for the TIC 402710280.01 preprint.

Everything is regenerated from the stored analysis products and light curves, so the record matches the paper:
  data/transit_windows/   each transit window (and the unused sector-66 window), cadence by cadence
  data/ffi_photometry/    our light curves of sectors 100-103 (no public FFI light curves exist for them)
  data/results/           transit times, parameters, predictions, alias and secondary-eclipse tests, blend test,
                          TRICERATOPS runs, difference-image centroids, sector-66 diagnostics
  code/                   the scripts that reproduce the paper (paper_modules(): every module the paper's scripts
                          import, transitively), at the commit given in README.md
  manuscript/             LaTeX source, bibliography, figures and PDF
plus LICENSE-DATA (CC-BY-4.0), LICENSE-CODE (MIT), README.md and SHA256SUMS.txt.

  python build_release.py
Output: release/tic402710280_zenodo/ and release/tic402710280_zenodo.zip
"""
import csv
import hashlib
import json
import math
import shutil
import subprocess
import zipfile

import env  # noqa: F401
import numpy as np

import transit_fit as tf
from env import CATALOGS, RESULTS

ROOT = RESULTS.parent
REL = ROOT / "release" / "tic402710280_zenodo"
TIC = 402710280
C = RESULTS / "candidates"


def write_csv(path, header, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        w.writerows(rows)


def alias_and_secondary(lcs, P, T0, D, depth_ppm, depth_err_ppm, known):
    from transit_fit import box_depth
    alias_rows, sec_rows = [], []
    use = [lc for lc in lcs if lc["meta"]["sector"] != 66]
    for k in range(2, 37):                   # every integer fraction of P down to 1.0 d
        Pa = P / k
        for lc in use:
            t = lc["time"]
            for n in range(math.ceil((t.min() - T0) / Pa), math.floor((t.max() - T0) / Pa) + 1):
                tc = T0 + n * Pa
                if min(abs(tc - x) for x in known) < 0.3:
                    continue
                b = box_depth(lc, tc, D)
                if b and b["coverage"] >= 0.6:
                    alias_rows.append([f"P/{k}", f"{Pa:.5f}", lc["meta"]["sector"], f"{tc + 2457000:.4f}",
                                       f"{b['depth'] * 1e6:.0f}", f"{b['err'] * 1e6:.0f}",
                                       f"{(depth_ppm - b['depth'] * 1e6) / math.hypot(b['err'] * 1e6, depth_err_ppm):.1f}"])
    for lc in use:
        t = lc["time"]
        for n in range(math.ceil((t.min() - T0 - P / 2) / P), math.floor((t.max() - T0 - P / 2) / P) + 1):
            tc = T0 + P / 2 + n * P
            b = box_depth(lc, tc, D)
            if b and b["coverage"] >= 0.6:
                sec_rows.append([lc["meta"]["sector"], f"{tc + 2457000:.4f}", f"{b['depth'] * 1e6:.0f}", f"{b['err'] * 1e6:.0f}"])
    return alias_rows, sec_rows


INTERMEDIATE = ("audit_402710280.json", "qlp_recheck_402710280.json", "qlp_timebase_402710280.json", "shape_free.json",
                "second_signal_402710280.json", "second_signal_chain.npy", "fpp_tic402710280.json",
                "fpp_tic402710280_v2.json", "fpp_tic402710280_v2_box5.json", "fpp_tic402710280_3.78d.json",
                "alias_test.csv", "secondary_eclipse.csv", "archive_comparison.csv",
                "prereg_predictions.json")

PAPER_SCRIPTS = ("ffi_new", "centroid_vet", "audit_402710280", "fit_new", "shape_chains", "timing_diagnostics_402710280",
                 "secondary_scan_402710280", "blend_shape_402710280", "fpp_triceratops", "second_signal_402710280",
                 "fpp_second_signal", "qlp_recheck_402710280", "qlp_timebase_402710280", "archive_comparison", "paper_figures",
                 "paper_numbers", "build_release", "build_prereg")


def paper_modules():
    """The paper's scripts and every pipeline module they import, directly or not (read from their source)."""
    import ast
    pipe = ROOT / "pipeline"
    local = {p.stem for p in pipe.glob("*.py")}
    seen, todo = set(), list(PAPER_SCRIPTS)
    while todo:
        m = todo.pop()
        if m in seen:
            continue
        seen.add(m)
        for n in ast.walk(ast.parse((pipe / f"{m}.py").read_text())):
            names = [a.name for a in n.names] if isinstance(n, ast.Import) else [n.module] if isinstance(n, ast.ImportFrom) and n.module else []
            todo += [x.split(".")[0] for x in names if x.split(".")[0] in local]
    return sorted(seen)


def main():
    from lcaccess import lcs_for
    if REL.exists():
        shutil.rmtree(REL)
    REL.mkdir(parents=True)
    fit = json.load(open(RESULTS / "report" / "transit_fits.json"))["TIC402710280"]
    ex = json.load(open(RESULTS / "report" / "extras.json"))["TIC402710280"]
    au = json.load(open(C / "audit_402710280.json"))
    lcs = [lc for lc in lcs_for(TIC) if lc["time"].size > 100]
    by = {lc["meta"]["sector"]: lc for lc in lcs}
    P, T0 = fit["period_d"][0], fit["tref_btjd"][0]
    times = {t["sector"]: t for t in fit["times"]}
    s66 = T0 + (20 - fit["ref_epoch"]) * P

    # transit windows: flux normalised by the out-of-transit quadratic, as used in the shape fit
    for s, tc, use in ((39, times[39]["tc"], "used"), (93, times[93]["tc"], "used"), (100, times[100]["tc"], "used"),
                       (102, times[102]["tc"], "used"), (66, s66, "not used (stray light)")):
        lc = by[s]
        w = tf.windows_from_lcs([lc], [tc], 0.8, 5.3)[0]
        write_csv(REL / "data" / "transit_windows" / f"tic402710280_s{s:03d}_window.csv",
                  ["bjd_tdb", "flux_normalised", "sector", "product", "exposure_s", "centre_bjd_tdb", "status"],
                  [[f"{x + 2457000:.6f}", f"{y:.7f}", s, lc["meta"].get("product"), f"{w[3]:.0f}", f"{tc + 2457000:.6f}", use]
                   for x, y in zip(w[0], w[1])])
    # our FFI photometry, sectors 100-103
    for s in (100, 101, 102, 103):
        z = np.load(CATALOGS.parent / "cache" / "ffi_lc" / f"{TIC}_s{s}.npz", allow_pickle=True)
        meta = json.loads(str(z["meta"]))
        rows = [[f"{t + 2457000:.6f}", f"{f:.7f}", f"{e:.7f}", f"{sa:.7f}", f"{b:.3f}", f"{x:.4f}", f"{y:.4f}"]
                for t, f, e, sa, b, x, y in zip(z["time"], z["flux"], z["flux_err"], z["sap_flux"], z["bkg"], z["cx"], z["cy"])]
        write_csv(REL / "data" / "ffi_photometry" / f"tic402710280_s{s:03d}_ffi.csv",
                  ["bjd_tdb", "flux", "flux_err", "sap_flux", "background", "centroid_col", "centroid_row"], rows)
        json.dump(meta, open(REL / "data" / "ffi_photometry" / f"tic402710280_s{s:03d}_meta.json", "w"), indent=1)
    res = REL / "data" / "results"
    rms = math.sqrt(np.mean([t["o_c_min"] ** 2 for t in fit["times"]]))
    write_csv(res / "transit_times.csv",
              ["sector", "light_curve", "epoch", "tc_bjd_tdb", "tc_err_d", "o_minus_c_min", "depth_ppm", "depth_err_ppm"],
              [[t["sector"], {"tspoc": "TESS-SPOC", "qlp": "QLP", "ffi": "this work"}[t["product"]], t["epoch"],
                f"{t['tc'] + 2457000:.5f}", f"{t['err_min'] / 1440:.5f}", f"{t['o_c_min']:.1f}",
                f"{t['depth_ppm']:.0f}", f"{t['depth_err_ppm']:.0f}"] for t in fit["times"]])
    params = {k: fit[k] for k in ("rp_rs", "b", "rho_solar", "log10_rho", "a_rs", "t14_h", "u", "n_points", "chi2_min",
                                  "chain_summary", "ephemeris_chi2", "depth_transits_ppm")}
    params.update(period_d=fit["period_d"], t0_bjd_tdb=[fit["tref_btjd"][0] + 2457000] + fit["tref_btjd"][1:],
                  timing_rms_min=rms, derived={k: ex[k] for k in ("a_au", "s_earth", "teq", "mass_assumed", "k_ms")},
                  notes="Values are [median, +1sigma, -1sigma]. depth_transits_ppm is [weighted mean, error] of the four "
                        "per-transit model depths at mid-transit. rho in solar units (1.408 g/cm^3).")
    json.dump(params, open(res / "parameters.json", "w"), indent=1, default=float)
    write_csv(res / "predicted_transits.csv",
              ["epoch", "tc_bjd_tdb", "sigma_ephemeris_min", "sigma_prediction_min", "utc_la_silla", "fraction_at_night_la_silla"],
              [[x["epoch"], f"{x['btjd'] + 2457000:.4f}", f"{x['sigma_min']:.1f}", f"{math.hypot(x['sigma_min'], rms):.1f}",
                x["utc"], f"{x['night_frac']:.2f}"] for x in ex["all_next"][:30]])
    D = fit["t14_h"][0] / 24
    al, sec = alias_and_secondary(lcs, P, T0, D, *fit["depth_transits_ppm"], [t["tc"] for t in fit["times"]])
    write_csv(res / "alias_test.csv", ["alias", "period_d", "sector", "predicted_bjd_tdb", "depth_ppm", "err_ppm", "exclusion_sigma"], al)
    write_csv(res / "secondary_eclipse.csv", ["sector", "phase_half_bjd_tdb", "depth_ppm", "err_ppm"], sec)
    for name in ("alias_test.csv", "secondary_eclipse.csv"):       # paper_numbers.py reads them from results/candidates
        shutil.copy(res / name, C / name)
    shutil.copy(C / "archive_comparison.csv", res / "archive_comparison.csv")
    for src, dst in (("blend_shape_402710280.json", "blend_test.json"), ("fpp_tic402710280_v2.json", "triceratops.json"),
                     ("fpp_tic402710280_v2_box5.json", "triceratops_5x5_apertures.json"),
                     ("shape_free.json", "shape_free_limb_darkening.json"),
                     ("timing_diagnostics_402710280.json", "timing_diagnostics.json"),
                     ("secondary_scan_402710280.json", "secondary_scan.json"),
                     ("second_signal_402710280.json", "signal_3.78d.json"),
                     ("fpp_tic402710280_3.78d.json", "triceratops_3.78d.json"),
                     ("paper_numbers.json", "paper_numbers.json"),
                     ("qlp_recheck_402710280.json", "qlp_comparison.json"),
                     ("qlp_timebase_402710280.json", "qlp_timebase.json")):
        if (C / src).exists():
            shutil.copy(C / src, res / dst)
    np.save(res / "shape_chain.npy", np.load(C / "shape_chain_fixed.npy"))
    json.dump({"combined": au["centroid_noS66"], "per_transit": [e for e in json.load(open(C / "new_candidates.json"))
                                                                  ["402710280"]["centroid"]["per_event"] if e.get("sector") != 66]},
              open(res / "centroids.json", "w"), indent=1, default=float)
    json.dump(au.get("s66_raw", {}), open(res / "sector66_diagnostics.json", "w"), indent=1, default=float)
    # code and manuscript
    commit = subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    (REL / "code" / "pipeline").mkdir(parents=True)
    mods = paper_modules()
    for m in mods:
        shutil.copy(ROOT / "pipeline" / f"{m}.py", REL / "code" / "pipeline" / f"{m}.py")
    shutil.copy(ROOT / "pipeline" / "prereg_protocol.template.tex", REL / "code" / "pipeline" / "prereg_protocol.template.tex")
    print(len(mods), "modules in code/pipeline:", " ".join(mods))
    shutil.copy(ROOT / "requirements.txt", REL / "code" / "requirements.txt")
    # the intermediate results the manuscript is filled from, so that paper_numbers.py runs from this record alone
    inter = REL / "code" / "results"
    (inter / "report").mkdir(parents=True)
    (inter / "candidates").mkdir()
    for name in ("transit_fits.json", "extras.json"):
        json.dump({"TIC402710280": json.load(open(RESULTS / "report" / name))["TIC402710280"]},
                  open(inter / "report" / name, "w"), indent=1, default=float)
    for name in INTERMEDIATE:
        shutil.copy(C / name, inter / "candidates" / name)
    (REL / "code" / "paper").mkdir()
    shutil.copy(ROOT / "paper" / "tic402710280.template.tex", REL / "code" / "paper" / "tic402710280.template.tex")
    man = REL / "manuscript"
    man.mkdir()
    for f in ("tic402710280.tex", "refs.bib", "tic402710280.pdf", "aastex701.cls", "aasjournalv7.bst"):
        shutil.copy(ROOT / "paper" / f, man / f)
    shutil.copytree(ROOT / "paper" / "figures", man / "figures")
    (REL / "LICENSE-DATA").write_text(
        "Data products, manuscript and figures in this record are licensed under the Creative Commons Attribution 4.0\n"
        "International licence (CC BY 4.0): https://creativecommons.org/licenses/by/4.0/\n")
    (REL / "LICENSE-CODE").write_text(
        "MIT License\n\nCopyright (c) 2026 Rishit Vora\n\nPermission is hereby granted, free of charge, to any person "
        "obtaining a copy\nof this software and associated documentation files (the \"Software\"), to deal\nin the Software "
        "without restriction, including without limitation the rights\nto use, copy, modify, merge, publish, distribute, "
        "sublicense, and/or sell\ncopies of the Software, and to permit persons to whom the Software is\nfurnished to do so, "
        "subject to the following conditions:\n\nThe above copyright notice and this permission notice shall be included in "
        "all\ncopies or substantial portions of the Software.\n\nTHE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF "
        "ANY KIND, EXPRESS OR\nIMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,\nFITNESS FOR A "
        "PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE\nAUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY "
        "CLAIM, DAMAGES OR OTHER\nLIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,\nOUT OF OR "
        "IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE\nSOFTWARE.\n")
    (REL / "README.md").write_text(README.format(commit=commit, P=P, T0=T0 + 2457000))
    # checksums and zip
    files = sorted(p for p in REL.rglob("*") if p.is_file())
    with open(REL / "SHA256SUMS.txt", "w") as fh:
        for p in files:
            fh.write(f"{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.relative_to(REL)}\n")
    zp = REL.parent / "tic402710280_zenodo.zip"
    with zipfile.ZipFile(zp, "w", zipfile.ZIP_DEFLATED) as z:
        for p in sorted(REL.rglob("*")):
            if p.is_file():
                z.write(p, p.relative_to(REL.parent))
    print(len(files), "files;", zp, f"{zp.stat().st_size / 1e6:.1f} MB")


README = """# TIC 402710280.01: data products and analysis code

Supporting record for the article *A 36.34-day Sub-Neptune Candidate Transiting the Nearby K Dwarf TIC 402710280*
(R. Vora, 2026), submitted to the Astronomical Journal.

Version 4. DOIs (all versions): this record 10.5281/zenodo.23390273; the preprint 10.5281/zenodo.23390316
(arXiv:2610.26056). Earlier versions: this record 10.5281/zenodo.23390274 (v1), 10.5281/zenodo.23466491 (v2) and
10.5281/zenodo.23538511 (v3); the preprint 10.5281/zenodo.23390317 (v1) and 10.5281/zenodo.23466634 (v2). The pre-registered test of the 2027
transits, 10.5281/zenodo.23390325, has one version and will not be revised before the 2027 data are analysed.

Ephemeris: P = {P:.5f} d, T0 = {T0:.4f} BJD_TDB. All times in this record are BJD_TDB.

## Contents

| Path | Description |
|---|---|
| `data/transit_windows/` | The four transits used (Sectors 39, 93, 100, 102) and the unused Sector 66 window, cadence by cadence, flux normalised by the out-of-transit quadratic baseline used in the shape fit. |
| `data/ffi_photometry/` | Our light curves of Sectors 100-103 from the calibrated TESS full-frame images (quality-flagged cadences removed; times barycentric for the target), with per-sector extraction metadata (camera, CCD, circular aperture radius in pixels around the target's position, scatter). QLP light curves of these sectors exist at MAST (released 2026 March-July) but are absent from the AWS copy of the archive that the pipeline reads; version 1 stated wrongly that none existed. |
| `data/results/transit_times.csv` | Individual mid-times, uncertainties (red-noise scaled), O-C and depths (Table 2 of the paper). |
| `data/results/parameters.json` | Ephemeris and transit-shape posterior summaries, convergence diagnostics and derived quantities (Table 3). |
| `data/results/predicted_transits.csv` | The next 30 predicted transits, with ephemeris-only and total (including the observed timing scatter) uncertainties, both including the covariance of T0 and P, and topocentric UTC and visibility at La Silla (Table 4). |
| `data/results/alias_test.csv` | Every transit predicted by the integer fractions P/2 to P/36 inside the data, with the measured depth and the exclusion significance: the observed depth minus the measured one, divided by the two uncertainties in quadrature. |
| `data/results/archive_comparison.csv` | The NASA Exoplanet Archive comparison sample of the introduction, with the query and its date in the header (`archive_comparison.py`). |
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
| `data/results/qlp_comparison.json` | Version 2: the Sector 100 and 102 transits timed in the QLP light curves, with the QLP time-stamp offsets measured against the full-frame image headers; our photometry of the Sector 39 and 93 transits as a control; and the Sector 66 event fitted with the final shape for three baselines (`qlp_recheck_402710280.py`). |
| `data/results/qlp_timebase.json` | Version 3: for each QLP light curve, the departure of TIME from a straight line in cadence number, the drift of the QLP offset, and the rate of change of the target's barycentric correction at the transit (`qlp_timebase_402710280.py`). |
| `code/results/`, `code/paper/` | The intermediate results the manuscript is filled from and the manuscript template, so that `paper_numbers.py` runs from this record alone. |
| `code/pipeline/` | The scripts that reproduce this paper and everything they import (Python), from git commit `{commit}` of the author's analysis repository. |
| `manuscript/` | LaTeX source (AASTeX 7.0.1 class and bibliography style included), bibliography, figures and PDF of the submitted article. |

## Reproducing the paper

**From this record alone.** With the packages in `code/requirements.txt`, run `python paper_numbers.py` in
`code/pipeline/`. It reads `code/results/` and writes `code/paper/tic402710280.tex`, the manuscript with the computed
values filled in, which matches `manuscript/tic402710280.tex`. No network access is needed.

**From the TESS data.** With the packages in `code/requirements.txt` (and TRICERATOPS 1.1 for the false-positive calculation), from
`code/pipeline/`:

1. `ffi_new.py` builds the Sector 100-103 light curves from the calibrated full-frame images on the MAST AWS bucket.
2. `audit_402710280.py` measures preliminary transit times, the Sector 66 diagnostics and the combined centroid.
3. `fit_new.py` fits the transit shape (each window centred on its measured mid-time), re-times the
   transits with that shape, and fits the ephemeris; `shape_chains.py free` repeats the shape fit with free limb
   darkening.
4. `timing_diagnostics_402710280.py`, `secondary_scan_402710280.py`, `blend_shape_402710280.py` and
   `fpp_triceratops.py` run the timing, secondary-eclipse, blend and TRICERATOPS tests;
   `second_signal_402710280.py` and `fpp_second_signal.py` characterise the 3.78-day signal.
5. `qlp_recheck_402710280.py` makes the comparison with the QLP light curves and the Sector 66 refit;
   `qlp_timebase_402710280.py` tests how the QLP time stamps were corrected to the barycenter.
6. `paper_numbers.py` fills the manuscript from the result files, `paper_figures.py` makes the figures, and
   `build_release.py` and `build_prereg.py` build this record and the pre-registration.

The per-transit difference-image centroids (`data/results/centroids.json`) come from `centroid_vet.py` (method in
`centroids.py`), run as `python centroid_vet.py 402710280:5.27:0.0021:39@2384.609,93@3838.185,100@4092.585,102@4128.945`.
The search that first found the two archival transits (a box-filter search of public full-frame light curves for
individual transit-like events) belongs to a wider survey and is not part of this record; the paper describes it.

Light curves are read from the public MAST archive (TESS-SPOC and QLP high-level science products, SPOC target
pixel files and calibrated full-frame images).

## Versions

- Version 4: the manuscript as submitted to the Astronomical Journal. The uncertainties of predicted transit times
  now include the covariance of T0 and P (correlation 0.47; 34 instead of 30 minutes at epoch 57). The alias test
  covers every integer fraction of the period down to 1.0 d (P/2 to P/36), and its exclusion significances include
  the uncertainty of the observed depth. Adds the Exoplanet Archive comparison sample (`archive_comparison.py`) and
  the intermediate results under `code/results/`, from which `paper_numbers.py` regenerates the manuscript without
  network access. The light curves, transit times, fits and ephemeris are unchanged.
- Version 3: adds `qlp_timebase_402710280.py` and its output, `data/results/qlp_timebase.json`. They show that QLP
  applied one barycentric correction to all cadences of Sectors 100 and 102 instead of computing it for each cadence,
  which explains the offset in the QLP time stamps and its drift. Also adds the manuscript submitted to the
  Astronomical Journal and the figures with embedded fonts. The light curves, transit times, fits, ephemeris and
  predictions are unchanged.
- Version 2 (2026 October 19): adds `data/results/qlp_comparison.json` and `qlp_recheck_402710280.py` (the QLP light
  curves of Sectors 100 and 102 reproduce the timing excess between those transits; their time stamps are 121 and
  139 s early at the transits), the Sector 66 refit, and the version 2 manuscript. The light curves, transit times,
  fits, ephemeris and predictions are unchanged.
- Version 1 (2026 October 15): first release.

## Licences

Data, manuscript and figures: CC BY 4.0 (`LICENSE-DATA`). Code: MIT (`LICENSE-CODE`).

## Integrity

`SHA256SUMS.txt` lists the SHA-256 checksum of every file.
"""


if __name__ == "__main__":
    main()
