"""Build the pre-registration record for TIC 402710280.01 (release/prereg_tic402710280/) from the paper's fit.

  1. writes the registered quantities (ephemeris, timing scatter, transit shape, depth window) into the frozen test,
     prereg_test.py, between its two marker lines;
  2. writes predictions.json (epochs 57-66);
  3. validates the test on the four archival transits of the data record (each must come out "recovered"), on
     transit-free Sector 101 photometry shifted onto a predicted mid-time ("excluded"), and on the same data with
     the prediction inside a gap ("insufficient data"); results in validation.json;
  4. fills PROTOCOL.tex from prereg_protocol.template.tex and compiles PROTOCOL.pdf;
  5. writes SHA256SUMS.txt and release/prereg_tic402710280.zip.

  python build_prereg.py      (after fit_new.py and build_release.py)
"""
import csv
import hashlib
import importlib.util
import json
import math
import re
import shutil
import subprocess
import tempfile
import zipfile
from pathlib import Path

import env  # noqa: F401
import numpy as np

from env import RESULTS

ROOT = RESULTS.parent
REC = ROOT / "release" / "prereg_tic402710280"
DATA = ROOT / "release" / "tic402710280_zenodo" / "data"
TPL = Path(__file__).with_name("prereg_protocol.template.tex")
BEGIN = "# ---- registered quantities (from the preprint, Table 3; do not edit) "
END = "WINDOW = 0.8                     # d, half-width of the analysis window around the prediction\n"


def registered(f):
    oc = [t["o_c_min"] for t in f["times"]]
    rms = math.sqrt(np.mean(np.square(oc)))
    dep = f["depth_transits_ppm"][0]
    return {"P": f["period_d"][0], "SIG_P": f["period_d"][1], "T0": f["tref_btjd"][0] + 2457000,
            "SIG_T0": f["tref_btjd"][1], "RMS": rms, "REF": f["ref_epoch"], "K": f["rp_rs"][0], "B": f["b"][0],
            "ARS": f["a_rs"][0], "U": tuple(f["u"]), "T14": f["t14_h"][0], "DEPTH": dep}


def write_test(r):
    p = REC / "prereg_test.py"
    s = p.read_text()
    a = s.index(BEGIN)
    b = s.index(END)
    block = (BEGIN + "-" * 10 + "\n"
             f"P = {r['P']!r}            # d\n"
             f"T0 = {r['T0']!r}        # BJD_TDB, epoch {r['REF']} (the Sector 39 transit is epoch 0)\n"
             f"SIG_P = {r['SIG_P']!r}    # d\n"
             f"SIG_T0 = {r['SIG_T0']!r}    # d\n"
             f"SIG_RMS = {r['RMS']:.2f} / 1440           # d, rms of the four timing residuals\n"
             f"REF_EPOCH = {r['REF']}\n"
             f"RP_RS, B, A_RS = {r['K']!r}, {r['B']!r}, {r['ARS']!r}\n"
             f"U = ({r['U'][0]:.5f}, {r['U'][1]:.5f})\n"
             f"T14 = {r['T14']!r} / 24    # d\n"
             f"DEPTH = {r['DEPTH']:.1f}                   # ppm, expected depth at mid-transit\n"
             "DEPTH_LO, DEPTH_HI = 0.5 * DEPTH, 1.5 * DEPTH\n"
             "DCHI2_DETECT = 25.0\n")
    p.write_text(s[:a] + block + s[b:])


def load_test():
    spec = importlib.util.spec_from_file_location("prereg_test", REC / "prereg_test.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def csv_lc(path, t, f):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["bjd_tdb", "flux"])
        w.writerows(zip(t, f))


def validate(m):
    out = {"archival": [], "null": []}
    tmp = Path(tempfile.mkdtemp())
    for s, ep, exp in ((39, 0, 600.0), (93, 40, 200.0), (100, 47, 200.0), (102, 48, 200.0)):
        rows = list(csv.DictReader(open(DATA / "transit_windows" / f"tic402710280_s{s:03d}_window.csv")))
        t = np.array([float(x["bjd_tdb"]) for x in rows])
        f = np.array([float(x["flux_normalised"]) for x in rows])
        csv_lc(tmp / "lc.csv", t, f)
        r = m.main(str(tmp / "lc.csv"), ep, exp)
        r["sector"] = s
        out["archival"].append(r)
    rows = list(csv.DictReader(open(DATA / "ffi_photometry" / "tic402710280_s101_ffi.csv")))
    t = np.array([float(x["bjd_tdb"]) for x in rows])
    f = np.array([float(x["flux"]) for x in rows])
    tc57, _ = m.predicted(57)
    for label, centre in (("transit-free, offset A", 2461106.0), ("transit-free, offset B", 2461115.0)):
        csv_lc(tmp / "lc.csv", t - centre + tc57, f)
        r = m.main(str(tmp / "lc.csv"), 57, 200.0)
        r["case"] = label
        out["null"].append(r)
    gaps = np.where(np.diff(t) > 0.5)[0]
    g = gaps[0]
    mid = 0.5 * (t[g] + t[g + 1])
    csv_lc(tmp / "lc.csv", t - mid + tc57, f)
    r = m.main(str(tmp / "lc.csv"), 57, 200.0)
    r["case"] = "prediction inside a data gap"
    out["null"].append(r)
    shutil.rmtree(tmp)
    for r in out["archival"]:
        assert r["outcome"] == "recovered", r
    assert [r["outcome"] for r in out["null"]] == ["excluded", "excluded", "insufficient data"], out["null"]
    return out


def fill_protocol(r, preds, val, ex):
    s = TPL.read_text()
    arch = "; ".join(f"Sector~{x['sector']}: $O-C = {x['o_minus_c_min']:+.1f}$~min, depth ${x['depth_ppm']:.0f} \\pm "
                     f"{x['depth_err_ppm']:.0f}$~ppm" for x in val["archival"])
    nul = val["null"]
    vtext = (f"Applied to the four archival transits of the data record, the code returns \\emph{{recovered}} for each "
             f"({arch}). Applied to transit-free data from Sector~101 shifted onto a predicted mid-time, it returns "
             f"\\emph{{excluded}} (depths ${nul[0]['depth_ppm']:.0f} \\pm {nul[0]['depth_err_ppm']:.0f}$ and "
             f"${nul[1]['depth_ppm']:.0f} \\pm {nul[1]['depth_err_ppm']:.0f}$~ppm at two offsets), and "
             f"\\emph{{insufficient data}} when the predicted transit falls in a data gap. Software: Python~3.13, "
             f"numpy~{np.__version__.rsplit('.', 1)[0]}, batman~2.5.1.")
    utc59 = next(x["utc"][11:] for x in ex["all_next"] if x["epoch"] == 59)
    rep = {"P": f"{r['P']:.5f}", "EP": f"{r['SIG_P']:.5f}", "T0": f"{r['T0']:.5f}", "ET0": f"{r['SIG_T0']:.5f}",
           "RMS": f"{r['RMS']:.1f}", "TC57": f"{preds[57][0]:.4f}", "S57": f"{preds[57][1]:.1f}",
           "TC59": f"{preds[59][0]:.4f}", "S59": f"{preds[59][1]:.1f}", "UTC59": utc59,
           "DEPTH": f"{r['DEPTH']:.0f}", "T14": f"{r['T14']:.2f}", "K": f"{r['K']:.5f}", "B": f"{r['B']:.3f}",
           "ARS": f"{r['ARS']:.2f}", "U1": f"{r['U'][0]:.3f}", "U2": f"{r['U'][1]:.3f}",
           "DLO": f"{0.5 * r['DEPTH']:.0f}", "DHI": f"{1.5 * r['DEPTH']:.0f}", "THREE_RMS": f"{3 * r['RMS']:.0f}",
           "VALIDATION": vtext}
    for k, v in rep.items():
        s = s.replace(f"@@{k}@@", v)
    left = re.findall(r"@@([A-Z0-9_]+)@@", s)
    assert not left, left
    (REC / "PROTOCOL.tex").write_text(s)
    subprocess.run(["latexmk", "-pdf", "-interaction=nonstopmode", "-quiet", "PROTOCOL.tex"], cwd=REC,
                   capture_output=True)
    subprocess.run(["latexmk", "-c"], cwd=REC, capture_output=True)
    assert (REC / "PROTOCOL.pdf").stat().st_mtime > (REC / "PROTOCOL.tex").stat().st_mtime - 1


def main():
    f = json.load(open(RESULTS / "report" / "transit_fits.json"))["TIC402710280"]
    ex = json.load(open(RESULTS / "report" / "extras.json"))["TIC402710280"]
    r = registered(f)
    write_test(r)
    m = load_test()
    preds = {e: m.predicted(e) for e in range(57, 67)}
    preds = {e: (tc, sig * 1440) for e, (tc, sig) in preds.items()}
    json.dump({"period_d": r["P"], "t0_bjd_tdb": r["T0"], "reference_epoch": r["REF"], "sigma_period_d": r["SIG_P"],
               "sigma_t0_d": r["SIG_T0"], "sigma_rms_min": round(r["RMS"], 2), "expected_depth_ppm": round(r["DEPTH"], 1),
               "expected_t14_h": round(r["T14"], 4), "registered_epochs": [57, 59],
               "predictions": [{"epoch": e, "tc_bjd_tdb": round(tc, 5), "sigma_pred_min": round(sg, 1)}
                               for e, (tc, sg) in preds.items()]},
              open(REC / "predictions.json", "w"), indent=1)
    val = validate(m)
    json.dump(val, open(REC / "validation.json", "w"), indent=1, default=float)
    fill_protocol(r, preds, val, ex)
    readme = (REC / "README.md").read_text()
    if "validation.json" not in readme:
        readme = readme.replace("| `predictions.json` |", "| `validation.json` | Output of the frozen code on the archival transits and on transit-free data (Section 7 of the protocol). |\n| `predictions.json` |")
    readme = readme.replace("Companion records: the preprint and its data/code record (R. Vora 2026, Zenodo).",
                            "Companion records: the preprint *A 36.34-day Sub-Neptune Candidate Transiting the Nearby K Dwarf "
                            "TIC 402710280, Identified in TESS Full-Frame Images* and its data/code record (R. Vora 2026, Zenodo).")
    (REC / "README.md").write_text(readme)
    files = sorted(p for p in REC.iterdir() if p.is_file() and p.name != "SHA256SUMS.txt")
    (REC / "SHA256SUMS.txt").write_text("".join(f"{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.name}\n" for p in files))
    zp = REC.parent / "prereg_tic402710280.zip"
    with zipfile.ZipFile(zp, "w", zipfile.ZIP_DEFLATED) as z:
        for p in sorted(REC.iterdir()):
            z.write(p, p.relative_to(REC.parent))
    print("registered:", {k: (round(v, 6) if isinstance(v, float) else v) for k, v in r.items()})
    print("epoch 57: %.4f +- %.1f min; epoch 59: %.4f +- %.1f min" % (*preds[57], *preds[59]))
    print("validation:", [(x["sector"], x["outcome"], round(x["o_minus_c_min"], 1)) for x in val["archival"]],
          [(x["case"], x["outcome"]) for x in val["null"]])


if __name__ == "__main__":
    main()
