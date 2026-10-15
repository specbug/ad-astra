"""Fill the TIC 402710280.01 manuscript (paper/tic402710280.template.tex -> paper/tic402710280.tex) with the numbers
of the stored analysis products, so that every quoted value traces to a result file.

Inputs: results/report/transit_fits.json (shape, times, ephemeris; fit_new.py), results/report/extras.json
(transit predictions with site times; predictions.py), results/candidates/audit_402710280.json (per-window
scatter), shape_free.json (free limb darkening; shape_chains.py), second_signal_402710280.json and
fpp_tic402710280_3.78d.json (the 3.78-day signal).

  python paper_numbers.py
"""
import json
import math
import re

import env  # noqa: F401
import numpy as np

from env import RESULTS

ROOT = RESULTS.parent
TPL = ROOT / "paper" / "tic402710280.template.tex"
OUT = ROOT / "paper" / "tic402710280.tex"
NUM = RESULTS / "candidates" / "paper_numbers.json"
R_STAR, E_R = 0.817, 0.047
M_STAR, E_M = 0.88, 0.11
TEFF, E_T = 5182.0, 121.0
RMS_EPOCHS = (57, 58, 59)
PRODUCT = {"tspoc": "TESS-SPOC", "qlp": "QLP", "ffi": "This work"}
CADENCE = {39: 600, 93: 200, 100: 200, 102: 200}


def pm(v, digits, sym_tol=0.15):
    """Median with asymmetric errors, as LaTeX; symmetric form when the two errors agree within sym_tol."""
    m, up, lo = v
    if abs(up - lo) <= sym_tol * max(up, lo):
        return f"{m:.{digits}f} \\pm {0.5 * (up + lo):.{digits}f}"
    return f"{m:.{digits}f}^{{+{up:.{digits}f}}}_{{-{lo:.{digits}f}}}"


def main():
    f = json.load(open(RESULTS / "report" / "transit_fits.json"))["TIC402710280"]
    ex = json.load(open(RESULTS / "report" / "extras.json"))["TIC402710280"]
    au = json.load(open(RESULTS / "candidates" / "audit_402710280.json"))
    sig = {x["sector"]: x["sigma_ppm"] for x in au["times_noS102_shape"]}
    P, eP = f["period_d"][0], f["period_d"][1]
    T0, eT = f["tref_btjd"][0] + 2457000, f["tref_btjd"][1]
    ref = f["ref_epoch"]
    oc = [t["o_c_min"] for t in f["times"]]
    rms = math.sqrt(np.mean(np.square(oc)))
    errs = [t["err_min"] for t in f["times"]]
    dep, edep = f["depth_transits_ppm"]
    k = f["rp_rs"]
    rp = k[0] * R_STAR * 109.08
    erp = rp * math.hypot(E_R / R_STAR, 0.5 * (k[1] + k[2]) / k[0])
    a_au = (M_STAR * (P / 365.25) ** 2) ** (1 / 3)
    ea = a_au * E_M / M_STAR / 3
    s = (R_STAR ** 2) * (TEFF / 5772) ** 4 / a_au ** 2
    es = s * math.sqrt((2 * E_R / R_STAR) ** 2 + (2 * ea / a_au) ** 2 + (4 * E_T / TEFF) ** 2)
    teq = TEFF * math.sqrt(R_STAR * 0.00465047 / (2 * a_au))
    eteq = teq * math.sqrt((E_T / TEFF) ** 2 + (E_R / R_STAR / 2) ** 2 + (ea / a_au / 2) ** 2)
    by = {t["sector"]: t for t in f["times"]}
    s100, s102 = by[100]["tc"], by[102]["tc"]
    rows = []
    for t in sorted(f["times"], key=lambda x: x["sector"]):
        sct = t["sector"]
        rows.append(f"{sct} & {PRODUCT[t['product']]} & {CADENCE[sct]} & {sig[sct]:.0f} & {t['beta']:.2f} & {t['epoch']} & "
                    f"${t['tc'] + 2457000:.4f} \\pm {t['err_min'] / 1440:.4f}$ & ${t['o_c_min']:+.1f}$ & "
                    f"${t['depth_ppm']:.0f} \\pm {t['depth_err_ppm']:.0f}$ \\\\")
    preds = []
    nxt = {x["epoch"]: x for x in ex["all_next"]}
    notes = {57: "\\tess\\ Sector~114", 58: "daytime at La Silla", 59: "full transit at night at La Silla"}
    for e in RMS_EPOCHS:
        n = e - ref
        tc = T0 + n * P
        sp = math.sqrt(eT ** 2 + (n * eP) ** 2 + (rms / 1440) ** 2) * 1440
        x = nxt[e]
        assert abs(x["btjd"] + 2457000 - tc) < 1e-4, "extras.json predates the current ephemeris"
        y, mo, d = x["utc"][:10].split("-")
        mon = ["Jan.", "Feb.", "Mar.", "Apr.", "May", "June", "July", "Aug.", "Sep.", "Oct.", "Nov.", "Dec."][int(mo) - 1]
        preds.append(f"{e} & {tc:.3f} & {sp:.0f} & {y} {mon} {int(d)}, {x['utc'][11:]} & {notes[e]} \\\\")
    rep = {
        "P": f"{P:.5f} \\pm {eP:.5f}", "T0": f"{T0:.4f} \\pm {eT:.4f}", "RMS": f"{rms:.1f}",
        "DEPTH": f"{dep:.0f} \\pm {edep:.0f}", "DEPTH_CHI2": f"{f['depth_chi2'][0]:.1f}",
        "T14": pm(f["t14_h"], 2), "K": pm(k, 4), "K_FULL": pm(k, 5), "B": pm(f["b"], 2), "RHO": pm(f["rho_solar"], 2),
        "ARS": pm(f["a_rs"], 1), "RP": f"{rp:.2f} \\pm {erp:.2f}", "A_AU": f"{a_au:.3f} \\pm {ea:.3f}",
        "A_AU_SHORT": f"{a_au:.3f}", "S": f"{s:.1f} \\pm {es:.1f}", "TEQ": f"{teq:.0f} \\pm {eteq:.0f}",
        "B95": f"{f['b_95']:.2f}", "NTAU": f"{int(f['chain_summary']['steps_over_tau'] // 10 * 10)}",
        "REDCHI2": f"{f['chi2_min'] / f['n_points']:.2f}", "NPTS": f"{f['n_points']}",
        "OC_LIST": ", ".join(f"${x:+.1f}$" for x in oc[:-1]) + f" and ${oc[-1]:+.1f}$",
        "ERR_RANGE": f"{min(errs):.1f}--{max(errs):.1f}", "EPH_CHI2": f"{f['ephemeris_chi2'][0]:.0f}",
        "INTERVAL_EXCESS": f"{((s102 - s100) - P) * 1440:.1f}",
        "TIMES_ROWS": "\n".join(rows), "PRED_ROWS": "\n".join(preds),
    }
    rel = ROOT / "release" / "tic402710280_zenodo" / "data" / "results"
    import csv
    al = list(csv.DictReader(open(rel / "alias_test.csv")))
    assert abs(float(al[0]["period_d"]) - P / 2) < 1e-4, "alias_test.csv predates the current ephemeris"
    sig_of = {a: [float(r["exclusion_sigma"]) for r in al if r["alias"] == a] for a in ("P/2", "P/3", "P/4", "P/5")}
    p2_secs = [r["sector"] for r in al if r["alias"] == "P/2"]
    rep["ALIAS_MIN"] = f"{min(min(v) for v in sig_of.values()):.1f}"
    rep["P2_MIN"], rep["P3_MIN"] = f"{min(sig_of['P/2']):.1f}", f"{min(sig_of['P/3']):.1f}"
    words = {2: "two", 3: "three", 4: "four", 5: "five", 6: "six", 7: "seven", 8: "eight", 9: "nine", 10: "ten"}
    rep["N_P2"] = words.get(len(sig_of["P/2"]), str(len(sig_of["P/2"])))
    rep["N_P3"] = words.get(len(sig_of["P/3"]), str(len(sig_of["P/3"])))
    rep["N_P45"] = str(len(sig_of["P/4"]) + len(sig_of["P/5"]))
    rep["P2_SECTORS"] = ", ".join(p2_secs[:-1]) + (" and~" if len(p2_secs) > 1 else "") + p2_secs[-1]
    sec = list(csv.DictReader(open(rel / "secondary_eclipse.csv")))
    d_ = np.array([float(r["depth_ppm"]) for r in sec])
    e_ = np.array([float(r["err_ppm"]) for r in sec])
    w_ = 1 / e_ ** 2
    sd, se = float(np.sum(w_ * d_) / w_.sum()), float(1 / math.sqrt(w_.sum()))
    rep["SEC_DEPTH"], rep["SEC_3SIG"] = f"{sd:.0f} \\pm {se:.0f}", f"{sd + 3 * se:.0f}"
    s66 = au["s66_fit"]
    rep["S66_DEPTH"] = f"{s66['depth_ppm']:.0f} \\pm {s66['depth_err_ppm']:.0f}"
    fl = json.load(open(RESULTS / "candidates" / "shape_free.json"))
    rep["K_FREELD"] = "$" + pm(fl["rp_rs"], 4) + "$"
    rep["B_FREELD"] = "$" + pm(fl["b"], 2) + "$"
    ss = json.load(open(RESULTS / "candidates" / "second_signal_402710280.json"))
    sf = ss["fit_free_rho"]
    rep["SS_P"] = f"{ss['bls_joint']['P']:.5f}"
    rep["SS_DEPTH"] = f"{ss['depth_all'][0]:.0f} \\pm {ss['depth_all'][1]:.0f}"
    rep["SS_T14"] = f"{sf['t14_h'][0]:.1f}"
    dd = np.array([r["depth_ppm"] for r in ss["sectors_depth"]])
    de = np.array([r["err_ppm"] for r in ss["sectors_depth"]])
    mw = np.sum(dd / de ** 2) / np.sum(1 / de ** 2)
    g2 = ss["epoch_groups"]["2"]
    rep["SS_ODD"], rep["SS_EVEN"] = f"{g2[1][0]:.0f} \\pm {g2[1][1]:.0f}", f"{g2[0][0]:.0f} \\pm {g2[0][1]:.0f}"
    sc = ss["secondary"]["phase0.5_ppm"]
    rep["SS_SEC"] = f"{sc[0]:.0f} \\pm {sc[1]:.0f}"
    rep["SS_NEIGH"] = f"{ss['neighbour']['required_eclipse_depth'] * 100:.0f}"
    rep["SS_SECTOR_CHI2"] = f"{np.sum(((dd - mw) / de) ** 2):.1f}"
    rep["SS_FIT"] = (f"$R_p/R_\\star = {pm(sf['rp_rs'], 4)}$, $b = {pm(sf['b'], 2)}$ and "
                     f"$\\rho_\\star = {pm(sf['rho_solar'], 1)}\\,\\rho_\\odot$")
    rr = ss["radius_earth_if_on_target"]
    chs = np.load(RESULTS / "candidates" / "second_signal_chain.npy")
    low = 10 ** chs[:, 3] < 1.61 + 2 * 0.38
    rep["SS_INTERP"] = SS_INTERP.format(rp=rr[0], erp=rr[1], rho=pm(sf["rho_solar"], 1, sym_tol=0.0),
                                        frac=100 * low.mean(), bmin=float(np.percentile(chs[low, 2], 2.5)))
    fp = RESULTS / "candidates" / "fpp_tic402710280_3.78d.json"
    if fp.exists():
        d = json.load(open(fp))
        pl = 0.0
        for r in d["runs"]:
            pl += sum(v for k_, v in r["scenarios"].items() if k_.split("@")[0] in ("TP", "PTP", "STP")) / len(d["runs"])
        rep["SS_FPP"] = (f"\\texttt{{TRICERATOPS}}, applied as in \\S\\ref{{sec:fpp}}, assigns {100 * pl:.0f}\\% of the "
                         f"probability to a planet around the target or an unresolved companion, but gives "
                         f"${{\\rm FPP}} = {d['FPP_mean']:.2f} \\pm {d['FPP_std']:.2f}$ and ${{\\rm NFPP}} = "
                         f"{d['NFPP_mean']:.3f} \\pm {d['NFPP_std']:.3f}$, above the limit of 0.001 below which "
                         f"\\citet{{Giacalone2021}} consider nearby false positives excluded.")
    else:
        rep["SS_FPP"] = ""
    tex = open(TPL).read()
    for key, val in rep.items():
        tex = tex.replace(f"@@{key}@@", val)
    left = sorted(set(re.findall(r"@@([A-Z0-9_]+)@@", tex)))
    assert not left, f"unfilled: {left}"
    open(OUT, "w").write(tex)
    json.dump(rep, open(NUM, "w"), indent=1)
    print("wrote", OUT)
    for key in ("P", "T0", "RMS", "DEPTH", "K_FULL", "B", "B95", "RHO", "RP", "TEQ", "OC_LIST", "K_FREELD", "B_FREELD",
                "SS_DEPTH", "SS_FIT", "SS_FPP"):
        print(f"  {key}: {rep[key]}")


SS_INTERP = ("On the target, the signal would correspond to a planet of ${rp:.2f} \\pm {erp:.2f}\\,\\rearth$ on an inner "
             "orbit. The fitted density, ${rho}\\,\\rho_\\odot$, exceeds the target's $1.61 \\pm 0.38\\,\\rho_\\odot$: "
             "{frac:.0f}\\% of the posterior lies below $2.37\\,\\rho_\\odot$ (the target's value plus $2\\sigma$), almost all of "
             "it at $b > {bmin:.2f}$, so a transit across the target is consistent with the short duration only if it is "
             "nearly grazing.")


if __name__ == "__main__":
    main()
