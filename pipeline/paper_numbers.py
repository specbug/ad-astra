"""Fill the TIC 402710280.01 manuscript (paper/tic402710280.template.tex -> paper/tic402710280.tex) with the numbers
of the stored analysis products, so that every quoted value traces to a result file.

Inputs: results/report/transit_fits.json (shape, times, ephemeris; fit_new.py), results/report/extras.json
(transit predictions with site times; predictions.py), results/candidates/audit_402710280.json (per-window
scatter), shape_free.json (free limb darkening; shape_chains.py), second_signal_402710280.json and
fpp_tic402710280_3.78d.json (the 3.78-day signal); qlp_recheck_402710280.json (version 2: the QLP light curves of
Sectors 100 and 102, our photometry of Sectors 39 and 93, and Sector 66 with the final shape).

  python paper_numbers.py
"""
import json
import math
import re

import env  # noqa: F401
import numpy as np

import transit_fit as tf
from scipy.optimize import brentq

from env import RESULTS

ROOT = RESULTS.parent
TPL = ROOT / "paper" / "tic402710280.template.tex"
OUT = ROOT / "paper" / "tic402710280.tex"
NUM = RESULTS / "candidates" / "paper_numbers.json"
R_STAR, E_R = 0.817, 0.047
M_STAR, E_M = 0.88, 0.11
RHO_TIC, E_RHO = 1.61, 0.38  # solar units, from the TIC mass and radius (Section 5.2)
TEFF, E_T = 5182.0, 121.0
J_MAG = 8.528  # 2MASS J of TIC 402710280, from TIC v8.2
RMS_EPOCHS = (57, 58, 59)
PRODUCT = {"tspoc": "TESS-SPOC", "qlp": "QLP", "ffi": "This work"}
CADENCE = {39: 600, 93: 200, 100: 200, 102: 200}


def longest_transit_h(rho_sun, period_d, k, q=None):
    """Longest transit (b = 0) of a planet with radius ratio k on a period_d orbit around a star of density rho_sun:
    centred at apastron of an eccentric orbit with periastron at q stellar radii, or circular if q is None."""
    a_r = (6.674e-11 * rho_sun * 1408.9 * (period_d * 86400) ** 2 / (3 * math.pi)) ** (1 / 3)
    e = 0.0 if q is None else max(0.0, 1 - q / a_r)
    mean_anom = np.linspace(math.pi - 1.5, math.pi + 1.5, 300001)
    ecc = mean_anom.copy()
    for _ in range(60):
        ecc -= (ecc - e * np.sin(ecc) - mean_anom) / (1 - e * np.cos(ecc))
    r = a_r * (1 - e * np.cos(ecc))
    nu = 2 * np.arctan2(math.sqrt(1 + e) * np.sin(ecc / 2), math.sqrt(1 - e) * np.cos(ecc / 2))
    x = np.abs(r * np.sin(nu - math.pi))                 # sky-projected distance from the star's centre, apastron at 0
    return float(np.sum(x < 1 + k) * (mean_anom[1] - mean_anom[0]) * period_d * 24 / (2 * math.pi))


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
        sp = tf.sigma_pred_min(f, e, rms)   # with the T0-P covariance
        x = nxt[e]
        assert abs(x["btjd"] + 2457000 - tc) < 1e-4, "extras.json predates the current ephemeris"
        y, mo, d = x["utc"][:10].split("-")
        mon = ["Jan.", "Feb.", "Mar.", "Apr.", "May", "June", "July", "Aug.", "Sep.", "Oct.", "Nov.", "Dec."][int(mo) - 1]
        preds.append(f"{e} & {tc:.4f} & {sp:.0f} & {y} {mon} {int(d)}, {x['utc'][11:]} & {notes[e]} \\\\")
    rep = {
        "P": f"{P:.5f} \\pm {eP:.5f}", "T0": f"{T0:.4f} \\pm {eT:.4f}", "RMS": f"{rms:.1f}",
        "DEPTH": f"{dep:.0f} \\pm {edep:.0f}", "DEPTH_CHI2": f"{f['depth_chi2'][0]:.1f}",
        "T14": pm(f["t14_h"], 2), "K": pm(k, 4), "K_FULL": pm(k, 5), "B": pm(f["b"], 2), "RHO": pm(f["rho_solar"], 2),
        "ARS": pm(f["a_rs"], 1), "RP": f"{rp:.2f} \\pm {erp:.2f}", "A_AU": f"{a_au:.3f} \\pm {ea:.3f}",
        "A_AU_SHORT": f"{a_au:.3f}", "S": f"{s:.1f} \\pm {es:.1f}", "TEQ": f"{teq:.0f} \\pm {eteq:.0f}",
        "B95": f"{f['b_95']:.2f}", "NTAU": f"{int(f['chain_summary']['steps_over_tau'] // 10 * 10)}",
        "REDCHI2": f"{f['chi2_min'] / f['n_points']:.2f}", "NPTS": f"{f['n_points']}",
        "OC_LIST": ", ".join(f"${x:+.1f}$" for x in oc[:-1]) + f", and ${oc[-1]:+.1f}$",
        "ERR_RANGE": f"{min(errs):.1f}--{max(errs):.1f}", "EPH_CHI2": f"{f['ephemeris_chi2'][0]:.0f}",
        "INTERVAL_EXCESS": f"{((s102 - s100) - P) * 1440:.1f}",
        "TIMES_ROWS": "\n".join(rows), "PRED_ROWS": "\n".join(preds),
        "DEP93": f"{by[93]['depth_ppm']:.0f} \\pm {by[93]['depth_err_ppm']:.0f}",
    }
    # transmission spectroscopy metric (Kempton et al. 2018; scale factor 1.28 for 2.75-4 Earth radii), with the
    # Chen & Kipping (2017) mass of about 14 Earth masses quoted in the Discussion and the 2MASS J magnitude from the TIC
    rep["TSM"] = f"{round(1.28 * rp ** 3 * teq / (14 * R_STAR ** 2) * 10 ** (-J_MAG / 5), -1):.0f}"
    # the T0-P correlation of the linear ephemeris (reference epoch 40; the weighted mean epoch is earlier)
    _, _, cov = tf.ephemeris(f)
    w_ep = 1 / np.square(errs)
    rep["EPOCH_WMEAN"] = f"{np.sum(w_ep * np.array([t['epoch'] for t in f['times']])) / w_ep.sum():.0f}"
    rep["T0P_CORR"] = f"{cov[0, 1] / math.sqrt(cov[0, 0] * cov[1, 1]):.2f}"
    # the registered test (doi:10.5281/zenodo.23390325; a copy of its predictions.json) searches +-3 sigma_pred, with
    # sigma_pred computed without the T0-P covariance
    reg = json.load(open(RESULTS / "candidates" / "prereg_predictions.json"))
    s_reg = next(x["sigma_pred_min"] for x in reg["predictions"] if x["epoch"] == 57)
    s_now = tf.sigma_pred_min(f, 57, rms)
    rep["PREREG_SIG57"], rep["SIG57"], rep["PREREG_NSIG"] = f"{s_reg:.1f}", f"{s_now:.1f}", f"{3 * s_reg / s_now:.1f}"
    rel = RESULTS / "candidates"            # written there by build_release.py
    import csv
    al = list(csv.DictReader(open(rel / "alias_test.csv")))
    assert abs(float(al[0]["period_d"]) - P / 2) < 1e-4, "alias_test.csv predates the current ephemeris"
    sig_of = {a: [float(r["exclusion_sigma"]) for r in al if r["alias"] == a] for a in ("P/2", "P/3", "P/4", "P/5")}
    # P/6 to P/36: each alias, all its covered windows combined
    hi = {}
    for a in sorted({r["alias"] for r in al}, key=lambda x: int(x[2:])):
        if int(a[2:]) < 6:
            continue
        rows = [r for r in al if r["alias"] == a]
        dd = np.array([float(r["depth_ppm"]) for r in rows])
        ee = np.array([float(r["err_ppm"]) for r in rows])
        ww = 1 / ee ** 2
        hi[a] = (len(rows), (dep - np.sum(ww * dd) / ww.sum()) / math.hypot(1 / math.sqrt(ww.sum()), edep))
    kmax = max(int(a[2:]) for a in hi)
    rep["ALIAS_KMAX"], rep["ALIAS_PMIN"] = str(kmax), f"{P / kmax:.1f}"
    rep["ALIAS_HI_N"] = f"{min(v[0] for v in hi.values())} to {max(v[0] for v in hi.values())}"
    rep["ALIAS_HI_SIG"] = f"{math.floor(min(v[1] for v in hi.values()))}"
    # periods below the last tested fraction: the longest transit they allow, against the observed duration
    p_sub, t14 = P / (kmax + 1), f["t14_h"][0]
    rep["SUB_TMAX"] = f"{longest_transit_h(RHO_TIC, p_sub, k[0], q=1.0):.1f}"
    rho_lim = brentq(lambda r_: longest_transit_h(r_, p_sub, k[0], q=1.0) - t14, 0.1, RHO_TIC)
    rep["SUB_RHO_SIG"] = f"{(RHO_TIC - rho_lim) / E_RHO:.1f}"
    rep["BETA_MAX"] = f"{max(t['beta'] for t in f['times']):.1f}"
    p2_secs = [r["sector"] for r in al if r["alias"] == "P/2"]
    rep["ALIAS_MIN"] = f"{min(min(v) for v in sig_of.values()):.1f}"   # P/2 to P/5, every window individually
    rep["P2_MIN"], rep["P3_MIN"] = f"{min(sig_of['P/2']):.1f}", f"{min(sig_of['P/3']):.1f}"
    words = {2: "two", 3: "three", 4: "four", 5: "five", 6: "six", 7: "seven", 8: "eight", 9: "nine", 10: "ten"}
    rep["N_P2"] = words.get(len(sig_of["P/2"]), str(len(sig_of["P/2"])))
    rep["N_P3"] = words.get(len(sig_of["P/3"]), str(len(sig_of["P/3"])))
    rep["N_P45"] = str(len(sig_of["P/4"]) + len(sig_of["P/5"]))
    rep["P2_SECTORS"] = (", ".join(p2_secs[:-1]) + ("," if len(p2_secs) > 2 else "") + (" and~" if len(p2_secs) > 1 else "")
                         + p2_secs[-1])
    sec = list(csv.DictReader(open(rel / "secondary_eclipse.csv")))
    d_ = np.array([float(r["depth_ppm"]) for r in sec])
    e_ = np.array([float(r["err_ppm"]) for r in sec])
    w_ = 1 / e_ ** 2
    sd, se = float(np.sum(w_ * d_) / w_.sum()), float(1 / math.sqrt(w_.sum()))
    rep["SEC_DEPTH"], rep["SEC_3SIG"] = f"{sd:.0f} \\pm {se:.0f}", f"{sd + 3 * se:.0f}"
    # version 2: Sector 66 with the final shape, and the QLP light curves of Sectors 100 and 102 (qlp_recheck_402710280.py)
    qc = json.load(open(RESULTS / "candidates" / "qlp_recheck_402710280.json"))
    s66 = qc["s66"]
    fb = s66["free_time_by_baseline"]
    assert not any(v["at_grid_edge"] for v in fb.values())
    rep["S66_FIRST_H"] = f"{-s66['first_quality0_rel_mid_h']:.1f}"
    rep["S66_PRE_INGRESS_H"] = f"{-s66['first_quality0_rel_ingress_h']:.1f}"
    r10 = lambda x: f"{10 * round(x / 10):.0f}"   # depths and errors of ~100 ppm: nearest 10 ppm
    rep["S66_DEPTHS"] = ", ".join(r10(fb[b]["depth_ppm"]) for b in ("constant", "linear")) + f", and {r10(fb['quadratic']['depth_ppm'])}"
    rep["S66_DERR"] = f"{r10(min(v['depth_err_ppm'] for v in fb.values()))}--{r10(max(v['depth_err_ppm'] for v in fb.values()))}"
    ocs = [fb[b]["o_c_min"] for b in ("constant", "linear", "quadratic")]
    rep["S66_OC_RANGE"] = f"${min(ocs):+.0f}$ to ${max(ocs):+.0f}$"
    rq = qc["retime_qlp"]
    q100, q102 = rq["100"], rq["102"]
    rep["QLP_OFF100"] = f"{-q100['time_offset']['offset_at_transit_s']:.0f}"
    rep["QLP_OFF102"] = f"{-q102['time_offset']['offset_at_transit_s']:.0f}"
    rep["QLP_OFF93"] = f"{math.ceil(abs(rq['93']['time_offset']['offset_at_transit_s'])):.0f}"
    assert q100["minus_published_min"] < 0 and q102["minus_published_min"] < 0
    rep["QLP_EARLY100"], rep["QLP_EARLY102"] = f"{-q100['minus_published_min']:.1f}", f"{-q102['minus_published_min']:.1f}"
    rep["QLP_SIG100"], rep["QLP_SIG102"] = f"{-q100['sigma_of_difference']:.1f}", f"{-q102['sigma_of_difference']:.1f}"
    rep["QLP_INTERVAL"] = f"{((q102['tc'] - q100['tc']) - P) * 1440:.1f}"
    rep["QLP_DEP100"] = f"{q100['depth_ppm']:.0f} \\pm {q100['depth_err_ppm']:.0f}"
    rep["QLP_DEP102"] = f"{q102['depth_ppm']:.0f} \\pm {q102['depth_err_ppm']:.0f}"
    assert q100["depth_minus_published_sigma"] < 0 and q102["depth_minus_published_sigma"] < 0
    rep["QLP_DSIG100"], rep["QLP_DSIG102"] = f"{-q100['depth_minus_published_sigma']:.1f}", f"{-q102['depth_minus_published_sigma']:.1f}"
    rep["OWN_MAX"] = f"{max(abs(v['minus_published_min']) for v in qc['own_control'].values()):.1f}"
    # sensitivity of the depth, radius, period and timing scatter to adopting the QLP measurements of Sectors 100 and 102
    q_ = qc["retime_qlp"]
    dd_ = {s_: (by[s_]["depth_ppm"], by[s_]["depth_err_ppm"]) for s_ in (39, 93)}
    dd_.update({s_: (q_[str(s_)]["depth_ppm"], q_[str(s_)]["depth_err_ppm"]) for s_ in (100, 102)})
    v_ = np.array([x[0] for x in dd_.values()])
    w_ = 1 / np.array([x[1] for x in dd_.values()]) ** 2
    mq = float(np.sum(w_ * v_) / w_.sum())
    rep["QMIX_DEPTH"] = f"{mq:.0f} \\pm {1 / math.sqrt(w_.sum()):.0f}"
    rep["QMIX_CHI2"] = f"{float(np.sum(w_ * (v_ - mq) ** 2)):.0f}"
    rep["QMIX_RP"] = f"{rp * math.sqrt(mq / dep):.2f}"
    ep_ = np.array([by[s_]["epoch"] for s_ in (39, 93, 100, 102)], float) - ref
    tt_ = np.array([by[39]["tc"], by[93]["tc"], q_["100"]["tc"], q_["102"]["tc"]])
    st_ = np.array([by[39]["err_min"], by[93]["err_min"], q_["100"]["err_min"], q_["102"]["err_min"]]) / 1440
    A_ = np.vstack([np.ones(4), ep_]).T
    W_ = np.diag(1 / st_ ** 2)
    pq = np.linalg.solve(A_.T @ W_ @ A_, A_.T @ W_ @ tt_)
    rep["QMIX_DP"] = f"{abs(pq[1] - P) * 86400:.0f}"
    rep["QMIX_RMS"] = f"{math.sqrt(np.mean(((tt_ - A_ @ pq) * 1440) ** 2)):.1f}"
    fv2 = json.load(open(RESULTS / "candidates" / "fpp_tic402710280_v2.json"))
    runs = [r_["FPP"] for r_ in fv2["runs"]]
    rep["FPP_RANGE"] = f"{min(runs):.3f}--{max(runs):.3f}"
    rep["FPP_BOX5"] = f"{json.load(open(RESULTS / 'candidates' / 'fpp_tic402710280_v2_box5.json'))['FPP_mean']:.3f}"
    rep["FPP_LINFOLD"] = f"{json.load(open(RESULTS / 'candidates' / 'fpp_tic402710280.json'))['FPP_mean']:.3f}"
    from archive_comparison import load as load_archive
    adate, arows = load_archive()
    words_ = {1: "one", 2: "two", 3: "three", 4: "four", 5: "five", 6: "six", 7: "seven", 8: "eight", 9: "nine"}
    longp = [x for x in arows if float(x["pl_orbper"]) > 30]
    assert len(longp) == 1 and longp[0]["pl_name"] == "TOI-2076 d", "the introduction names TOI-2076 d"
    rep["ARCH_N"] = words_.get(len(arows), str(len(arows)))
    rep["ARCH_P"] = f"{float(longp[0]['pl_orbper']):.1f}"
    y, m, d = (int(x) for x in adate.split("-"))
    rep["ARCH_DATE"] = f"{y} {['January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'December'][m - 1]}~{d}"
    tb = json.load(open(RESULTS / "candidates" / "qlp_timebase_402710280.json"))
    # epoch at which the target's own barycentric correction equals the single correction QLP applied to the sector:
    # QLP - true = C - ltt(t), and C = ltt(tc) + offset at the transit
    from qlp_timebase_402710280 import ltt
    zero = []
    for s_ in ("100", "102"):
        q = qc["retime_qlp"][s_]
        c = ltt(q["tc"]) + q["time_offset"]["offset_at_transit_s"]
        grid = np.arange(q["tc"] - 40, q["tc"] + 0.01, 0.05)
        lt = np.array([ltt(t) for t in grid])
        zero.append(grid[np.argmin(np.abs(lt - c))])
    rep["QLP_ZERO100"], rep["QLP_ZERO102"] = f"{zero[0]:.0f}", f"{zero[1]:.0f}"
    lin = max(tb[k]["departure_from_line_s"]["range"] for k in ("100", "102"))
    rep["QLP_LIN"] = f"{math.ceil(lin * 100) / 100:.2f}"
    rep["QLP_SLOPES"] = f"${tb['100']['offset_drift_s_per_d']:.1f}$ and ${tb['102']['offset_drift_s_per_d']:.1f}$"
    rep["BARY_RATES"] = f"{tb['100']['bary_rate_s_per_d']:.1f} and {tb['102']['bary_rate_s_per_d']:.1f}"
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
    rep["SS_FIT"] = (f"$R_p/R_\\star = {pm(sf['rp_rs'], 4)}$, $b = {pm(sf['b'], 2)}$, and "
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
                         f"probability to a planet around the target or an unresolved companion but gives "
                         f"${{\\rm FPP}} = {d['FPP_mean']:.2f} \\pm {d['FPP_std']:.2f}$ and ${{\\rm NFPP}} = "
                         f"{d['NFPP_mean']:.3f} \\pm {d['NFPP_std']:.3f}$, above the NFPP limit of 0.001 below which "
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
             "orbit. The fitted density, ${rho}\\,\\rho_\\odot$, exceeds the target's $1.61 \\pm 0.38\\,\\rho_\\odot$. Only "
             "{frac:.0f}\\% of the posterior lies below $2.37\\,\\rho_\\odot$ (the target's value plus $2\\sigma$), almost all of "
             "it at $b > {bmin:.2f}$, so a transit of the target matches the short duration only if it is nearly grazing.")


if __name__ == "__main__":
    main()
