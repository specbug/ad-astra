"""Transit fits of the two candidates found by the period hunt, for the report (results/report/transit_fits.json,
keys TIC402710280 and TOI2493), and their derived quantities and upcoming transits (results/report/extras.json).

TIC 402710280.01 (v2, after the independent audit of 15 October 2026): four transits (S39 TESS-SPOC, S93 QLP, S100
and S102 from our own full-frame photometry, ffi_new.py). The S66 window lies on a stray-light ramp with data
starting 0.4 h before ingress (audit_402710280.py): shown, not used. The mid-times scatter by tens of minutes about a
straight line, so the shape is fitted with each window centred on its own measured mid-time
(transit_fit.fit_stacked_windows) and the ephemeris is a weighted line through the times with errors inflated by
sqrt(chi2_red). v1 fitted S39, S66, S93 and S100 on one linear ephemeris, which put the S66 ramp into the depth and
smeared the shape (the low density once blamed on the S102 ramp came from that smear).
TOI-2493.01: three transits in SPOC 2-min data (S33, S98, S107).

  python fit_new.py
"""
import json

import env  # noqa: F401
import numpy as np

import report_extras as rx
import ticdb
import transit_fit as tf
from dossier import all_sectors, best_lcs
from env import RESULTS
from ffi_new import new_lcs

FITS = RESULTS / "report" / "transit_fits.json"
EXTRAS = RESULTS / "report" / "extras.json"


def lcs_for(tic):
    st = ticdb.query([tic], ["ra", "dec"]).iloc[0]
    lcs = best_lcs(tic, all_sectors(tic))
    have = {lc["meta"]["sector"] for lc in lcs}
    return lcs + [lc for lc in new_lcs(tic, float(st.ra), float(st.dec), build_missing=False) if lc["meta"]["sector"] not in have]


NSTEP = 30000     # autocorrelation times reach ~450 steps (b-rho-k degeneracy); 30,000 steps is >50 of them


def main():
    fits = json.load(open(FITS))
    ex = json.load(open(EXTRAS))
    # TIC 402710280.01 (v2, after the independent audit of 15 October 2026; audit_402710280.py). The sector-66
    # window lies on a stray-light ramp (data start 0.4 h before ingress) and is shown but not used; the shape is
    # fitted to sectors 39, 93 and 100 (sector 102 sits on a ramp at the start of an orbit and distorts the shape);
    # the ephemeris comes from the four measured mid-times, its errors inflated by sqrt(chi2_red) because the times
    # scatter by tens of minutes about a straight line.
    tic = 402710280
    lcs = lcs_for(tic)
    au = json.load(open(RESULTS / "candidates" / "audit_402710280.json"))
    meas = {x["sector"]: x for x in au["times_noS102_shape"]}
    ep_of = {39: 0, 66: 20, 93: 40, 100: 47, 102: 48}
    from audit_402710280 import beta_factor
    data = []
    for w in tf.windows_from_lcs(lcs, [meas[s]["tc"] for s in (39, 93, 100, 102)], 0.8, 5.3):
        oot = np.abs(w[0] - w[4]) > 0.75 * 5.3 / 24
        data.append((w[0], w[1], w[2] * beta_factor(w[0][oot], w[1][oot] - 1, w[2]), w[3], w[4]))
    out, ch = tf.fit_stacked_windows(data, 36.3399, 5182.0, nstep=NSTEP)
    ref = 40
    out["ref_epoch"] = ref
    n = np.array([ep_of[s] for s in (39, 93, 100, 102)], float)
    tt = np.array([meas[s]["tc"] for s in (39, 93, 100, 102)])
    st = np.array([meas[s]["err_min"] for s in (39, 93, 100, 102)]) / 1440
    A = np.vstack([np.ones_like(n), n - ref]).T
    W = np.diag(1 / st ** 2)
    cov = np.linalg.inv(A.T @ W @ A)
    p = cov @ A.T @ W @ tt
    chi2 = float(np.sum(((tt - A @ p) / st) ** 2))
    sc = np.sqrt(max(chi2 / 2, 1.0))
    eT, eP = float(np.sqrt(cov[0, 0]) * sc), float(np.sqrt(cov[1, 1]) * sc)
    out["tref_btjd"] = [float(p[0]), eT, eT]
    out["period_d"] = [float(p[1]), eP, eP]
    out["ephemeris_chi2"] = [chi2, 2]
    out["times"] = [{"sector": s, "epoch": ep_of[s], "tc": meas[s]["tc"], "err_min": meas[s]["err_min"],
                     "depth_ppm": meas[s]["depth_ppm"], "depth_err_ppm": meas[s]["depth_err_ppm"],
                     "product": meas[s]["product"],
                     "o_c_min": float((meas[s]["tc"] - p[0] - (ep_of[s] - ref) * p[1]) * 1440)} for s in (39, 93, 100, 102)]
    out["depth_transits_ppm"] = au["per_transit_depth_mean"]
    med = np.median(ch, axis=0)
    aR = ((10 ** med[2]) * tf.RHO_SUN * tf.G * (p[1] * 86400) ** 2 / (3 * np.pi)) ** (1 / 3)
    inc = np.degrees(np.arccos(med[1] / aR))
    centres = [meas[39]["tc"], float(p[0] + (20 - ref) * p[1]), meas[93]["tc"], meas[100]["tc"], meas[102]["tc"]]
    pdata = tf.windows_from_lcs(lcs, centres, 0.8, 5.3)
    out["plot"] = tf.plot_data([(w[0], w[1], w[2]) for w in pdata], centres, float(p[1]), med[0], aR, inc, out["u"])
    out["sectors"] = [39, 66, 93, 100, 102]
    out["fit_windows"] = [39, 93, 100, 102]
    fits["TIC402710280"] = out
    print("TIC 402710280", {k: v for k, v in out.items() if k != "plot"}, flush=True)
    # TOI-2493.01
    tic = 123664207
    lcs = lcs_for(tic)
    data = tf.windows_from_lcs(lcs, [2225.82, 4043.534, 4280.627], 0.8, 5.0)
    out, ch = tf.fit_periodic_windows(data, [0, 23, 26], 79.03019, 4242.0, nstep=NSTEP)
    out["sectors"] = [33, 98, 107]
    fits["TOI2493"] = out
    print("TOI-2493", {k: v for k, v in out.items() if k != "plot"}, flush=True)
    fresh = json.load(open(FITS))           # other fits may have been written meanwhile (fit_headline.py)
    fresh.update({k: fits[k] for k in ("TIC402710280", "TOI2493")})
    json.dump(fresh, open(FITS, "w"), indent=1, default=float)
    ex = json.load(open(EXTRAS))
    ex["TIC402710280"] = rx.periodic(fits, "TIC402710280", 0.88, 0.817, 5182.0, 248.524, -63.55, now=4325.5)
    ex["TOI2493"] = rx.periodic(fits, "TOI2493", 0.66, 0.652, 4242.0, 89.74609, -21.57558, now=4325.5)
    json.dump(ex, open(EXTRAS, "w"), indent=1, default=float)
    for k in ("TIC402710280", "TOI2493"):
        print(k, {a: b for a, b in ex[k].items() if a not in ("next", "all_next")})
        for x in ex[k]["next"]:
            print("   ", x["utc"], x["note"])


if __name__ == "__main__":
    main()
