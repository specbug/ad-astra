"""Transit fit of TIC 402710280.01 (results/report/transit_fits.json, key TIC402710280), and its derived quantities
and upcoming transits (results/report/extras.json).

TIC 402710280.01 (v2, after the independent audit of 15 October 2026): four transits (S39 TESS-SPOC, S93 QLP, S100
and S102 from our own full-frame photometry, ffi_new.py). The S66 window lies on a stray-light ramp with data
starting 0.4 h before ingress (audit_402710280.py): shown, not used. The mid-times scatter by tens of minutes about a
straight line, so the shape is fitted with each window centred on its own measured mid-time
(transit_fit.fit_stacked_windows) and the ephemeris is a weighted line through the times with errors inflated by
sqrt(chi2_red). v1 fitted S39, S66, S93 and S100 on one linear ephemeris, which put the S66 ramp into the depth and
smeared the shape (the low density once blamed on the S102 ramp came from that smear).

  python fit_new.py
"""
import json

import env  # noqa: F401
import numpy as np

import predictions
import transit_fit as tf
from env import RESULTS
from lcaccess import lcs_for

FITS = RESULTS / "report" / "transit_fits.json"
EXTRAS = RESULTS / "report" / "extras.json"


NSTEP = 30000     # autocorrelation times reach ~450 steps (b-rho-k degeneracy); 30,000 steps is >50 of them


def main():
    fits = json.load(open(FITS))
    ex = json.load(open(EXTRAS))
    # TIC 402710280.01 (v2, after the independent audit of 15 October 2026; audit_402710280.py). The sector-66
    # window lies on a stray-light ramp (data start 0.4 h before ingress) and is shown but not used; the shape is
    # fitted to the four transits (39, 93, 100, 102), each window centred on its measured mid-time; the ephemeris
    # comes from the four mid-times re-measured with that shape, its errors inflated by sqrt(chi2_red) because the
    # times scatter by tens of minutes about a straight line.
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
    np.save(RESULTS / "candidates" / "shape_chain_fixed.npy", ch.astype(np.float32))
    # re-time each transit with this shape (the windows above were centred on times measured with a preliminary
    # shape, audit_402710280.py); these times give the ephemeris and Table 2 of the paper
    import math
    from audit_402710280 import raw_window, time_one
    by = {lc["meta"]["sector"]: lc for lc in lcs}
    rp_, b_, aR_ = out["rp_rs"][0], out["b"][0], out["a_rs"][0]
    inc_ = math.degrees(math.acos(b_ / aR_))
    pre = meas
    meas = {}
    for s in (39, 93, 100, 102):
        t, f = raw_window(by[s], pre[s]["tc"])
        x = time_one(t, f, pre[s]["tc"], 36.3399, rp_, aR_, inc_, out["u"], float(np.median(np.diff(t)) * 86400), 5.3 / 24)
        x.update(sector=s, product=pre[s]["product"], shift_from_preliminary_min=(x["tc"] - pre[s]["tc"]) * 1440)
        meas[s] = x
    print("re-timed", {s: (round(x["tc"], 5), round(x["err_min"], 2), round(x["shift_from_preliminary_min"], 2))
                       for s, x in meas.items()}, flush=True)
    dep = np.array([meas[s]["depth_ppm"] for s in meas])
    derr = np.array([meas[s]["depth_err_ppm"] for s in meas])
    w = 1 / derr ** 2
    dmean = float(np.sum(w * dep) / w.sum())
    au["per_transit_depth_mean"] = [dmean, float(1 / np.sqrt(w.sum()))]
    out["depth_chi2"] = [float(np.sum(((dep - dmean) / derr) ** 2)), len(dep) - 1]
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
    out["t0_period_cov_d2"] = float(cov[0, 1] * sc ** 2)
    out["ephemeris_chi2"] = [chi2, 2]
    out["times"] = [{"sector": s, "epoch": ep_of[s], "tc": meas[s]["tc"], "err_min": meas[s]["err_min"],
                     "depth_ppm": meas[s]["depth_ppm"], "depth_err_ppm": meas[s]["depth_err_ppm"],
                     "product": meas[s]["product"], "beta": meas[s]["beta"], "err_white_min": meas[s]["err_white_min"],
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
    fresh = json.load(open(FITS))           # other fits may have been written to the file meanwhile
    fresh["TIC402710280"] = fits["TIC402710280"]
    json.dump(fresh, open(FITS, "w"), indent=1, default=float)
    write_extras()


def write_extras():
    """Orbit, derived quantities and upcoming transits (results/report/extras.json) from the stored fit."""
    fresh = json.load(open(FITS))
    ex = json.load(open(EXTRAS))
    ex["TIC402710280"] = predictions.periodic(fresh, "TIC402710280", 0.88, 0.817, 5182.0, 248.524, -63.55, now=4325.5)
    json.dump(ex, open(EXTRAS, "w"), indent=1, default=float)
    print("TIC402710280", {a: b for a, b in ex["TIC402710280"].items() if a not in ("next", "all_next")})

if __name__ == "__main__":
    main()
