"""TIC 402710280.01: re-time the Sector 100 and 102 transits in the QLP light curves released after the preprint
(the preprint timed them in our own full-frame photometry, ffi_new.py), and re-measure the Sector 66 window.

Timing uses the preprint's method unchanged (audit_402710280.time_one): shape fixed at the published median model,
mid-time scanned on a 0.5-minute grid over +-2 h with a quadratic baseline refitted at every step, 1-sigma error from
delta-chi2 = 1 times the red-noise factor beta. Sector 93 (QLP, already in the preprint) is re-run as a control.

QLP time stamps: the high-passed QLP flux is cross-correlated with our full-frame photometry (whose times match
each source FFI's file time to within 1 s) at whole-cadence lags, orbit by orbit. A peak at lag +1 means QLP labels
the flux one 200-s cadence early; the QLP mid-times are then corrected by that lag. Sectors 39 and 93 peak at lag 0.

Sector 66: TESS-SPOC PDCSAP, quality-flag-zero cadences only, around the predicted mid-time; the depth is fitted at
the predicted time and at the best time, with three baselines, to show how much the one-sided baseline moves it.

  python qlp_recheck_402710280.py
Output: results/candidates/qlp_recheck_402710280.json
"""
import json
import math

import env  # noqa: F401
import numpy as np

import s3io
import transit_fit as tf
from audit_402710280 import raw_window, time_one
from env import RESULTS

TIC = 402710280
FITS = RESULTS / "report" / "transit_fits.json"
OUT = RESULTS / "candidates" / "qlp_recheck_402710280.json"
D = 5.3 / 24


def shape():
    d = json.load(open(FITS))["TIC402710280"]
    rp, aR, b = d["rp_rs"][0], d["a_rs"][0], d["b"][0]
    return d, rp, aR, math.degrees(math.acos(b / aR)), d["u"]


def load_qlp(sector):
    """QLP light curve from the AWS mirror, or from MAST directly for sectors the mirror does not hold yet."""
    lc = s3io.load("qlp", TIC, sector)
    if lc is None:
        uri = "mast:HLSP/" + s3io.qlp_key(TIC, sector).split("mast/hlsp/", 1)[1]
        raw = s3io.fetch(f"https://mast.stsci.edu/api/v0.1/Download/file?uri={uri}")
        lc = s3io.parse_lc(raw, "qlp")
    return lc


def cadence_lag(s, tc, cad=200 / 86400):
    """Whole-cadence lag of QLP against our photometry in the orbit holding the transit at tc (+1: QLP early)."""
    z = np.load(f"{env.WORK}/cache/ffi_lc/{TIC}_s{s}.npz", allow_pickle=True)
    q = load_qlp(s)
    gap = q["time"][np.argmax(np.diff(q["time"]))]
    sel = (lambda t: t > gap) if tc > gap else (lambda t: t <= gap)
    t0 = q["time"][0]

    def grid(t, f):
        k = np.round((t - t0) / cad).astype(int)
        g = np.full(k.max() + 10, np.nan)
        g[k[k >= 0]] = f[k >= 0]
        m = np.isfinite(g)
        sm = np.convolve(np.where(m, g, 0), np.ones(31), "same") / np.maximum(np.convolve(m.astype(float), np.ones(31), "same"), 1)
        return g - sm

    a = grid(z["time"][sel(z["time"])], z["sap_flux"][sel(z["time"])])
    b = grid(q["time"][sel(q["time"])], q["sap_flux"][sel(q["time"])])
    n = min(a.size, b.size)
    r = {}
    for lag in range(-2, 3):
        x, y = a[max(0, lag):n + min(0, lag)], b[max(0, -lag):n - max(0, lag)]
        m = np.isfinite(x) & np.isfinite(y)
        r[lag] = float(np.corrcoef(x[m], y[m])[0, 1])
    return max(r, key=r.get), r


def retime(d, rp, aR, inc, u):
    pub = {x["sector"]: x for x in d["times"]}
    P, T0, ref = d["period_d"][0], d["tref_btjd"][0], d["ref_epoch"]
    out = {}
    for s in (93, 100, 102):
        lc = load_qlp(s)
        t, f = raw_window(lc, pub[s]["tc"])
        x = time_one(t, f, pub[s]["tc"], P, rp, aR, inc, u, float(np.median(np.diff(t)) * 86400), D)
        x.pop("chi2_min", None)
        x["o_c_min"] = (x["tc"] - T0 - (pub[s]["epoch"] - ref) * P) * 1440
        x["published_tc"], x["published_err_min"], x["published_product"] = pub[s]["tc"], pub[s]["err_min"], pub[s]["product"]
        x["minus_published_min"] = (x["tc"] - pub[s]["tc"]) * 1440
        x["sigma_of_difference"] = x["minus_published_min"] / math.hypot(x["err_min"], pub[s]["err_min"])
        if s in (100, 102):
            lag, corr = cadence_lag(s, pub[s]["tc"])
            x["qlp_cadence_lag"], x["lag_correlations"] = lag, corr
            x["tc_lag_corrected"] = x["tc"] + lag * 200 / 86400
            x["lag_corrected_minus_published_min"] = (x["tc_lag_corrected"] - pub[s]["tc"]) * 1440
            x["lag_corrected_sigma"] = x["lag_corrected_minus_published_min"] / math.hypot(x["err_min"], pub[s]["err_min"])
        out[s] = x
    out["interval_100_102"] = {
        "qlp_minus_P_min": ((out[102]["tc"] - out[100]["tc"]) - P) * 1440,
        "published_minus_P_min": ((pub[102]["tc"] - pub[100]["tc"]) - P) * 1440}
    return out


def s66(d, rp, aR, inc, u):
    from astropy.io import fits
    import io
    P, T0, ref = d["period_d"][0], d["tref_btjd"][0], d["ref_epoch"]
    tc_pred = T0 + (20 - ref) * P
    h = fits.open(io.BytesIO(s3io.fetch(s3io.key_for("tspoc", TIC, 66))))[1].data
    t, q, f = (np.asarray(h[k], float) for k in ("TIME", "QUALITY", "PDCSAP_FLUX"))
    ok = np.isfinite(t) & np.isfinite(f) & (q == 0) & (np.abs(t - tc_pred) < 0.8)
    t, f = t[ok], f[ok] / np.nanmedian(f[ok])
    exp_s = float(np.median(np.diff(t)) * 86400)
    res = {"tc_pred_btjd": tc_pred, "n_quality0": int(ok.sum()),
           "first_quality0_rel_ingress_h": float((t.min() - (tc_pred - 5.27 / 48)) * 24)}
    m = tf._model(t, tc_pred, P, rp, aR, inc, u, exp_s=exp_s)
    depth0 = 1 - m.min()
    oot = np.abs(t - tc_pred) > 0.75 * D
    sig = 1.4826 * np.median(np.abs(np.diff(f[oot]))) / math.sqrt(2)
    fits_ = {}
    for name, deg in (("constant", 0), ("linear", 1), ("quadratic", 2)):
        dt = t - tc_pred
        A = np.vstack([dt ** k for k in range(deg + 1)] + [-(1 - m)]).T
        c, *_ = np.linalg.lstsq(A, f, rcond=None)
        cov = np.linalg.inv(A.T @ A) * sig ** 2
        fits_[name] = {"depth_ppm": float(c[-1] * depth0 * 1e6), "err_ppm": float(math.sqrt(cov[-1, -1]) * depth0 * 1e6)}
    res["depth_at_predicted_time"] = fits_
    intr = np.abs(t - tc_pred) < 0.4 * D
    post = (t - tc_pred > 0.75 * D) & (t - tc_pred < 0.75 * D + 0.25)
    res["box_in_minus_post_ppm"] = float((np.median(f[post]) - np.median(f[intr])) * 1e6)
    x = time_one(t, f, tc_pred, P, rp, aR, inc, u, exp_s, D)
    res["free_time"] = {"tc": x["tc"], "err_min": x["err_min"], "o_c_min": (x["tc"] - tc_pred) * 1440,
                        "depth_ppm": x["depth_ppm"], "depth_err_ppm": x["depth_err_ppm"], "beta": x["beta"]}
    return res


def main():
    d, rp, aR, inc, u = shape()
    out = {"retime_qlp": retime(d, rp, aR, inc, u), "s66": s66(d, rp, aR, inc, u)}
    json.dump(out, open(OUT, "w"), indent=1, default=float)
    r = out["retime_qlp"]
    for s in (93, 100, 102):
        x = r[s]
        print(f"S{s} QLP tc={x['tc']:.5f} +- {x['err_min']:.2f} min (beta {x['beta']:.2f}); O-C {x['o_c_min']:+.1f} min; "
              f"minus published ({x['published_product']}) {x['minus_published_min']:+.1f} min "
              f"= {x['sigma_of_difference']:+.1f} sigma; depth {x['depth_ppm']:.0f} +- {x['depth_err_ppm']:.0f} ppm")
        if "qlp_cadence_lag" in x:
            print(f"     QLP cadence lag {x['qlp_cadence_lag']:+d} (r = {x['lag_correlations'][x['qlp_cadence_lag']]:.2f}); "
                  f"lag-corrected minus published {x['lag_corrected_minus_published_min']:+.1f} min = {x['lag_corrected_sigma']:+.1f} sigma")
    print("S100->S102 interval minus P: QLP %+.1f min, published %+.1f min" % (
        r["interval_100_102"]["qlp_minus_P_min"], r["interval_100_102"]["published_minus_P_min"]))
    print("S66:", json.dumps(out["s66"], indent=1, default=float))


if __name__ == "__main__":
    main()
