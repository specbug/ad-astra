"""TIC 402710280.01: re-time the Sector 100 and 102 transits in the QLP light curves (released at MAST in 2026 March-July
but absent from the AWS mirror that s3io reads, so version 1 of the preprint used only our own full-frame photometry,
ffi_new.py, for these sectors), and re-measure the Sector 66 window.

Timing uses the preprint's method unchanged (audit_402710280.time_one): shape fixed at the published median model,
mid-time scanned on a 0.5-minute grid over +-2 h with a quadratic baseline refitted at every step, 1-sigma error from
delta-chi2 = 1 times the red-noise factor beta. Sector 93 (QLP, already in the preprint) is re-run as a control.

QLP time stamps: each QLP cadence carries the SPOC cadence number (CADENCENO), which is the FFIINDEX keyword of the
full-frame image it came from. For cadences at the transit and 0.11 and 0.5 d either side, the image header gives the
exposure mid-time (TSTART + TSTOP)/2, which equals the TIME of the TESScut cube; corrected to the target's barycentric
position (ffi_new.target_times), it is the time that cadence should carry. QLP minus that time is fitted with a line
across the window and removed before the QLP light curve is timed. The offset is +2.5 s in Sector 93 and about -2 min
in Sectors 100 and 102. Our photometry carries exactly the target-corrected times (checked cadence by cadence).

Control: our own photometry of the Sector 39 and 93 transits, timed the same way, against the published TESS-SPOC
and QLP times.

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


def qlp_raw(sector):
    """CADENCENO and TIME columns of the QLP file (s3io.parse_lc drops CADENCENO)."""
    import io
    from astropy.io import fits
    k = s3io.qlp_key(TIC, sector)
    raw = s3io.fetch(k) or s3io.fetch("https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:HLSP/" + k.split("mast/hlsp/", 1)[1])
    d = fits.open(io.BytesIO(raw))[1].data
    return np.asarray(d["CADENCENO"], int), np.asarray(d["TIME"], float)


def ffi_header(name):
    """Primary header of a calibrated full-frame image (first 20 FITS blocks only)."""
    import requests
    from astropy.io import fits
    sec, cam, ccd = name.split("-")[1], name.split("-")[2], name.split("-")[3]
    url = f"https://stpubdata.s3.amazonaws.com/tess/public/ffi/{sec}/{name[4:8]}/{name[8:11]}/{cam}-{ccd}/{name}"
    r = requests.get(url, headers={"Range": "bytes=0-57599"}, timeout=60)
    r.raise_for_status()
    return fits.Header.fromstring(r.content[: 2880 * (len(r.content) // 2880)].decode("ascii", "ignore"))


CUBE = {93: (2, 4), 100: (4, 1), 102: (2, 2)}


def qlp_time_offset(s, tc):
    """QLP TIME minus the correct target-barycentric mid-exposure time, at cadences across the transit window."""
    import os
    from astropy.io import fits
    import pixels
    from ffi_new import target_times
    ra, dec = 248.5236, -63.5497
    path = pixels.cutout(ra, dec, s, *CUBE[s], size=3)
    with fits.open(path) as h:
        d = h[1].data
        ct, corr, ff = np.asarray(d["TIME"], float), np.asarray(d["TIMECORR"], float), np.asarray(d["FFI_FILE"])
    os.remove(path)
    qn, qt = qlp_raw(s)
    rows = []
    for dt in (-0.5, -0.11, 0.0, 0.11, 0.5):
        i = int(np.nanargmin(np.abs(ct - (tc + dt))))
        name = str(ff[i]).strip()
        name = name if name.endswith(".fits") else name.rstrip("cr") + "c.fits"
        hdr = ffi_header(name)
        mid = 0.5 * (hdr["TSTART"] + hdr["TSTOP"])
        assert abs(mid - ct[i]) * 86400 < 0.5, "cube TIME is not the image mid-exposure time"
        right = float(target_times(np.array([ct[i]]), np.array([corr[i]]), ra, dec)[0])
        j = np.where(qn == hdr["FFIINDEX"])[0]
        if j.size:
            rows.append({"ffiindex": int(hdr["FFIINDEX"]), "t": right, "qlp_minus_true_s": float((qt[j[0]] - right) * 86400)})
    a = np.polyfit([r["t"] - tc for r in rows], [r["qlp_minus_true_s"] for r in rows], 1)
    return {"anchors": rows, "slope_s_per_d": float(a[0]), "offset_at_transit_s": float(a[1])}, a


def retime(d, rp, aR, inc, u):
    pub = {x["sector"]: x for x in d["times"]}
    P, T0, ref = d["period_d"][0], d["tref_btjd"][0], d["ref_epoch"]
    out = {}
    for s in (93, 100, 102):
        lc = load_qlp(s)
        t, f = raw_window(lc, pub[s]["tc"])
        exp_s = float(np.median(np.diff(t)) * 86400)
        raw = time_one(t, f, pub[s]["tc"], P, rp, aR, inc, u, exp_s, D)
        off, a = qlp_time_offset(s, pub[s]["tc"])
        t_fix = t - np.polyval(a, t - pub[s]["tc"]) / 86400
        x = time_one(t_fix, f, pub[s]["tc"], P, rp, aR, inc, u, exp_s, D)
        x.pop("chi2_min", None)
        x["time_offset"] = off
        x["tc_uncorrected"] = raw["tc"]
        x["o_c_min"] = (x["tc"] - T0 - (pub[s]["epoch"] - ref) * P) * 1440
        x["published_tc"], x["published_err_min"], x["published_product"] = pub[s]["tc"], pub[s]["err_min"], pub[s]["product"]
        x["minus_published_min"] = (x["tc"] - pub[s]["tc"]) * 1440
        x["sigma_of_difference"] = x["minus_published_min"] / math.hypot(x["err_min"], pub[s]["err_min"])
        x["depth_minus_published_sigma"] = (x["depth_ppm"] - pub[s]["depth_ppm"]) / math.hypot(x["depth_err_ppm"], pub[s]["depth_err_ppm"])
        out[s] = x
    out["interval_100_102"] = {
        "qlp_minus_P_min": ((out[102]["tc"] - out[100]["tc"]) - P) * 1440,
        "published_minus_P_min": ((pub[102]["tc"] - pub[100]["tc"]) - P) * 1440}
    return out


def own_control(d, rp, aR, inc, u):
    """Our full-frame photometry (ffi_new.py, same settings as Sectors 100 and 102) of the Sector 39 and 93 transits,
    timed as above, against the published TESS-SPOC and QLP times: a test for an offset specific to our reduction."""
    import ffi_new
    pub = {x["sector"]: x for x in d["times"]}
    P = d["period_d"][0]
    out = {}
    for s, cam, ccd in ((39, 2, 4), (93, 2, 4)):
        lc = ffi_new.build(TIC, 248.5236, -63.5497, s, cam, ccd, target_time=True)
        t, f = raw_window(lc, pub[s]["tc"])
        x = time_one(t, f, pub[s]["tc"], P, rp, aR, inc, u, float(np.median(np.diff(t)) * 86400), D)
        out[s] = {"tc": x["tc"], "err_min": x["err_min"], "depth_ppm": x["depth_ppm"], "depth_err_ppm": x["depth_err_ppm"],
                  "published_product": pub[s]["product"], "minus_published_min": (x["tc"] - pub[s]["tc"]) * 1440,
                  "sigma": (x["tc"] - pub[s]["tc"]) * 1440 / math.hypot(x["err_min"], pub[s]["err_min"])}
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
    res["first_quality0_rel_mid_h"] = float((t.min() - tc_pred) * 24)
    x = time_one(t, f, tc_pred, P, rp, aR, inc, u, exp_s, D)
    res["free_time"] = {"tc": x["tc"], "err_min": x["err_min"], "o_c_min": (x["tc"] - tc_pred) * 1440,
                        "depth_ppm": x["depth_ppm"], "depth_err_ppm": x["depth_err_ppm"], "beta": x["beta"]}
    # mid-time and depth with the time free, for each baseline order (time_one uses the quadratic)
    grid = tc_pred + np.arange(-2, 2 + 1e-9, 0.5 / 60) / 24
    res["free_time_by_baseline"] = {}
    for name, deg in (("constant", 0), ("linear", 1), ("quadratic", 2)):
        chi = []
        for tc in grid:
            m = tf._model(t, tc, P, rp, aR, inc, u, exp_s=exp_s)
            A = np.vstack([m * (t - tc) ** k for k in range(deg + 1)]).T
            c, *_ = np.linalg.lstsq(A, f, rcond=None)
            chi.append(float(np.sum(((f - A @ c) / sig) ** 2)))
        chi = np.array(chi)
        i = int(np.argmin(chi))
        lo, hi = max(i - 6, 0), min(i + 7, grid.size)
        a2, a1, _ = np.polyfit(grid[lo:hi] - grid[i], chi[lo:hi], 2)
        tb = grid[i] - a1 / (2 * a2)
        m = tf._model(t, tb, P, rp, aR, inc, u, exp_s=exp_s)
        A = np.vstack([(t - tb) ** k for k in range(deg + 1)] + [-(1 - m)]).T
        c, *_ = np.linalg.lstsq(A, f, rcond=None)
        cov = np.linalg.inv(A.T @ A) * sig ** 2
        dm = 1 - m.min()
        res["free_time_by_baseline"][name] = {
            "o_c_min": float((tb - tc_pred) * 1440), "err_white_min": float(1440 / math.sqrt(a2)),
            "at_grid_edge": bool(i in (0, grid.size - 1)),
            "depth_ppm": float(c[-1] * dm * 1e6), "depth_err_ppm": float(math.sqrt(cov[-1, -1]) * dm * 1e6)}
    return res


def main():
    d, rp, aR, inc, u = shape()
    out = {"retime_qlp": retime(d, rp, aR, inc, u), "own_control": own_control(d, rp, aR, inc, u),
           "s66": s66(d, rp, aR, inc, u)}
    json.dump(out, open(OUT, "w"), indent=1, default=float)
    r = out["retime_qlp"]
    for s in (93, 100, 102):
        x = r[s]
        print(f"S{s} QLP tc={x['tc']:.5f} +- {x['err_min']:.2f} min (beta {x['beta']:.2f}); O-C {x['o_c_min']:+.1f} min; "
              f"minus published ({x['published_product']}) {x['minus_published_min']:+.1f} min "
              f"= {x['sigma_of_difference']:+.1f} sigma; depth {x['depth_ppm']:.0f} +- {x['depth_err_ppm']:.0f} ppm")
        print(f"     QLP time offset at transit {x['time_offset']['offset_at_transit_s']:+.1f} s "
              f"(slope {x['time_offset']['slope_s_per_d']:+.1f} s/d); uncorrected minus published {(x['tc_uncorrected'] - x['published_tc']) * 1440:+.1f} min; "
              f"depth minus published {x['depth_minus_published_sigma']:+.1f} sigma")
    print("S100->S102 interval minus P: QLP %+.1f min, published %+.1f min" % (
        r["interval_100_102"]["qlp_minus_P_min"], r["interval_100_102"]["published_minus_P_min"]))
    for s, x in out["own_control"].items():
        print(f"control S{s}: our photometry minus published {x['published_product']} {x['minus_published_min']:+.1f} min ({x['sigma']:+.1f} sigma)")
    print("S66:", json.dumps(out["s66"], indent=1, default=float))


if __name__ == "__main__":
    main()
