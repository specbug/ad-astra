"""TIC 402710280.01 after the independent audit (15 October 2026): the sector-66 window, individual transit times
with red-noise errors, the ephemeris they imply, and a transit fit without sector 66.

  - S66: the raw TESS-SPOC and QLP files around the predicted time (quality flags, SAP and PDCSAP flux, background);
  - transit times: for each transit the shape is fixed at the joint fit and the mid-time is scanned on a grid, the
    window's quadratic baseline refitted at every step (so the baseline cannot be set by the transit itself); the
    1-sigma error from delta-chi2 = 1, multiplied by the red-noise factor beta (residual rms in 10-60 minute bins
    against the white-noise expectation, Winn et al. 2008);
  - ephemeris: weighted straight line through the times, errors scaled by sqrt(chi2_red) when above 1;
  - joint fits (transit_fit.fit_periodic_windows) on S39 + S93 + S100 + S102, without S102, and on S39 + S93 alone
    (the public light curves), each with per-window errors scaled by beta.

  python audit_402710280.py
Output: results/candidates/audit_402710280.json, results/report/cases/TIC402710280_audit_s66.png
"""
import io
import json
import math

import env  # noqa: F401
import numpy as np

import transit_fit as tf
from env import RESULTS

TIC = 402710280
P0, T00 = 36.33972, 3474.8014
TRANSITS = {39: 2384.607, 93: 3838.214, 100: 4092.595, 102: 4128.936}
S66_TC = 3111.411
TEFF = 5180.0
OUT = RESULTS / "candidates" / "audit_402710280.json"
FIG = RESULTS / "report" / "cases" / "TIC402710280_audit_s66.png"


def beta_factor(t, r, sig):
    """Red-noise factor: rms of residuals binned over 10-60 minutes against sig/sqrt(N) (N points per bin)."""
    cad = float(np.median(np.diff(t)))
    betas = []
    for m in (10, 20, 30, 45, 60):
        n = max(int(round(m / 1440 / cad)), 1)
        if n < 2:
            continue
        k = len(r) // n
        if k < 4:
            continue
        rb = r[: k * n].reshape(k, n).mean(axis=1)
        betas.append(float(np.std(rb) / (sig / math.sqrt(n)) * math.sqrt(k / (k - 1))))
    return max(1.0, float(np.median(betas))) if betas else 1.0


def raw_window(lc, tc, win=0.8):
    t, f = lc["time"], lc["flux"]
    ok = np.isfinite(t) & np.isfinite(f) & (np.abs(t - tc) < win)
    return t[ok], f[ok] / np.nanmedian(f[ok])


def time_one(t, f, tc0, P, rp, aR, inc, u, exp_s, D, span_h=2.0):
    """Mid-time on a grid with a quadratic baseline refitted at each step; returns tc, err (white x beta), beta,
    depth scale, and the chi2 profile."""
    grid = tc0 + np.arange(-span_h, span_h + 1e-9, 0.5 / 60) / 24
    oot = np.abs(t - tc0) > 0.75 * D
    sig = 1.4826 * np.median(np.abs(np.diff(f[oot]))) / math.sqrt(2)
    chi = []
    for tc in grid:
        m = tf._model(t, tc, P, rp, aR, inc, u, exp_s=exp_s)
        dt = t - tc
        A = np.vstack([m, m * dt, m * dt ** 2]).T
        c, *_ = np.linalg.lstsq(A, f, rcond=None)
        chi.append(float(np.sum(((f - A @ c) / sig) ** 2)))
    chi = np.array(chi)
    i = int(np.argmin(chi))
    # parabola through the minimum for sub-grid precision, and delta-chi2 = 1 width
    lo, hi = max(i - 6, 0), min(i + 7, len(grid))
    a2, a1, a0 = np.polyfit(grid[lo:hi] - grid[i], chi[lo:hi], 2)
    tbest = grid[i] - a1 / (2 * a2) if a2 > 0 else grid[i]
    err_w = 1 / math.sqrt(a2) if a2 > 0 else float("nan")
    m = tf._model(t, tbest, P, rp, aR, inc, u, exp_s=exp_s)
    dt = t - tbest
    # depth scale: free amplitude of the transit (1 - m) with the baseline
    A = np.vstack([np.ones_like(t), dt, dt ** 2, -(1 - m)]).T
    c, *_ = np.linalg.lstsq(A, f, rcond=None)
    resid = f - A @ c
    beta = beta_factor(t, resid, sig)
    cov = np.linalg.inv(A.T @ A) * sig ** 2 * beta ** 2
    model_depth = float(1 - m.min())
    return {"tc": float(tbest), "err_white_min": err_w * 1440, "beta": beta, "err_min": err_w * beta * 1440,
            "sigma_ppm": sig * 1e6, "depth_scale": float(c[3]), "depth_scale_err": float(math.sqrt(cov[3, 3])),
            "depth_ppm": float(c[3] * model_depth * 1e6), "depth_err_ppm": float(math.sqrt(cov[3, 3]) * model_depth * 1e6),
            "chi2_min": float(chi[i]), "n": int(len(t)), "dchi2_no_transit": float(np.sum(((f - np.vstack(
                [np.ones_like(t), t - tbest, (t - tbest) ** 2]).T @ np.linalg.lstsq(np.vstack(
                    [np.ones_like(t), t - tbest, (t - tbest) ** 2]).T, f, rcond=None)[0]) / sig) ** 2) - chi[i])}


def ephemeris(times):
    n = np.array([x["epoch"] for x in times], float)
    tt = np.array([x["tc"] for x in times])
    st = np.array([x["err_min"] for x in times]) / 1440
    A = np.vstack([np.ones_like(n), n]).T
    W = np.diag(1 / st ** 2)
    cov = np.linalg.inv(A.T @ W @ A)
    p = cov @ A.T @ W @ tt
    res = tt - A @ p
    chi2 = float(np.sum((res / st) ** 2))
    dof = max(len(tt) - 2, 1)
    s = math.sqrt(max(chi2 / dof, 1.0))
    return {"T0_btjd": float(p[0]), "T0_err_min": float(math.sqrt(cov[0, 0]) * s * 1440), "P": float(p[1]),
            "P_err": float(math.sqrt(cov[1, 1]) * s), "chi2": chi2, "dof": dof,
            "residuals_min": {int(x["sector"]): float(r * 1440) for x, r in zip(times, res)},
            "residuals_sigma": {int(x["sector"]): float(r / e) for x, r, e in zip(times, res, st)}}


def s66_raw():
    import s3io
    from astropy.io import fits
    out = {}
    for prod in ("tspoc", "qlp"):
        d = fits.open(io.BytesIO(s3io.fetch(s3io.key_for(prod, TIC, 66))))[1].data
        t, q = np.asarray(d["TIME"], float), np.asarray(d["QUALITY"])
        w = np.isfinite(t) & (np.abs(t - S66_TC) < 0.5)
        good = w & (q == 0)
        out[prod] = {"n_within_12h": int(w.sum()), "n_quality0_within_12h": int(good.sum()),
                     "first_quality0_h": float((t[good].min() - S66_TC) * 24) if good.any() else None,
                     "flags": {int(k): int(v) for k, v in zip(*np.unique(q[w], return_counts=True))}}
        if prod == "tspoc":
            out[prod]["arrays"] = {k: np.asarray(d[k], float)[np.isfinite(t) & (np.abs(t - S66_TC) < 0.75)].tolist()
                                   for k in ("TIME", "SAP_FLUX", "PDCSAP_FLUX", "SAP_BKG", "QUALITY")}
    return out


def plot_s66(s66):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    a = {k: np.array(v) for k, v in s66["tspoc"]["arrays"].items()}
    h = (a["TIME"] - S66_TC) * 24
    q0 = a["QUALITY"] == 0
    fig, ax = plt.subplots(3, 1, figsize=(8, 7), sharex=True)
    sap = a["SAP_FLUX"] / np.nanmedian(a["SAP_FLUX"][q0 & (h > 10)])
    pdc = a["PDCSAP_FLUX"] / np.nanmedian(a["PDCSAP_FLUX"][q0 & (h > 10)])
    ax[0].plot(h[~q0], (sap[~q0] - 1) * 1e6, ".", ms=2, color="0.7", label="flagged (stray light / manual exclude)")
    ax[0].plot(h[q0], (sap[q0] - 1) * 1e6, ".", ms=2, color="k", label="quality 0")
    ax[0].set_ylabel("SAP (ppm)")
    ax[0].legend(fontsize=7, loc="lower left")
    ax[1].plot(h[q0], (pdc[q0] - 1) * 1e6, ".", ms=2, color="k")
    ax[1].set_ylabel("PDCSAP (ppm)")
    ax[2].plot(h, a["SAP_BKG"], "-", color="#b5402a")
    ax[2].set_ylabel("background (e-/s)")
    ax[2].set_xlabel("hours from predicted mid-transit (BTJD 3111.411)")
    for x in ax:
        x.axvspan(-2.6, 2.6, color="#2a5db5", alpha=0.08)
    ax[0].set_title("TIC 402710280, sector 66 (TESS-SPOC): the predicted transit falls on a stray-light ramp", fontsize=9)
    fig.tight_layout()
    fig.savefig(FIG, dpi=110)
    plt.close(fig)


def centroid_no_s66():
    """Difference-image centroids of the four used transits (new_candidates.json, eb_pixels), combined, and compared
    with every Gaia source able to make the depth, all moved to the mean epoch of the transits."""
    import ticdb
    from centroids import centroid_verdict, combine_centroids, gaia_field, neighbour_offsets
    d = json.load(open(RESULTS / "candidates" / "new_candidates.json"))
    ev = [dict(e) for e in d[str(TIC)]["centroid"]["per_event"] if e.get("sector") != 66 and "dra" in e]
    st = ticdb.query([TIC], ["ra", "dec", "Tmag"]).iloc[0]
    gf = gaia_field(float(st.ra), float(st.dec), float(st.Tmag))
    epoch = float(np.mean([2016.0 + (e["tc"] + 2457000 - 2457389.0) / 365.25 for e in ev]))
    cen = combine_centroids(ev)
    src = neighbour_offsets(gf, 0.0021, 150.0, epoch=epoch)
    return {"centroid": cen, **centroid_verdict(cen, src), "n_candidates": len(src), "epoch": epoch}


def main():
    from lcaccess import lcs_for
    lcs = [lc for lc in lcs_for(TIC) if lc["time"].size > 100]
    by_sec = {lc["meta"]["sector"]: lc for lc in lcs}
    out = {"centroid_noS66": centroid_no_s66()}
    out |= {"products": {int(s): by_sec[s]["meta"].get("product") for s in (39, 66, 93, 100, 102) if s in by_sec}}
    out["s66_raw"] = s66 = s66_raw()
    plot_s66(s66)
    for s in s66:
        s66[s].pop("arrays", None)
    epochs = {s: int(round((tc - T00) / P0)) for s, tc in TRANSITS.items()}
    D = 5.2 / 24
    # joint fits with beta-scaled errors
    fits = {}
    for name, secs in (("all4", [39, 93, 100, 102]), ("no_S102", [39, 93, 100]), ("public_S39_S93", [39, 93])):
        data = tf.windows_from_lcs([by_sec[s] for s in secs], [TRANSITS[s] for s in secs], 0.8, 5.2)
        scaled = []
        for (t, f, sig, exp_s, tc) in data:
            oot = np.abs(t - tc) > 0.75 * D
            b = beta_factor(t[oot], f[oot] - 1, sig)
            scaled.append((t, f, sig * b, exp_s, tc))
        r, _ = tf.fit_periodic_windows(scaled, [epochs[s] for s in secs], P0, TEFF, nstep=4000)
        r.pop("plot", None)
        fits[name] = r
        print(name, "depth %.0f" % r["depth_ppm"][0], "rp %.4f" % r["rp_rs"][0], "rho %.2f" % r["rho_solar"][0],
              "b %.2f" % r["b"][0], "P %.5f" % r["period_d"][0], flush=True)
    out["fits"] = fits
    # individual times with the all-four shape
    f4 = fits["all4"]
    rp, aR, b = f4["rp_rs"][0], f4["a_rs"][0], f4["b"][0]
    inc = math.degrees(math.acos(b / aR))
    u = f4["u"]
    times = []
    for s, tc in TRANSITS.items():
        t, f = raw_window(by_sec[s], tc)
        exp_s = float(np.median(np.diff(t)) * 86400)
        x = time_one(t, f, tc, P0, rp, aR, inc, u, exp_s, D)
        x.update(sector=s, epoch=epochs[s], product=by_sec[s]["meta"].get("product"))
        times.append(x)
        print(s, "tc %.4f" % x["tc"], "err %.1f min (white %.1f, beta %.2f)" % (x["err_min"], x["err_white_min"], x["beta"]),
              "depth %.0f +- %.0f" % (x["depth_ppm"], x["depth_err_ppm"]), "dchi2 %.0f" % x["dchi2_no_transit"], flush=True)
    out["times"] = times
    out["ephemeris_all4"] = eph = ephemeris(times)
    out["ephemeris_no_S102"] = ephemeris([x for x in times if x["sector"] != 102])
    print("ephemeris", json.dumps(eph, indent=1))
    # S66 with the same machinery, quality-0 TESS-SPOC data only, for the record
    t, f = raw_window(by_sec[66], S66_TC) if 66 in by_sec else (None, None)
    if t is not None and len(t) > 50:
        k = round((S66_TC - eph["T0_btjd"]) / eph["P"])
        tc_pred = eph["T0_btjd"] + k * eph["P"]
        x = time_one(t, f, tc_pred, P0, rp, aR, inc, u, float(np.median(np.diff(t)) * 86400), D)
        x.update(tc_predicted=tc_pred, first_point_h=float((t.min() - tc_pred) * 24), last_point_h=float((t.max() - tc_pred) * 24))
        out["s66_fit"] = x
        print("S66", json.dumps({k: x[k] for k in ("tc", "err_min", "depth_ppm", "depth_err_ppm", "first_point_h", "dchi2_no_transit")}))
    # times again with the shape of the fit without sector 102 (whose window sits on a ramp and distorts the shape);
    # these centre the windows of the final stacked fit (fit_new.py), which re-times the transits with its own shape
    f3 = fits["no_S102"]
    rp, aR, b, u = f3["rp_rs"][0], f3["a_rs"][0], f3["b"][0], f3["u"]
    inc = math.degrees(math.acos(b / aR))
    times = []
    for s, tc in TRANSITS.items():
        t, f = raw_window(by_sec[s], tc)
        x = time_one(t, f, tc, P0, rp, aR, inc, u, float(np.median(np.diff(t)) * 86400), D)
        x.update(sector=s, epoch=epochs[s], product=by_sec[s]["meta"].get("product"))
        times.append(x)
    out["times_noS102_shape"] = times
    out["ephemeris_noS102_shape_all4"] = ephemeris(times)
    dep = np.array([x["depth_ppm"] for x in times])
    err = np.array([x["depth_err_ppm"] for x in times])
    w = 1 / err ** 2
    out["per_transit_depth_mean"] = [float(np.sum(w * dep) / w.sum()), float(1 / math.sqrt(w.sum()))]
    print("times (no-S102 shape)", [(x["sector"], round(x["tc"], 5), round(x["err_min"], 1)) for x in times], flush=True)
    json.dump(out, open(OUT, "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
