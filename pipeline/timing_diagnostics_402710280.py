"""Diagnostics behind the timing scatter of TIC 402710280.01, and the star's variability.

  1. each mid-time re-measured with constant, linear and quadratic baselines and window half-widths 0.6, 0.8, 1.0 d;
  2. the sector-93 transit (QLP) timed on each of its flux columns (detrended, raw aperture, small and large
     apertures), with the background and centroid traces across the window;
  3. Lomb-Scargle periodogram of every sector (rotation), and a flare count (>= 3 consecutive cadences > 3 sigma);
  4. sector-103 coverage around the predicted epoch-49 transit.

  python timing_diagnostics_402710280.py
Output: results/candidates/timing_diagnostics_402710280.json
"""
import json
import math

import env  # noqa: F401
import numpy as np
from astropy.timeseries import LombScargle

import transit_fit as tf
from audit_402710280 import beta_factor
from env import RESULTS
from lcaccess import lcs_for

TIC = 402710280
P, T0, REF = 36.33989590709609, 2460838.20190249 - 2457000, 40
RP, B, AR = 0.04208122638534462, 0.21934461732459248, 53.68913012236519
U = (0.39272, 0.20682)
INC = math.degrees(math.acos(B / AR))
D = 5.2723 / 24
TCS = {39: 2384.6094, 93: 3838.1847, 100: 4092.5851, 102: 4128.9449}
EPOCH = {39: 0, 93: 40, 100: 47, 102: 48}
OUT = RESULTS / "candidates" / "timing_diagnostics_402710280.json"


def time_order(t, f, tc0, exp_s, order, span_h=2.0):
    grid = tc0 + np.arange(-span_h, span_h + 1e-9, 0.5 / 60) / 24
    oot = np.abs(t - tc0) > 0.75 * D
    sig = 1.4826 * np.median(np.abs(np.diff(f[oot]))) / math.sqrt(2)
    chi = []
    for tc in grid:
        m = tf._model(t, tc, P, RP, AR, INC, U, exp_s=exp_s)
        A = np.vstack([m * (t - tc) ** j for j in range(order + 1)]).T
        c, *_ = np.linalg.lstsq(A, f, rcond=None)
        chi.append(float(np.sum(((f - A @ c) / sig) ** 2)))
    chi = np.array(chi)
    i = int(np.argmin(chi))
    lo, hi = max(i - 6, 0), min(i + 7, len(grid))
    a2, a1, _ = np.polyfit(grid[lo:hi] - grid[i], chi[lo:hi], 2)
    tb = grid[i] - a1 / (2 * a2) if a2 > 0 else grid[i]
    m = tf._model(t, tb, P, RP, AR, INC, U, exp_s=exp_s)
    A = np.vstack([m * (t - tb) ** j for j in range(order + 1)]).T
    c, *_ = np.linalg.lstsq(A, f, rcond=None)
    beta = beta_factor(t, f - A @ c, sig)
    return float(tb), float(beta / math.sqrt(a2) * 1440) if a2 > 0 else float("nan"), float(beta)


def window(lc, tc, win, col="flux"):
    t, f = lc["time"], lc[col]
    ok = np.isfinite(t) & np.isfinite(f) & (np.abs(t - tc) < win)
    return t[ok], f[ok] / np.nanmedian(f[ok]), ok


def main():
    lcs = lcs_for(TIC)
    by = {lc["meta"]["sector"]: lc for lc in lcs}
    out = {"baseline_and_window": {}, "s93_columns": {}, "variability": {}}
    for s, tc in TCS.items():
        lc = by[s]
        exp_s = float(lc["meta"]["cadence_s"])
        pred = T0 + (EPOCH[s] - REF) * P
        rows = []
        for win in (0.6, 0.8, 1.0):
            t, f, _ = window(lc, tc, win)
            for order in (0, 1, 2):
                tb, err, beta = time_order(t, f, tc, exp_s, order)
                rows.append({"win_d": win, "order": order, "tc": tb, "err_min": err, "beta": beta,
                             "o_minus_c_min": (tb - pred) * 1440})
        oc = [r["o_minus_c_min"] for r in rows]
        out["baseline_and_window"][s] = {"rows": rows, "o_minus_c_range_min": [min(oc), max(oc)]}
        print(s, [round(x, 1) for x in oc], flush=True)
    # sector 93 flux columns, background and centroid in the window
    lc = by[93]
    tc = TCS[93]
    pred = T0 + (EPOCH[93] - REF) * P
    for col in ("flux", "sap_flux", "det_flux", "flux_sml", "flux_lag"):
        if col not in lc:
            continue
        t, f, _ = window(lc, tc, 0.8, col)
        if t.size < 50:
            continue
        tb, err, beta = time_order(t, f, tc, 200.0, 2)
        out["s93_columns"][col] = {"tc": tb, "err_min": err, "o_minus_c_min": (tb - pred) * 1440}
    t, _, ok = window(lc, tc, 0.8)
    inn = np.abs(t - tc) < D / 2
    for col in ("bkg", "cx", "cy"):
        v = np.asarray(lc[col], float)[ok]
        out["s93_columns"][col + "_in_minus_out"] = float(np.nanmedian(v[inn]) - np.nanmedian(v[~inn]))
        out["s93_columns"][col + "_oot_scatter"] = float(np.nanstd(v[~inn]))
    dt = np.diff(t)
    out["s93_columns"]["max_gap_in_window_min"] = float(dt.max() * 1440)
    out["s93_columns"]["n_cadences_window"] = int(t.size)
    out["s93_columns"]["expected_cadences_window"] = int(1.6 * 86400 / 200)
    out["s93_columns"]["orbit_boundary_d"] = float(np.min(np.abs(lc["time"][np.r_[False, np.diff(lc["orbit"]) != 0]] - tc))) \
        if "orbit" in lc and np.any(np.diff(lc["orbit"]) != 0) else None
    print("s93", out["s93_columns"], flush=True)
    # variability
    for s, lc in sorted(by.items()):
        t, f = lc["time"], lc["flux"]
        ok = np.isfinite(t) & np.isfinite(f)
        t, f = t[ok], f[ok] / np.nanmedian(f[ok])
        mask = np.abs(((t - T0 + P / 2) % P) - P / 2) > 0.3
        t, f = t[mask], f[mask]
        fr = np.linspace(1 / 13, 2, 20000)
        pw = LombScargle(t, f).power(fr)
        i = int(np.argmax(pw))
        ph = (t * fr[i]) % 1
        amp = float(np.ptp([np.median(f[(ph >= a) & (ph < a + 0.1)]) for a in np.arange(0, 1, 0.1)]) * 1e6)
        r = f - np.convolve(f, np.ones(15) / 15, mode="same")
        sig = 1.4826 * np.median(np.abs(r))
        hi = r > 3 * sig
        runs = np.diff(np.r_[0, hi.astype(int), 0])
        starts, ends = np.where(runs == 1)[0], np.where(runs == -1)[0]
        n_fl = int(np.sum((ends - starts) >= 3))
        out["variability"][s] = {"best_period_d": float(1 / fr[i]), "power": float(pw[i]), "fap": float(LombScargle(t, f).false_alarm_probability(pw[i])),
                                 "peak_to_peak_ppm": amp, "rms_ppm": float(np.std(f) * 1e6), "n_flare_candidates": n_fl,
                                 "product": lc["meta"]["product"]}
        print(s, out["variability"][s], flush=True)
    # sector 103 around epoch 49
    lc = by[103]
    t = lc["time"][np.isfinite(lc["flux"])]
    pred49 = T0 + (49 - REF) * P
    gaps = [(float(a), float(b)) for a, b in zip(t[:-1], t[1:]) if b - a > 0.1]
    out["s103"] = {"data_start": float(t.min()), "data_end": float(t.max()), "gaps": gaps, "pred_epoch49": pred49,
                   "n_points_within_0.5d": int(np.sum(np.abs(t - pred49) < 0.5))}
    print("s103", out["s103"], flush=True)
    json.dump(out, open(OUT, "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
