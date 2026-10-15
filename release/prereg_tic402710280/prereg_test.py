"""Pre-registered test of TIC 402710280.01 at a predicted transit (frozen; see PROTOCOL.pdf).

Input: a light curve covering one predicted transit, as CSV with columns bjd_tdb, flux (any normalisation), and
optionally other columns, which are ignored; only cadences with valid (finite) flux are used, so quality-flagged cadences must be removed
before the file is given to this script. Output: the outcome (recovered / inconsistent / excluded /
insufficient data), the fitted mid-time, O-C, depth and detection statistic, as JSON on stdout.

  python prereg_test.py lightcurve.csv EPOCH [exposure_seconds]

Nothing in this file may be changed after registration; the SHA-256 checksum in SHA256SUMS.txt identifies it.
"""
import csv
import json
import math
import sys

import batman
import numpy as np

# ---- registered quantities (from the preprint, Table 3; do not edit) ----------
P = 36.33986263847596            # d
T0 = 2460838.2008036566        # BJD_TDB, epoch 40 (the Sector 39 transit is epoch 0)
SIG_P = 0.0005105373388072621    # d
SIG_T0 = 0.011922049057609426    # d
SIG_RMS = 21.77 / 1440           # d, rms of the four timing residuals
REF_EPOCH = 40
RP_RS, B, A_RS = 0.04209035139686697, 0.2243329623513426, 53.629504154358976
U = (0.39272, 0.20682)
T14 = 5.272711765676274 / 24    # d
DEPTH = 2092.0                   # ppm, expected depth at mid-transit
DEPTH_LO, DEPTH_HI = 0.5 * DEPTH, 1.5 * DEPTH
DCHI2_DETECT = 25.0
WINDOW = 0.8                     # d, half-width of the analysis window around the prediction


def predicted(epoch):
    n = epoch - REF_EPOCH
    tc = T0 + n * P
    sig = math.sqrt(SIG_T0 ** 2 + (n * SIG_P) ** 2 + SIG_RMS ** 2)
    return tc, sig


def model(t, tc, exp_s):
    pm = batman.TransitParams()
    pm.t0, pm.per, pm.rp, pm.a = tc, P, RP_RS, A_RS
    pm.inc = math.degrees(math.acos(B / A_RS))
    pm.ecc, pm.w, pm.u, pm.limb_dark = 0.0, 90.0, list(U), "quadratic"
    m = batman.TransitModel(pm, t, supersample_factor=5, exp_time=exp_s / 86400).light_curve(pm)
    return m


def fit_at(t, f, w, tc, exp_s):
    """Weighted linear fit of flux = baseline(quadratic) - amplitude * (1 - model); returns chi2, amplitude, error."""
    dt = t - tc
    shape = 1 - model(t, tc, exp_s)
    if np.count_nonzero(shape > 0) < 3:          # transit falls in a data gap at this trial time: no constraint
        return float("nan"), float("nan"), float("nan")
    A = np.vstack([np.ones_like(t), dt, dt ** 2, -shape]).T
    Aw = A * w[:, None]
    c, *_ = np.linalg.lstsq(Aw, f * w, rcond=None)
    r = (f - A @ c) * w
    cov = np.linalg.inv(Aw.T @ Aw)
    depth_scale = 1 - model(np.array([tc]), tc, exp_s)[0]   # model depth at mid-transit for unit amplitude
    return float(r @ r), float(c[3] * depth_scale * 1e6), float(math.sqrt(cov[3, 3]) * depth_scale * 1e6)


def main(path, epoch, exp_s=120.0):
    rows = list(csv.DictReader(open(path)))
    t = np.array([float(r["bjd_tdb"]) for r in rows])
    f = np.array([float(r["flux"]) for r in rows])
    ok = np.isfinite(t) & np.isfinite(f) & (f > 0)
    tc0, sig = predicted(epoch)
    ok &= np.abs(t - tc0) < WINDOW
    t, f = t[ok], f[ok]
    out = {"epoch": epoch, "predicted_tc": tc0, "sigma_pred_min": sig * 1440, "n_points": int(t.size)}
    half = 3 * sig + T14 / 2
    span = t[np.abs(t - tc0) <= half]
    cad = float(np.median(np.diff(t))) if t.size > 1 else exp_s / 86400
    coverage = min(1.0, span.size * cad / (2 * half)) if t.size else 0.0
    out["coverage_of_window"] = coverage
    if t.size < 30:
        out["outcome"] = "insufficient data"
        return out
    f = f / np.median(f)
    oot = np.abs(t - tc0) > half
    base = f[oot] if oot.sum() > 20 else f
    sig_pt = 1.4826 * np.median(np.abs(np.diff(base))) / math.sqrt(2)
    # registered choice: uniform weights from the out-of-window point-to-point scatter (pipeline errors, if given,
    # are not used, so that all data sources are treated alike)
    w = np.full(t.size, 1 / sig_pt)
    grid = tc0 + np.arange(-3 * sig, 3 * sig + 1e-9, 0.5 / 1440)
    res = [fit_at(t, f, w, tc, exp_s) for tc in grid]
    chi = np.array([r[0] for r in res])
    dep = np.array([r[1] for r in res])
    der = np.array([r[2] for r in res])
    if not np.isfinite(chi).any():
        out["outcome"] = "insufficient data"
        return out
    i = int(np.nanargmin(chi))
    # no-transit chi2 (quadratic baseline only)
    dt = t - tc0
    A0 = np.vstack([np.ones_like(t), dt, dt ** 2]).T * w[:, None]
    c0, *_ = np.linalg.lstsq(A0, f * w, rcond=None)
    chi0 = float(np.sum((f * w - A0 @ c0) ** 2))
    dchi2 = chi0 - chi[i]
    out.update(best_tc=float(grid[i]), o_minus_c_min=float((grid[i] - tc0) * 1440), depth_ppm=float(dep[i]),
               depth_err_ppm=float(der[i]), delta_chi2=float(dchi2), at_grid_edge=bool(i in (0, len(grid) - 1)))
    detected = dchi2 >= DCHI2_DETECT and dep[i] > 0 and not out["at_grid_edge"]
    if detected and DEPTH_LO <= dep[i] <= DEPTH_HI:
        out["outcome"] = "recovered"
    elif detected:
        out["outcome"] = "inconsistent"
    elif coverage >= 0.8 and np.isfinite(dep).all() and np.all(dep + 3 * der < DEPTH_LO):
        out["outcome"] = "excluded"
    else:
        out["outcome"] = "insufficient data"
    return out


if __name__ == "__main__":
    print(json.dumps(main(sys.argv[1], int(sys.argv[2]), float(sys.argv[3]) if len(sys.argv) > 3 else 120.0), indent=1))
