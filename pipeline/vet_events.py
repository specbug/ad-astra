"""Stage-2 vetting of individual dips (Tracks B and D), run only on survey candidates.

For each event the light curve is re-fetched and three things are measured:
  1. dip-vs-bump: the strongest *upward* excursion of the same duration elsewhere in the light curve.
     Noise and stellar variability make bumps as well as dips; a transit makes only a dip.
  2. a local fit of a quadratic baseline times a trapezoid (symmetric) transit, compared with a
     baseline-only model, an asymmetric trapezoid (separate ingress/egress) and a Gaussian dip.
     Gives depth, total duration T14, ingress fraction, asymmetry and fit-quality numbers.
  3. depth in the raw (SAP) flux with the same baseline treatment.
"""
import numpy as np
from scipy.optimize import least_squares

from lcprep import bin_lc, detrend_fast, robust_std
from dips import box_mean


def dip_vs_bump(t, y, t0, dur_h):
    D = dur_h / 24
    m, n = box_mean(t, y, D / 2)
    ok = n >= max(2, 0.6 * D / np.median(np.diff(t)))
    dev = (m - 1.0)
    sd = robust_std(dev[ok])
    far = np.abs(t - t0) > 1.5 * D
    up = np.nanmax(dev[ok & far]) / sd if (ok & far).any() else np.nan
    down_other = np.nanmax(-dev[ok & far]) / sd if (ok & far).any() else np.nan
    here = np.nanmax(-dev[ok & ~far]) / sd if (ok & ~far).any() else np.nan
    return {"snr_here": float(here), "max_bump_snr": float(up), "max_other_dip_snr": float(down_other)}


def _trap(t, t0, depth, t14, frac_in, frac_eg=None):
    """Trapezoid dip: ingress lasts frac_in*t14/2, egress frac_eg*t14/2 (symmetric if frac_eg None)."""
    if frac_eg is None:
        frac_eg = frac_in
    a = t0 - t14 / 2
    b = t0 + t14 / 2
    tin = max(frac_in * t14 / 2, 1e-5)
    teg = max(frac_eg * t14 / 2, 1e-5)
    x = np.zeros_like(t)
    m = (t > a) & (t < b)
    x[m] = 1.0
    r = (t > a) & (t < a + tin)
    x[r] = (t[r] - a) / tin
    f = (t < b) & (t > b - teg)
    x[f] = np.minimum(x[f], (b - t[f]) / teg)
    return depth * x


def fit_event(t, f, t0, dur_h, win_mult=3.0):
    """Local fits around one event. f is normalised (not flattened) flux; returns a dict of metrics."""
    D = dur_h / 24
    w = max(win_mult * D, 0.5)
    m = np.abs(t - t0) < w + D / 2
    tt, ff = t[m], f[m]
    out = {"n_fit": int(m.sum())}
    if tt.size < 20:
        return out
    x = tt - t0
    oot = np.abs(x) > D
    sig = robust_std(np.diff(ff[oot])) / np.sqrt(2) if oot.sum() > 5 else robust_std(ff)
    sig = sig if sig > 0 else 1e-3

    def resid_base(p):
        return (ff - (p[0] + p[1] * x + p[2] * x ** 2)) / sig

    def resid_trap(p):
        base = p[0] + p[1] * x + p[2] * x ** 2
        return (ff - base * (1 - _trap(x, p[3], p[4], p[5], p[6]))) / sig

    def resid_atrap(p):
        base = p[0] + p[1] * x + p[2] * x ** 2
        return (ff - base * (1 - _trap(x, p[3], p[4], p[5], p[6], p[7]))) / sig

    def resid_gauss(p):
        base = p[0] + p[1] * x + p[2] * x ** 2
        return (ff - base * (1 - p[4] * np.exp(-0.5 * ((x - p[3]) / p[5]) ** 2))) / sig

    med = np.median(ff[oot]) if oot.any() else 1.0
    depth0 = max(med - np.median(ff[np.abs(x) < D / 2]), 1e-5) if (np.abs(x) < D / 2).any() else 1e-3
    try:
        rb = least_squares(resid_base, [med, 0, 0])
        rt = least_squares(resid_trap, [med, 0, 0, 0, depth0, D, 0.3],
                           bounds=([0, -1, -1, -D / 2, 0, 0.2 * D, 0.01], [2, 1, 1, D / 2, 1, 3 * D, 1.0]))
        ra = least_squares(resid_atrap, list(rt.x) + [rt.x[6]],
                           bounds=([0, -1, -1, -D / 2, 0, 0.2 * D, 0.01, 0.01], [2, 1, 1, D / 2, 1, 3 * D, 1.0, 1.0]))
        rg = least_squares(resid_gauss, [med, 0, 0, 0, depth0, D / 3],
                           bounds=([0, -1, -1, -D / 2, 0, 0.02 * D], [2, 1, 1, D / 2, 1, 2 * D]))
    except Exception as e:
        out["fit_error"] = repr(e)[:80]
        return out
    chi_b, chi_t, chi_a, chi_g = (2 * r.cost for r in (rb, rt, ra, rg))
    p = rt.x
    out.update({"fit_depth": float(p[4]), "fit_t14_h": float(p[5] * 24), "fit_ingress_frac": float(p[6]),
                "fit_t0": float(t0 + p[3]), "dchi2_trap": float(chi_b - chi_t), "dchi2_asym": float(chi_t - chi_a),
                "dchi2_trap_vs_gauss": float(chi_g - chi_t), "redchi2_trap": float(chi_t / max(tt.size - 7, 1)),
                "fit_in_frac": float(ra.x[6]), "fit_eg_frac": float(ra.x[7]), "fit_sigma": float(sig)})
    return out


def vet_one(lc, t0, dur_h):
    """All stage-2 metrics for one event in one light-curve dict (from s3io.load)."""
    t, f = lc["time"], lc["flux"]
    o = np.argsort(t)
    t, f = t[o], f[o]
    cad = lc["meta"]["cadence_s"] / 60
    tb, fb, _, _ = bin_lc(t, f, 10.0) if (np.isfinite(cad) and cad < 9) else (t, f, None, None)
    win = 0.75 if dur_h <= 4.5 else 2.5
    yb, _ = detrend_fast(tb, fb, win, knot_step=win / 36)
    ok = np.isfinite(yb)
    res = dip_vs_bump(tb[ok], yb[ok], t0, dur_h)
    res.update(fit_event(t, f, t0, dur_h))
    if "sap_flux" in lc:
        s = lc["sap_flux"][o]
        okk = np.isfinite(s)
        fs = fit_event(t[okk], s[okk] / np.nanmedian(s[okk]), t0, dur_h)
        res["sap_fit_depth"] = fs.get("fit_depth", np.nan)
        res["sap_dchi2_trap"] = fs.get("dchi2_trap", np.nan)
    return res
