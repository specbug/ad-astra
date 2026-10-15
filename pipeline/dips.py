"""Matched-filter search for individual transit-like dips in one light-curve segment set.

For each trial duration D the flattened flux is averaged in a sliding box of width D. The noise of that
box average is measured empirically as the robust spread of the box-averaged series itself, so red
noise on the timescale D is included automatically (white-noise error bars badly overstate the
significance of hours-long dips in TESS data). Local maxima of depth/noise above a threshold become
events; events found at several durations are merged, keeping the most significant duration.
"""
import numpy as np

from lcprep import robust_std

DURATIONS_H = np.array([0.75, 1.0, 1.5, 2.0, 3.0, 4.5, 6.5, 9.0, 13.0, 18.0])
SHORT_MAX_H = 4.5  # durations up to this use the short-window detrending


def box_mean(t, y, half):
    """Mean of y within [t_i - half, t_i + half] for every i, and the number of points used."""
    cs = np.concatenate([[0.0], np.cumsum(y)])
    lo = np.searchsorted(t, t - half, side="left")
    hi = np.searchsorted(t, t + half, side="right")
    n = hi - lo
    with np.errstate(invalid="ignore", divide="ignore"):
        m = (cs[hi] - cs[lo]) / n
    return m, n


def window_stat(t, y, center, half, func=np.nanmean):
    lo, hi = np.searchsorted(t, [center - half, center + half])
    seg = y[lo:hi]
    return func(seg) if seg.size else np.nan


def find_dips(t, y_short, y_long, cadence_d, snr_min=6.0, max_events=12, aux=None, seg_edges=None):
    """Return a list of event dicts for one light curve.

    t        times (sorted, days), already binned/cleaned
    y_short  flux flattened with a short window (for D <= SHORT_MAX_H)
    y_long   flux flattened with a long window (for longer D)
    aux      optional dict of auxiliary series aligned with t (bkg, cx, cy, flux_sml, flux_lag)
    """
    events = []
    sig_by_dur = {}
    if t.size < 50:
        return events, sig_by_dur
    for D_h in DURATIONS_H:
        D = D_h / 24.0
        y = y_short if D_h <= SHORT_MAX_H else y_long
        m, n = box_mean(t, y, D / 2)
        expected = D / cadence_d
        ok = n >= max(2, 0.6 * expected)
        if ok.sum() < 20:
            continue
        dep = 1.0 - m
        sd = robust_std(dep[ok])
        if not np.isfinite(sd) or sd <= 0:
            continue
        sig_by_dur[float(D_h)] = float(sd)
        snr = np.where(ok, dep / sd, -np.inf)
        cand = np.where(snr > snr_min)[0]
        if cand.size == 0:
            continue
        # local maxima within +-D
        order = cand[np.argsort(-snr[cand])]
        taken = np.zeros(t.size, bool)
        for i in order:
            if taken[i]:
                continue
            lo, hi = np.searchsorted(t, [t[i] - D, t[i] + D])
            if snr[i] < np.max(snr[lo:hi]):
                continue
            taken[lo:hi] = True
            events.append({"t0": float(t[i]), "dur_h": float(D_h), "depth": float(dep[i]),
                           "snr": float(snr[i]), "noise": float(sd), "npts": int(n[i]),
                           "coverage": float(n[i] / expected)})
            if len(events) > 400:
                break
    # merge across durations: events whose centres are within half the longer duration
    events.sort(key=lambda e: -e["snr"])
    merged = []
    for e in events:
        dup = False
        for k in merged:
            if abs(e["t0"] - k["t0"]) < 0.5 * max(e["dur_h"], k["dur_h"]) / 24.0:
                k["n_durations"] += 1
                dup = True
                break
        if not dup:
            e["n_durations"] = 1
            merged.append(e)
        if len(merged) >= max_events:
            break
    for e in merged:
        describe(e, t, y_short if e["dur_h"] <= SHORT_MAX_H else y_long, aux, seg_edges)
    return merged, sig_by_dur


def describe(e, t, y, aux, seg_edges):
    """Shape, context and pixel-level descriptors of one event (cheap, from light-curve columns)."""
    t0, D = e["t0"], e["dur_h"] / 24.0
    half = D / 2
    inb = (t >= t0 - half) & (t <= t0 + half)
    near = (np.abs(t - t0) < 3 * D + 0.25) & ~inb
    # U/V shape: central half of the box vs the outer quarters
    core = np.abs(t - t0) <= half / 2
    wing = inb & ~core
    e["depth_core"] = float(1 - np.nanmean(y[core])) if core.any() else np.nan
    e["depth_wing"] = float(1 - np.nanmean(y[wing])) if wing.any() else np.nan
    first = inb & (t < t0)
    second = inb & (t > t0)
    e["depth_first"] = float(1 - np.nanmean(y[first])) if first.any() else np.nan
    e["depth_second"] = float(1 - np.nanmean(y[second])) if second.any() else np.nan
    e["min_flux"] = float(np.nanmin(y[inb])) if inb.any() else np.nan
    e["local_rms"] = float(robust_std(y[near])) if near.sum() > 5 else np.nan
    e["n_near"] = int(near.sum())
    # distance to nearest data edge (gap or segment end), in days
    if seg_edges is not None and len(seg_edges):
        e["edge_dist_d"] = float(np.min(np.abs(np.asarray(seg_edges) - t0)))
    # local noise: spread of box means of the same width within +-3 d, excluding the event itself
    loc = (np.abs(t - t0) < 3.0)
    if loc.sum() > 20:
        tl, yl = t[loc], y[loc]
        bm, bn = box_mean(tl, yl, half)
        far = np.abs(tl - t0) > D
        okb = far & (bn >= 2)
        if okb.sum() > 10:
            ls = robust_std(1 - bm[okb])
            e["local_noise"] = float(ls)
            e["local_snr"] = float(e["depth"] / ls) if ls > 0 else np.nan
    if not aux:
        return
    for name, v in aux.items():
        if v is None:
            continue
        if name in ("sap_flux", "det_flux", "flux_sml", "flux_lag"):
            e[f"depth_{name}"] = _local_depth(t, v, t0, D)
            continue
        a_in = np.nanmean(v[inb]) if inb.any() else np.nan
        a_out = np.nanmedian(v[near]) if near.any() else np.nan
        s_out = robust_std(v[near]) if near.sum() > 5 else np.nan
        e[f"d_{name}"] = float(a_in - a_out)
        e[f"z_{name}"] = float((a_in - a_out) / s_out * np.sqrt(max(inb.sum(), 1))) \
            if np.isfinite(s_out) and s_out > 0 else np.nan


def _local_depth(t, v, t0, D):
    """Dip depth in a raw-ish flux series: in-box mean vs a line fitted to +-(D..3D+0.25d) around it."""
    inb = np.abs(t - t0) <= D / 2
    side = (np.abs(t - t0) > D) & (np.abs(t - t0) < 3 * D + 0.25)
    ok_in = inb & np.isfinite(v)
    ok_out = side & np.isfinite(v)
    if ok_in.sum() < 2 or ok_out.sum() < 6:
        return np.nan
    tt, vv = t[ok_out] - t0, v[ok_out]
    # need data on both sides to fit a slope; otherwise use a constant
    if (tt < 0).sum() >= 3 and (tt > 0).sum() >= 3:
        A = np.vstack([np.ones_like(tt), tt]).T
        coef, *_ = np.linalg.lstsq(A, vv, rcond=None)
        base = coef[0] + coef[1] * (t[ok_in] - t0)
    else:
        base = np.full(ok_in.sum(), np.median(vv))
    with np.errstate(invalid="ignore", divide="ignore"):
        return float(1 - np.mean(v[ok_in] / base))
