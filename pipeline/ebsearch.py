"""Eclipsing-binary analysis for Track C (eclipsing triples and 2+2 quadruples).

For a light curve with a strong periodic eclipse signal:
  1. find the EB period with BLS (and decide between P and 2P from primary/secondary depths),
  2. build a phase-binned median template of the EB and subtract it,
  3. mask the primary/secondary eclipse phases and search the residuals for
       - a second periodic eclipse signal (BLS)        -> 2+2 quadruple candidate,
       - isolated extra eclipses (matched filter)       -> eclipsing triple / tertiary eclipse candidate.
"""
import numpy as np
from astropy.timeseries import BoxLeastSquares

from lcprep import bin_lc, detrend_fast, robust_std, segments
from dips import find_dips


def eb_period(t, y, pmin=0.2, pmax=None, durations_h=(1, 2, 3, 5, 8)):
    base = t.max() - t.min()
    pmax = pmax or min(15.0, base / 2.0)
    bls = BoxLeastSquares(t, y)
    dur = np.array(durations_h) / 24
    dur = dur[dur < 0.3 * pmin * 3]  # keep durations shorter than the shortest period scale
    if dur.size == 0:
        dur = np.array([0.03])
    # frequency grid: phase drift over the baseline < 1/4 of the shortest duration
    df = 0.25 * dur.min() / base ** 2
    freqs = np.arange(1 / pmax, 1 / pmin, df)
    if freqs.size > 300000:
        freqs = np.linspace(1 / pmax, 1 / pmin, 300000)
    periods = np.sort(1 / freqs)
    res = bls.power(periods, dur, objective="likelihood")
    k = int(np.argmax(res.power))
    return float(res.period[k]), float(res.transit_time[k]), float(res.duration[k]), float(res.depth[k])


def fold_template(phase, y, nbins=400):
    edges = np.linspace(0, 1, nbins + 1)
    idx = np.clip(np.digitize(phase, edges) - 1, 0, nbins - 1)
    med = np.full(nbins, np.nan)
    for i in range(nbins):
        v = y[idx == i]
        if v.size:
            med[i] = np.median(v)
    ok = np.isfinite(med)
    c = 0.5 * (edges[1:] + edges[:-1])
    if ok.sum() < 3:
        return np.ones_like(y)
    # periodic interpolation
    cc = np.concatenate([c[ok] - 1, c[ok], c[ok] + 1])
    mm = np.concatenate([med[ok]] * 3)
    return np.interp(phase, cc, mm)


def analyse_eb(lc, window=3.0):
    """Full Track C analysis of one light curve dict (from s3io). Returns a result dict."""
    t, f = lc["time"], lc["flux"]
    o = np.argsort(t)
    t, f = t[o], f[o]
    t, f, _, _ = bin_lc(t, f, 10.0)
    y, _ = detrend_fast(t, f, window, knot_step=window / 36)
    ok = np.isfinite(y)
    t, y = t[ok], y[ok]
    out = {"n": int(t.size)}
    if t.size < 200:
        return out
    P, t0, dur, dep = eb_period(t, y)
    # P vs 2P: compare depths of odd and even eclipses at P
    ph = ((t - t0) / P) % 1.0
    ine = (np.abs(((t - t0 + 0.5 * P) % P) - 0.5 * P) < dur / 2)
    cyc = np.floor((t - t0 + 0.5 * P) / P).astype(int)
    d_odd = 1 - np.median(y[ine & (cyc % 2 == 1)]) if (ine & (cyc % 2 == 1)).any() else np.nan
    d_even = 1 - np.median(y[ine & (cyc % 2 == 0)]) if (ine & (cyc % 2 == 0)).any() else np.nan
    sd = robust_std(y)
    if np.isfinite(d_odd) and np.isfinite(d_even) and abs(d_odd - d_even) > max(5 * sd / 3, 0.1 * max(d_odd, d_even)):
        P *= 2
        if d_odd > d_even:
            t0 += P / 2
    ph = ((t - t0) / P) % 1.0
    tmpl = fold_template(ph, y, nbins=int(np.clip(P / (10 / 1440) / 2, 100, 1000)))
    resid = y - tmpl + 1.0
    # eclipse phases: wherever the template dips by > 3 sigma of the residual scatter
    rs = robust_std(resid)
    ecl = (1 - tmpl) > 3 * rs
    # grow the eclipse mask by one duration
    if ecl.any():
        grow = np.zeros_like(ecl)
        tt = t[ecl]
        idx = np.searchsorted(tt, t)
        left = np.abs(t - tt[np.clip(idx - 1, 0, tt.size - 1)])
        right = np.abs(t - tt[np.clip(idx, 0, tt.size - 1)])
        grow = np.minimum(left, right) < 0.5 * dur
        ecl = ecl | grow
    out.update({"P1": P, "t0_1": t0, "dur1_h": dur * 24, "depth1": float(1 - tmpl.min()),
                "rms": float(sd), "resid_rms": float(rs), "ecl_frac": float(ecl.mean())})
    tr, rr = t[~ecl], resid[~ecl]
    if tr.size < 200:
        return out
    # second periodic signal in residuals
    try:
        P2, t02, dur2, dep2 = eb_period(tr, rr)
        bls = BoxLeastSquares(tr, rr)
        st = bls.compute_stats(P2, dur2, t02)
        intr = np.abs(((tr - t02 + 0.5 * P2) % P2) - 0.5 * P2) < dur2 / 2
        n_in = int(intr.sum())
        snr2 = dep2 / (robust_std(rr) / np.sqrt(max(n_in, 1)))
        ratio = P2 / P
        near_harm = min(abs(ratio * m - round(ratio * m)) for m in (1, 2, 3, 4)) < 0.01
        out.update({"P2": P2, "t0_2": t02, "dur2_h": dur2 * 24, "depth2": dep2, "snr2": float(snr2),
                    "n_tr2": int(np.sum(st["per_transit_count"] > 0)), "P2_harmonic": bool(near_harm),
                    "depth2_odd": float(st["depth_odd"][0]), "depth2_even": float(st["depth_even"][0])})
    except Exception:
        pass
    # isolated extra dips in residuals
    segs = segments(tr, 0.5)
    edges = [tr[a] for a, b in segs] + [tr[b - 1] for a, b in segs]
    ev, _ = find_dips(tr, rr, rr, 10 / 1440, snr_min=7.0, max_events=6, seg_edges=edges)
    out["extra_events"] = [{k: e[k] for k in ("t0", "dur_h", "depth", "snr", "edge_dist_d") if k in e} for e in ev]
    return out
