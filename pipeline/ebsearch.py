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
    # frequency grid: accumulated phase drift over the baseline below 1/3 of the shortest duration,
    # i.e. df <= dur/(3*base*P): a geometric grid in frequency
    q = dur.min() / (3 * base)
    nf = int(np.ceil(np.log(pmax / pmin) / q)) + 1
    periods = pmin * np.exp(np.linspace(0, np.log(pmax / pmin), min(nf, 200000)))
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


def analyse_eb(lc, window=3.0, return_arrays=False):
    """Full Track C analysis of one light curve dict (from s3io). Returns a result dict
    (and, with return_arrays, the binned time, flattened flux, EB template and residuals)."""
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
    # further periodic signals in the residuals (2+2 quadruples, sextuples): find, record, mask, repeat
    keep = np.ones(tr.size, bool)
    for j in (2, 3):
        if keep.sum() < 200:
            break
        try:
            tt, yy = tr[keep], rr[keep]
            Pj, t0j, durj, depj = eb_period(tt, yy)
            bls = BoxLeastSquares(tt, yy)
            st = bls.compute_stats(Pj, durj, t0j)
            intr = np.abs(((tt - t0j + 0.5 * Pj) % Pj) - 0.5 * Pj) < durj / 2
            snrj = depj / (robust_std(yy) / np.sqrt(max(int(intr.sum()), 1)))
            ratio = Pj / P
            near_harm = min(abs(ratio * m - round(ratio * m)) for m in (1, 2, 3, 4)) < 0.01
            out.update({f"P{j}": Pj, f"t0_{j}": t0j, f"dur{j}_h": durj * 24, f"depth{j}": depj,
                        f"snr{j}": float(snrj), f"n_tr{j}": int(np.sum(st["per_transit_count"] > 0)),
                        f"P{j}_harmonic": bool(near_harm),
                        f"depth{j}_odd": float(st["depth_odd"][0]), f"depth{j}_even": float(st["depth_even"][0])})
            # mask this signal's dips (and their phase-0.5 counterparts) before the next pass
            ph = np.abs(((tr - t0j + 0.5 * Pj) % Pj) - 0.5 * Pj)
            ph2 = np.abs(((tr - t0j) % Pj) - 0.5 * Pj)
            keep &= (ph > 0.75 * durj) & (ph2 > 0.75 * durj)
            if snrj < 8:
                break
        except Exception:
            break
    # isolated extra dips in residuals
    segs = segments(tr, 0.5)
    edges = [tr[a] for a, b in segs] + [tr[b - 1] for a, b in segs]
    ev, _ = find_dips(tr, rr, rr, 10 / 1440, snr_min=7.0, max_events=6, seg_edges=edges)
    out["extra_events"] = [{k: e[k] for k in ("t0", "dur_h", "depth", "snr", "edge_dist_d") if k in e} for e in ev]
    if return_arrays:
        return out, {"t": t, "y": y, "template": tmpl, "resid": resid, "eclipse_mask": ecl}
    return out
