"""Period aliases of duo-transits: two (or more) single transits separated by a long gap.

The true period is P = dt / n for some integer n. Every alias predicts transits at other times; if any
prediction falls inside data we have and the data show no dip there, that alias is excluded.
"""
import numpy as np

from lcprep import bin_lc, detrend_fast, robust_std


def flatten_sectors(lcs, window=2.5):
    """Concatenate flattened, 30-min binned light curves of all sectors (for alias tests)."""
    T, Y, E = [], [], []
    for lc in lcs:
        t, f = lc["time"], lc["flux"]
        o = np.argsort(t)
        t, f = t[o], f[o]
        tb, fb, _, _ = bin_lc(t, f, 30.0)
        yb, _ = detrend_fast(tb, fb, window, knot_step=window / 36)
        ok = np.isfinite(yb)
        T.append(tb[ok]); Y.append(yb[ok]); E.append(np.full(ok.sum(), robust_std(yb[ok])))
    if not T:
        return None
    t = np.concatenate(T); y = np.concatenate(Y); e = np.concatenate(E)
    o = np.argsort(t)
    return t[o], y[o], e[o]


def alias_table(t_events, depth, dur_h, flat, p_min=None, n_max=None, cover_frac=0.5, exclude_sigma=3.0):
    """Return a list of dicts, one per alias n, with status 'allowed' or 'excluded'.

    t_events: times of the observed transits (days); depth: fractional depth; dur_h: duration (hours)
    flat: (t, y, err) from flatten_sectors
    """
    t_ev = np.sort(np.asarray(t_events, float))
    dt = t_ev[-1] - t_ev[0]
    D = dur_h / 24.0
    t, y, e = flat
    p_min = p_min or max(2 * D, 1.0)
    n_max = n_max or int(dt / p_min)
    cy = np.concatenate([[0.0], np.cumsum(y)])
    ce = np.concatenate([[0.0], np.cumsum(e)])
    rows = []
    for n in range(1, max(n_max, 1) + 1):
        P = dt / n
        if P < p_min:
            break
        # every observed event must sit on the alias grid (multi-event case)
        k = np.round((t_ev - t_ev[0]) / P)
        if np.any(np.abs(t_ev - t_ev[0] - k * P) > 0.5 * D):
            rows.append({"n": n, "P": P, "status": "inconsistent"})
            continue
        kmin = int(np.floor((t.min() - t_ev[0]) / P)) - 1
        kmax = int(np.ceil((t.max() - t_ev[0]) / P)) + 1
        # all predicted transits at once (vectorised; same test as one transit at a time)
        tc = t_ev[0] + np.arange(kmin, kmax + 1) * P
        tc = tc[np.min(np.abs(t_ev[None, :] - tc[:, None]), axis=1) >= 0.5 * D]   # not an observed event
        lo = np.searchsorted(t, tc - D / 2, side="right")
        hi = np.searchsorted(t, tc + D / 2, side="left")
        cnt = hi - lo
        expected = D / (30 / 1440)
        ok = cnt >= max(cover_frac * expected, 1)
        worst = None
        n_cov = int(ok.sum())
        if n_cov:
            m = 1 - (cy[hi[ok]] - cy[lo[ok]]) / cnt[ok]
            err = (ce[hi[ok]] - ce[lo[ok]]) / cnt[ok] / np.sqrt(cnt[ok])
            # significance with which the data reject a dip of the full depth here
            z = (depth - m) / err
            i = int(np.argmax(z))
            worst = (float(z[i]), float(tc[ok][i]), float(m[i]), float(err[i]))
        status = "excluded" if (worst is not None and worst[0] > exclude_sigma) else "allowed"
        rows.append({"n": n, "P": P, "status": status, "n_covered": n_cov,
                     "max_reject_sigma": None if worst is None else float(worst[0]),
                     "worst_time": None if worst is None else float(worst[1])})
    return rows


def single_period_constraints(t0, depth, dur_h, flat, pmin=5.0, pmax=3000.0, cover_frac=0.5, exclude_sigma=3.0):
    """Periods allowed for a single transit at t0, given every other observed data point.

    For each trial period P the predicted transits t0 + k*P (k != 0) that land in data with enough
    coverage are tested; P is excluded if any of them shows no dip (depth rejected at > exclude_sigma).
    Returns the period grid and a boolean 'allowed' mask, plus the minimum period not excluded.
    """
    t, y, e = flat
    D = dur_h / 24.0
    grid = np.exp(np.arange(np.log(pmin), np.log(pmax), D / 4 / pmax * 0 + 0.0005))
    allowed = np.ones(grid.size, bool)
    expected = D / (30 / 1440)
    cy = np.concatenate([[0.0], np.cumsum(y)])
    ce = np.concatenate([[0.0], np.cumsum(e)])
    for i, P in enumerate(grid):
        kmin = int(np.floor((t.min() - t0) / P)) - 1
        kmax = int(np.ceil((t.max() - t0) / P)) + 1
        k = np.arange(kmin, kmax + 1)
        tc = t0 + k[k != 0] * P
        lo = np.searchsorted(t, tc - D / 2)
        hi = np.searchsorted(t, tc + D / 2)
        n = hi - lo
        ok = (n >= cover_frac * expected) & (n > 0)
        if not ok.any():
            continue
        m = 1 - (cy[hi[ok]] - cy[lo[ok]]) / n[ok]
        err = (ce[hi[ok]] - ce[lo[ok]]) / n[ok] / np.sqrt(n[ok])
        if np.any((depth - m) / err > exclude_sigma):
            allowed[i] = False
    pmin_allowed = float(grid[allowed][0]) if allowed.any() else None
    return grid, allowed, pmin_allowed
