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
        worst = None
        n_cov = 0
        for kk in range(kmin, kmax + 1):
            tc = t_ev[0] + kk * P
            if np.min(np.abs(t_ev - tc)) < 0.5 * D:
                continue  # an observed event
            inb = np.abs(t - tc) < D / 2
            expected = D / (30 / 1440)
            if inb.sum() < cover_frac * expected:
                continue
            n_cov += 1
            m = 1 - np.mean(y[inb])
            err = np.median(e[inb]) / np.sqrt(inb.sum())
            # significance with which the data reject a dip of the full depth here
            z = (depth - m) / err
            if worst is None or z > worst[0]:
                worst = (z, tc, m, err)
        status = "excluded" if (worst is not None and worst[0] > exclude_sigma) else "allowed"
        rows.append({"n": n, "P": P, "status": status, "n_covered": n_cov,
                     "max_reject_sigma": None if worst is None else float(worst[0]),
                     "worst_time": None if worst is None else float(worst[1])})
    return rows
