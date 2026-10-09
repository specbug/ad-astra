"""Light-curve preparation: segmenting, binning, robust detrending."""
import numpy as np
from numba import njit
from wotan import flatten


def robust_std(x):
    x = x[np.isfinite(x)]
    if x.size == 0:
        return np.nan
    return 1.4826 * np.median(np.abs(x - np.median(x)))


def ptp_std(x):
    """Point-to-point scatter: insensitive to slow variability."""
    x = x[np.isfinite(x)]
    if x.size < 3:
        return np.nan
    return 1.4826 * np.median(np.abs(np.diff(x))) / np.sqrt(2)


def segments(t, gap=0.5):
    """Index ranges of contiguous data split at gaps longer than `gap` days."""
    if t.size == 0:
        return []
    br = np.where(np.diff(t) > gap)[0] + 1
    edges = np.concatenate([[0], br, [t.size]])
    return [(edges[i], edges[i + 1]) for i in range(len(edges) - 1)]


def bin_lc(t, f, bin_min=10.0, extra=None):
    """Average onto a regular grid of `bin_min` minutes (empty bins dropped).

    Returns t, f, n (points per bin) and binned copies of any arrays in `extra`.
    """
    w = bin_min / 1440.0
    if t.size == 0:
        return t, f, np.zeros(0, int), {k: v for k, v in (extra or {}).items()}
    k = np.floor((t - t[0]) / w).astype(np.int64)
    uk, inv, cnt = np.unique(k, return_inverse=True, return_counts=True)
    tb = np.bincount(inv, t) / cnt
    fb = np.bincount(inv, f) / cnt
    ex = {}
    for name, v in (extra or {}).items():
        ok = np.isfinite(v)
        num = np.bincount(inv[ok], v[ok], minlength=uk.size)
        den = np.bincount(inv[ok], minlength=uk.size)
        with np.errstate(invalid="ignore", divide="ignore"):
            ex[name] = np.where(den > 0, num / np.maximum(den, 1), np.nan)
    return tb, fb, cnt, ex


def detrend(t, f, window, method="biweight"):
    """Robust sliding-window trend removal per contiguous segment; returns flattened flux and trend."""
    flat, trend = flatten(t, f, method=method, window_length=window, edge_cutoff=0.0,
                          break_tolerance=0.5, return_trend=True)
    return flat, trend


@njit(cache=True)
def _biweight_knots(t, f, knots, half, c=5.0, n_iter=5):
    out = np.full(knots.size, np.nan)
    lo = 0
    hi = 0
    n = t.size
    buf = np.empty(n)
    for k in range(knots.size):
        tk = knots[k]
        while lo < n and t[lo] < tk - half:
            lo += 1
        if hi < lo:
            hi = lo
        while hi < n and t[hi] <= tk + half:
            hi += 1
        m = hi - lo
        if m < 3:
            continue
        for i in range(m):
            buf[i] = f[lo + i]
        loc = np.median(buf[:m])
        for it in range(n_iter):
            # MAD scale
            for i in range(m):
                buf[i] = abs(f[lo + i] - loc)
            mad = np.median(buf[:m])
            if mad <= 0:
                break
            num = 0.0
            den = 0.0
            for i in range(m):
                u = (f[lo + i] - loc) / (c * mad)
                if abs(u) < 1.0:
                    w = (1.0 - u * u) ** 2
                    num += w * (f[lo + i] - loc)
                    den += w
            if den <= 0:
                break
            step = num / den
            loc += step
            if abs(step) < 1e-7:
                break
        out[k] = loc
    return out


def detrend_fast(t, f, window, knot_step=1.0 / 48, gap=0.5):
    """Tukey biweight (c=5, as wotan) evaluated on a knot grid and linearly interpolated.

    About 10x faster than evaluating at every cadence; for windows >= 0.5 d and a 30-min knot step the
    trend differs from wotan's by well under the white noise. Segments are handled independently.
    """
    t = np.ascontiguousarray(t, dtype=np.float64)  # FITS columns are big-endian; numba needs native
    f = np.ascontiguousarray(f, dtype=np.float64)
    trend = np.full(t.size, np.nan)
    for a, b in segments(t, gap):
        ts, fs = t[a:b], f[a:b]
        if ts.size < 5:
            continue
        knots = np.arange(ts[0], ts[-1] + knot_step, knot_step)
        tk = _biweight_knots(ts, fs, knots, window / 2)
        ok = np.isfinite(tk)
        if ok.sum() < 2:
            continue
        trend[a:b] = np.interp(ts, knots[ok], tk[ok])
    return f / trend, trend


def clip_upward(f, sigma=4.0):
    """Remove flares/upward outliers only (downward points could be transits)."""
    sd = robust_std(f)
    med = np.nanmedian(f)
    return np.isfinite(f) & (f - med < sigma * sd)
