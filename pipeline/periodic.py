"""Multi-sector periodic transit search (Track A: small planets).

load_star       fetch every available sector (SPOC 2-min > TESS-SPOC > QLP) from S3, one product per sector
prepare_star    per-sector cleaning, noisy-sector rejection, biweight flattening, 30-min binning
bls_search      iterative coherent BLS over the full multi-year baseline (find, mask, repeat)
"""
import numpy as np
from astropy.timeseries import BoxLeastSquares

import s3io
from lcprep import bin_lc, detrend_fast, robust_std, ptp_std

PRODUCT_ORDER = ("spoc", "tspoc", "qlp")
DETREND_WINDOW = 0.75   # days; the reference project's injection tests favoured 0.75 d for 1-3 h transits
BIN_MIN = 30.0
PERIOD_MIN, PERIOD_MAX = 0.5, 30.0
MIN_TRANSITS = 3
OVERSAMPLE = 2   # period grid: transit drifts < D/2 (grazing) over the baseline
BLS_PHASE_OVERSAMPLE = 3


def load_star(tic, sectors_by_product):
    """sectors_by_product: {'spoc': [...], 'tspoc': [...], 'qlp': [...]} -> list of light-curve dicts."""
    taken, out = set(), []
    for prod in PRODUCT_ORDER:
        for s in sorted(sectors_by_product.get(prod, [])):
            if s in taken:
                continue
            lc = s3io.load(prod, tic, s)
            if lc is not None and lc["time"].size > 100:
                out.append(lc)
                taken.add(s)
    return out


def prepare_star(lcs, window=DETREND_WINDOW, bin_min=BIN_MIN, max_noise_ratio=3.0):
    """Flatten each sector, reject sectors far noisier than the star's typical one, bin and stitch."""
    parts, noise = [], {}
    for lc in lcs:
        t, f = lc["time"], lc["flux"]
        o = np.argsort(t)
        t, f = t[o], f[o]
        cad = lc["meta"]["cadence_s"] / 60
        if np.isfinite(cad) and cad < 9:
            t, f, _, _ = bin_lc(t, f, 10.0)
        flat, _ = detrend_fast(t, f, window, knot_step=window / 36)
        ok = np.isfinite(flat)
        t, flat = t[ok], flat[ok]
        if t.size < 50:
            continue
        sd = robust_std(flat)
        resid = flat - 1
        good = (resid < 3 * sd) & (resid > -10 * sd)   # clip flares; keep dips
        t, flat = t[good], flat[good]
        s = lc["meta"]["sector"]
        noise[s] = ptp_std(flat) * np.sqrt(10 / 60.0)  # per-hour point-to-point noise at 10-min bins
        parts.append((s, lc["meta"]["product"], t, flat))
    if not parts:
        return None
    typical = np.median(list(noise.values()))
    bad = {s for s, n in noise.items() if len(noise) >= 3 and n > max_noise_ratio * typical}
    T, Y, S = [], [], []
    for s, prod, t, y in parts:
        if s in bad:
            continue
        tb, yb, nb, _ = bin_lc(t, y, bin_min)
        T.append(tb); Y.append(yb); S.append(np.full(tb.size, s))
    t = np.concatenate(T); y = np.concatenate(Y); sec = np.concatenate(S)
    o = np.argsort(t)
    t, y, sec = t[o], y[o], sec[o]
    err = np.zeros_like(y)
    for s in np.unique(sec):
        m = sec == s
        err[m] = robust_std(y[m])
    return {"time": t, "flux": y, "flux_err": err, "sector": sec, "bad_sectors": sorted(bad),
            "noise_1h": {int(k): float(v) for k, v in noise.items()}}


def duration_days(period, rho_sun):
    """Central-transit duration of a circular orbit around a star of density rho (solar units)."""
    return 13.0 / 24 * (period / 365.25) ** (1 / 3) * rho_sun ** (-1 / 3)


def period_grid(baseline, rho_sun, pmin=PERIOD_MIN, pmax=PERIOD_MAX, oversample=OVERSAMPLE):
    """Frequency grid fine enough that a transit stays aligned (to D/oversample) across the baseline."""
    pmax = min(pmax, baseline / MIN_TRANSITS)
    f, fmin = 1 / pmin, 1 / pmax
    freqs = []
    while f > fmin:
        freqs.append(f)
        p = 1 / f
        f -= 0.5 * duration_days(p, rho_sun) / (oversample * baseline * p)
    return np.sort(1 / np.array(freqs))


def duration_grid(rho_sun, pmin=PERIOD_MIN, pmax=PERIOD_MAX, n=7):
    lo = 0.35 * duration_days(pmin, rho_sun)
    hi = min(1.3 * duration_days(pmax, rho_sun), 0.4 * pmin)
    return np.geomspace(max(lo, 0.025), max(hi, 0.05), n)


def sde(power, periods):
    """Signal detection efficiency after removing the slow trend of BLS power with period."""
    logp = np.log(periods)
    edges = np.linspace(logp.min(), logp.max(), 60)
    idx = np.clip(np.digitize(logp, edges) - 1, 0, len(edges) - 2)
    med = np.array([np.median(power[idx == i]) if np.any(idx == i) else np.nan for i in range(len(edges) - 1)])
    centers = 0.5 * (edges[1:] + edges[:-1])
    ok = np.isfinite(med)
    p = power - np.interp(logp, centers[ok], med[ok])
    return (p - np.median(p)) / np.std(p)


def transit_mask(t, period, t0, duration, factor=1.5):
    phase = ((t - t0 + 0.5 * period) % period) - 0.5 * period
    return np.abs(phase) < factor * duration / 2


def bls_search(lc, rho_sun, max_signals=3, invert=False):
    """Iterative BLS: find the best signal, record it, mask its transits, search again."""
    t, y, dy = lc["time"], lc["flux"].copy(), lc["flux_err"]
    if invert:
        y = 2.0 - y
    baseline = t.max() - t.min()
    if baseline < 3 * PERIOD_MIN:
        return []
    periods = period_grid(baseline, rho_sun)
    durations = duration_grid(rho_sun, pmax=periods.max())
    signals, mask = [], np.zeros(t.size, bool)
    for i in range(max_signals):
        if (~mask).sum() < 200:
            break
        bls = BoxLeastSquares(t[~mask], y[~mask], dy[~mask])
        res = bls.power(periods, durations, objective="snr", oversample=BLS_PHASE_OVERSAMPLE)
        s = sde(np.asarray(res.power), periods)
        k = int(np.nanargmax(s))
        P, t0, dur = float(res.period[k]), float(res.transit_time[k]), float(res.duration[k])
        st = bls.compute_stats(P, dur, t0)
        sig = {"signal": i + 1, "period": P, "t0": t0, "duration_hr": dur * 24,
               "depth_ppm": float(res.depth[k] * 1e6), "depth_err_ppm": float(res.depth_err[k] * 1e6),
               "snr": float(res.depth[k] / res.depth_err[k]), "sde": float(s[k]),
               "n_transits": int(np.sum(st["per_transit_count"] > 0)),
               "depth_odd_ppm": float(st["depth_odd"][0] * 1e6), "depth_odd_err_ppm": float(st["depth_odd"][1] * 1e6),
               "depth_even_ppm": float(st["depth_even"][0] * 1e6), "depth_even_err_ppm": float(st["depth_even"][1] * 1e6),
               "depth_half_ppm": float(st["depth_half"][0] * 1e6), "depth_half_err_ppm": float(st["depth_half"][1] * 1e6),
               "depth_phased_ppm": float(st["depth_phased"][0] * 1e6),
               "harmonic_dll": float(st["harmonic_delta_log_likelihood"]),
               "n_periods": int(periods.size), "baseline_d": float(baseline)}
        signals.append(sig)
        mask |= transit_mask(t, P, t0, dur)
    return signals
