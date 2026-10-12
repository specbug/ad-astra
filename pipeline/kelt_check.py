"""Periods for TESS single transits from archival KELT photometry (2006-2018, NASA Exoplanet Archive).

For each star: the KELT light curves within 15" (both telescope orientations, TFA-detrended versions) are found
through the archive's kelttimeseries table and the file URLs in its bulk wget scripts.

The test: a real period P, anchored on the TESS transit (t0), puts transits in the KELT years at t0 - kP. For every
trial period in the range the TESS data still allow, the KELT points within +-0.4 T14 of those times are averaged,
night by night, with each night weighted by its white noise and a red-noise floor measured from the same light
curve (KELT nights have correlated errors, so one bad night must not look like many points). The grid of periods
is fine enough that a third of a transit duration of drift over the KELT-to-TESS gap is resolved, and the sum is
done per night with difference arrays, so the ~10^5-10^6 trial periods take seconds.

Calibration, per star:
  - the same search anchored on 40 false transit times (t0 + offsets of a few days, where TESS saw nothing):
    the distribution of their best SNR is the false-alarm level of this star's light curves;
  - the TESS transit injected into the KELT data at three periods: is it recovered, and at what SNR?
    (says whether KELT could have seen the planet at all).
Stars with two or more TESS transits are searched only at the periods those transits allow.

  python kelt_check.py [tic ...]  (default: the bright single transits in results/candidates/kelt_targets.json)
Output: results/candidates/kelt_check.json
"""
import glob
import io
import json
import math
import re
import sys
import urllib.parse
import urllib.request

import env  # noqa: F401
import numpy as np
import pandas as pd

from env import CATALOGS, RESULTS

C = RESULTS / "candidates"
CACHE = CATALOGS.parent / "cache" / "kelt"
WGET = None
BTJD0 = 2457000.0


def url_index():
    """kelt light-curve file name -> URL, from the archive's bulk wget scripts."""
    global WGET
    if WGET is None:
        WGET = {}
        for f in glob.glob(str(CATALOGS / "kelt_wget" / "*.bat")):
            for m in re.finditer(r"wget -O '([^']+)' '([^']+)'", open(f).read()):
                WGET[m.group(1)] = m.group(2).replace("http://exoplanetarchive.ipac.caltech.edu:80", "https://exoplanetarchive.ipac.caltech.edu")
    return WGET


def tap(q):
    u = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync?" + urllib.parse.urlencode({"query": q, "format": "csv"})
    return pd.read_csv(io.BytesIO(urllib.request.urlopen(u, timeout=180).read()))


def kelt_time(t):
    """The archive's KELT-South files carry JD + 2450000 (an offset added twice); KELT-North files carry JD."""
    t = np.asarray(t, float)
    if np.median(t) > 4.5e6:
        t = t - 2450000.0
    elif np.median(t) < 1e5:
        t = t + (2400000.0 if np.median(t) > 3e4 else 2450000.0)
    return t


def kelt_lcs(ra, dec, r_arcsec=15.0):
    r = r_arcsec / 3600
    q = (f"select kelt_sourceid, kelt_orientation, ra, dec, kelt_mag, npts from kelttimeseries where "
         f"ra between {ra - r / math.cos(math.radians(dec))} and {ra + r / math.cos(math.radians(dec))} "
         f"and dec between {dec - r} and {dec + r}")
    d = tap(q)
    out = []
    idx = url_index()
    CACHE.mkdir(parents=True, exist_ok=True)
    seen = set()
    for (sid, orient, kmag) in zip(d.kelt_sourceid, d.kelt_orientation, d.kelt_mag):
        base = sid if sid.endswith(("_east", "_west")) else f"{sid}_{orient}"
        name = f"{base}_tfa_lc.tbl"
        url = idx.get(name)
        if not url or name in seen:   # the table can list one light curve twice
            continue
        seen.add(name)
        p = CACHE / name
        if not p.exists():
            p.write_bytes(urllib.request.urlopen(url, timeout=180).read())
        t, m, e = read_tbl(p)
        if t.size > 50:
            out.append({"name": name, "t": kelt_time(t), "mag": m, "err": e, "kelt_mag": float(kmag)})
    return out


def read_tbl(p):
    rows = []
    for line in open(p, errors="ignore"):
        if line.startswith(("\\", "|")) or not line.strip():
            continue
        parts = line.split()
        try:
            rows.append([float(x) for x in parts[:3]])
        except ValueError:
            continue
    a = np.array(rows) if rows else np.zeros((0, 3))
    return a[:, 0], a[:, 1], a[:, 2]


def prep(kelt, dur_h):
    """KELT light curves -> nights [(t, 1 - flux, white sigma, red sigma)] and a summary per light curve."""
    D = dur_h / 24
    nights, info = [], []
    for lc in kelt:
        o = np.argsort(lc["t"])
        t = lc["t"][o]
        f = 10 ** (-0.4 * (lc["mag"][o] - np.median(lc["mag"])))
        segs = np.split(np.arange(t.size), np.where(np.diff(t) > 0.25)[0] + 1)
        # night-by-night zero points, unless the transit is long enough to fill much of a night
        if dur_h <= 4:
            for s in segs:
                if s.size >= 5:
                    f[s] /= np.median(f[s])
        f /= np.median(f)
        mad = 1.4826 * np.median(np.abs(f - 1))
        ok = np.abs(f - 1) < 5 * mad
        sw = 1.4826 * np.median(np.abs(np.diff(f[ok]))) / math.sqrt(2)   # point-to-point (white) scatter
        y = 1 - f
        # red noise: scatter of means over windows as long as the in-transit window (0.8 T14)
        h = 0.4 * D
        lo, hi = np.searchsorted(t, t - h), np.searchsorted(t, t + h)
        yc = np.concatenate([[0], np.cumsum(np.where(ok, y, 0))])
        nc = np.concatenate([[0], np.cumsum(ok)])
        n = nc[hi] - nc[lo]
        m = (yc[hi] - yc[lo]) / np.maximum(n, 1)
        good = n >= 3
        sr2 = (1.4826 * np.median(np.abs(m[good]))) ** 2 - sw ** 2 / np.median(n[good]) if good.any() else 0.0
        sr = math.sqrt(max(sr2, 0.0))
        for s in segs:
            s = s[ok[s]]
            if s.size >= 3:
                nights.append((t[s], y[s], sw, sr))
        info.append({"name": lc["name"], "n": int(ok.sum()), "nights": len(segs), "kelt_mag": lc["kelt_mag"],
                     "span": [float(t.min()), float(t.max())], "white_ppm": sw * 1e6, "red_ppm": sr * 1e6})
    return nights, info


def search(nights, t0, D, p_min, p_max, r, inject=None):
    """Weighted mean depth, SNR and number of nights in transit for the periods p_min * r^i, anchored at t0
    (all KELT data precede t0). inject = (P, depth): a box transit added to the data first."""
    h = 0.4 * D
    lnr = math.log(r)
    nP = int(math.log(p_max / p_min) / lnr) + 1
    A = np.zeros(nP); B = np.zeros(nP); Nn = np.zeros(nP)
    step = D / 20
    for (tn, yn, sw, sr) in nights:
        if inject:
            Pi, di = inject
            ph = ((tn - t0) / Pi + 0.5) % 1 - 0.5
            yn = yn + di * (np.abs(ph * Pi) < D / 2)
        tau = np.arange(tn[0] - h, tn[-1] + h + step, step)
        lo, hi = np.searchsorted(tn, tau - h), np.searchsorted(tn, tau + h)
        cs = np.concatenate([[0], np.cumsum(yn)])
        n = hi - lo
        m = (cs[hi] - cs[lo]) / np.maximum(n, 1)
        w = 1 / (sw ** 2 / np.maximum(n, 1) + sr ** 2)
        da, db = t0 - tau[-1], t0 - tau[0]
        if da <= 0:
            continue
        kmin, kmax = max(int(math.ceil(da / p_max)), 1), int(math.floor(db / p_min))
        if kmax < kmin:
            continue
        ks = np.arange(kmin, kmax + 1)
        ia = np.maximum(np.ceil(np.log(da / ks / p_min) / lnr).astype(int), 0)
        ib = np.minimum(np.floor(np.log(db / ks / p_min) / lnr).astype(int), nP - 1)
        cnt = np.maximum(ib - ia + 1, 0)
        tot = int(cnt.sum())
        if tot == 0:
            continue
        kk = np.repeat(ks, cnt)
        ii = np.repeat(ia, cnt) + (np.arange(tot) - np.repeat(np.cumsum(cnt) - cnt, cnt))
        j = np.clip(np.rint((t0 - kk * p_min * np.exp(ii * lnr) - tau[0]) / step).astype(int), 0, tau.size - 1)
        g = n[j] >= 1   # KELT samples a field every ~20-40 min: one point per night can be all there is
        ii, j = ii[g], j[g]
        A += np.bincount(ii, weights=m[j] * w[j], minlength=nP)
        B += np.bincount(ii, weights=w[j], minlength=nP)
        Nn += np.bincount(ii, minlength=nP)
    P = p_min * np.exp(np.arange(nP) * lnr)
    with np.errstate(invalid="ignore", divide="ignore"):
        depth = np.where(B > 0, A / B, 0.0)
        snr = np.where(B > 0, A / np.sqrt(B), 0.0)
    return P, depth, snr, Nn


def allowed_by_tess(P, events_btjd, D):
    """Periods that put every TESS transit on one ephemeris (within a quarter of the duration)."""
    if not events_btjd or len(events_btjd) < 2:
        return np.ones(P.size, bool)
    ok = np.ones(P.size, bool)
    e0 = events_btjd[0]
    for e in events_btjd[1:]:
        dt = e - e0
        ok &= np.abs(dt / P - np.round(dt / P)) * P < 0.25 * D
    return ok


def peaks(P, depth, snr, Nn, mask, r, min_nights=3, n=8):
    s = np.where(mask & (Nn >= min_nights) & (depth > 0), snr, 0)
    out = []
    for i in np.argsort(-s):
        if s[i] <= 0 or len(out) >= n:
            break
        if all(abs(math.log(P[i] / q["P"])) > 3 * (r - 1) for q in out):
            out.append({"P": float(P[i]), "snr": float(s[i]), "depth": float(depth[i]), "nights": int(Nn[i])})
    return out


def nights_in_transit(nights, t0, P, D):
    """Per-night mean depth for one trial period, to look at by eye."""
    out = []
    for (tn, yn, sw, sr) in nights:
        ph = ((tn - t0) / P + 0.5) % 1 - 0.5
        inn = np.abs(ph * P) < 0.4 * D
        if inn.sum() >= 1:
            out.append({"t": float(np.mean(tn[inn])), "n": int(inn.sum()), "depth": float(np.mean(yn[inn])),
                        "err": float(math.sqrt(sw ** 2 / inn.sum() + sr ** 2))})
    return out


def check(tg, n_null=40, seed=1):
    kelt = kelt_lcs(tg["ra"], tg["dec"])
    if not kelt:
        return {"status": "no KELT data"}
    t0, D = tg["t0_bjd"], tg["dur_h"] / 24
    nights, info = prep(kelt, tg["dur_h"])
    if not nights:
        return {"status": "no usable KELT nights", "lcs": info}
    tk = np.concatenate([x[0] for x in nights])
    gap = t0 - np.median(tk)
    r = 1 + D / gap / 3
    p_min, p_max = max(tg.get("p_min") or 10.0, 2 * D), 1500.0
    events = tg.get("events")
    P, dep, snr, Nn = search(nights, t0, D, p_min, p_max, r)
    mask = allowed_by_tess(P, events, D)
    top = peaks(P, dep, snr, Nn, mask, r)
    for q in top[:3]:
        q["per_night"] = nights_in_transit(nights, t0, q["P"], D)
    # false-alarm level: the same search anchored where TESS saw nothing
    rng = np.random.default_rng(seed)
    null = []
    for off in rng.uniform(0.6, 6.0, n_null) * rng.choice([-1, 1], n_null):
        Pn, dn, sn, Nnn = search(nights, t0 + off, D, p_min, p_max, r)
        pk = peaks(Pn, dn, sn, Nnn, allowed_by_tess(Pn, events, D), r, n=1)
        null.append(pk[0]["snr"] if pk else 0.0)
    null = np.array(null)
    # sensitivity: the TESS transit injected at three periods
    inj = []
    for Pi in sorted({round(x, 4) for x in (1.3 * p_min, 2.5 * p_min, 6 * p_min) if x < p_max}):
        if events and len(events) >= 2 and not allowed_by_tess(np.array([Pi]), events, D)[0]:
            # pick the nearest period the TESS transits allow
            cand = P[mask]
            if cand.size == 0:
                continue
            Pi = float(cand[np.argmin(np.abs(cand - Pi))])
        Pj, dj, sj, Nj = search(nights, t0, D, p_min, p_max, r, inject=(Pi, tg["depth"]))
        i = int(np.argmin(np.abs(np.log(Pj / Pi))))
        w = slice(max(i - 2, 0), i + 3)
        best = peaks(Pj, dj, sj, Nj, mask, r, n=1)
        inj.append({"P": float(Pi), "snr_at_P": float(np.max(sj[w])), "nights": int(np.max(Nj[w])),
                    "top_is_injected": bool(best and abs(math.log(best[0]["P"] / Pi)) < 3 * (r - 1))})
    return {"status": "searched", "lcs": info, "n_points": int(tk.size), "n_nights": len(nights),
            "gap_days": float(gap), "grid_ratio": float(r), "n_trials": int(P.size), "n_trials_allowed": int(mask.sum()),
            "p_range": [float(p_min), p_max], "depth_expected": float(tg["depth"]), "top": top,
            "null_best_snr": {"median": float(np.median(null)), "p95": float(np.percentile(null, 95)),
                              "max": float(null.max())},
            "fap_top": float((null >= top[0]["snr"]).mean()) if top else None, "injections": inj}


def main(tics=None):
    targets = json.load(open(C / "kelt_targets.json"))
    out_p = C / "kelt_check.json"
    out = json.load(open(out_p)) if (tics and out_p.exists()) else {}
    for tg in targets:
        if tics and tg["tic"] not in tics:
            continue
        try:
            r = check(tg)
        except Exception as ex:
            r = {"status": "error", "error": repr(ex)[:200]}
        out[str(tg["tic"])] = r
        top = r.get("top", [{}])[0] if r.get("top") else {}
        print(tg["tic"], r["status"], r.get("n_nights"), "nights;",
              {k: round(v, 4) if isinstance(v, float) else v for k, v in top.items() if k != "per_night"},
              "null", r.get("null_best_snr"), "fap", r.get("fap_top"),
              "inj", [(round(x["P"], 2), round(x["snr_at_P"], 1), x["top_is_injected"]) for x in r.get("injections", [])],
              flush=True)
        json.dump(out, open(out_p, "w"), indent=1)


if __name__ == "__main__":
    main([int(a) for a in sys.argv[1:]] or None)
