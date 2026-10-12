"""Periods for single-transit TOIs and community TOIs (CTOIs) that have none, from every TESS sector.

Targets: ExoFOP TOIs and CTOIs (tables of 12 Oct 2026) whose period is missing or zero, minus false positives.
For each star:
  - every sector is loaded (SPOC 2-min > TESS-SPOC > QLP), flattened and binned to 30 min (duo.flatten_sectors);
  - transits of other planets already known on the star (TOIs/CTOIs with periods, confirmed transiting planets)
    are masked;
  - the catalogued transit is re-found, as a check that the data and detrending can see it;
  - every other part of the data is searched with a box of the catalogued duration for dips of similar depth
    (0.5-2x, at >= 7 sigma, with data on both sides of the dip);
  - two or more events: the periods dt/n not excluded by the data (duo.alias_table); one: the shortest period
    the data allow (duo.single_period_constraints).
Results: results/candidates/period_hunt.csv (one row per target) and period_hunt.json (details).

  python period_hunt.py [n_proc] [limit]
"""
import json
import math
import sys
import concurrent.futures as cf

import env  # noqa: F401
import numpy as np
import pandas as pd

from env import CATALOGS, RESULTS

E = CATALOGS / "exofop"
C = RESULTS / "candidates"
BJD0 = 2457000.0


def targets():
    toi = pd.read_csv(E / "toi_20261012.csv", low_memory=False)
    ctoi = pd.read_csv(E / "ctoi_20261012.csv", low_memory=False)
    tp = pd.to_numeric(toi["Period (days)"], errors="coerce")
    cp = pd.to_numeric(ctoi["Period (days)"], errors="coerce")
    rows = []
    for _, r in toi[(tp.isna() | (tp <= 0)) & ~toi["TFOPWG Disposition"].isin(["FP", "FA"])].iterrows():
        rows.append({"tic": int(r["TIC ID"]), "id": f"TOI-{r['TOI']}", "epoch": r["Epoch (BJD)"] - BJD0,
                     "depth_ppm": r["Depth (ppm)"], "dur_h": r["Duration (hours)"], "disp": r["TFOPWG Disposition"]})
    have = {(x["tic"], round(x["epoch"], 1)) for x in rows}
    for _, r in ctoi[(cp.isna() | (cp <= 0)) & (ctoi["User Disposition"] != "FP")].iterrows():
        ep = r["Transit Epoch (BJD)"] - BJD0
        if (int(r["TIC ID"]), round(ep, 1)) in have:
            continue
        rows.append({"tic": int(r["TIC ID"]), "id": f"CTOI {r['CTOI']}", "epoch": ep, "depth_ppm": r["Depth ppm"],
                     "dur_h": r["Duration (hrs)"], "disp": f"CTOI {r['User']}"})
    d = pd.DataFrame(rows)
    d = d[np.isfinite(d.epoch)]
    # other planets on the same star, to mask: TOIs/CTOIs with periods, confirmed transiting planets
    known = []
    for _, r in toi[tp > 0].iterrows():
        known.append((int(r["TIC ID"]), float(r["Period (days)"]), r["Epoch (BJD)"] - BJD0, r["Duration (hours)"]))
    for _, r in ctoi[cp > 0].iterrows():
        known.append((int(r["TIC ID"]), float(r["Period (days)"]), r["Transit Epoch (BJD)"] - BJD0, r["Duration (hrs)"]))
    return d.reset_index(drop=True), known


def box_search(flat, dur_h, exclude=()):
    """SNR of a box dip of the given duration centred on every 30-min bin; local maxima above 7 sigma."""
    t, y, e = flat
    D = dur_h / 24.0
    expected = max(D / (30 / 1440), 1)
    lo = np.searchsorted(t, t - D / 2)
    hi = np.searchsorted(t, t + D / 2)
    cs = np.concatenate([[0], np.cumsum(1 - y)])
    n = hi - lo
    depth = (cs[hi] - cs[lo]) / np.maximum(n, 1)
    snr = depth / (np.median(e) / np.sqrt(np.maximum(n, 1)))
    ok = n >= 0.6 * expected
    # data on both sides of the dip (rejects dips at segment edges)
    left = np.searchsorted(t, t - 1.5 * D)
    right = np.searchsorted(t, t + 1.5 * D)
    ok &= (lo - left >= 0.3 * expected) & (right - hi >= 0.3 * expected)
    for tc, w in exclude:
        ok &= np.abs(t - tc) > w
    snr = np.where(ok, snr, 0)
    peaks = []
    order = np.argsort(-snr)
    for i in order:
        if snr[i] < 7:
            break
        if all(abs(t[i] - p[0]) > D for p in peaks):
            peaks.append((float(t[i]), float(depth[i]), float(snr[i])))
        if len(peaks) >= 12:
            break
    return peaks


DURS = (1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 16.0, 20.0, 24.0)


def refine(flat, epoch, dur0):
    """Re-measure the catalogued transit: the duration and time that maximise the box SNR near the epoch
    (catalogued durations and depths are sometimes far off)."""
    t, y, e = flat
    best = None
    for d in DURS:
        D = d / 24
        win = max(D, dur0 / 24 if dur0 == dur0 else 0.25, 0.2)
        m = np.abs(t - epoch) < win
        if m.sum() < 3:
            continue
        for tc in t[m]:
            lo, hi = np.searchsorted(t, [tc - D / 2, tc + D / 2])
            n = hi - lo
            if n < max(0.6 * d * 2, 1):
                continue
            dep = 1 - np.mean(y[lo:hi])
            snr = dep / (np.median(e) / np.sqrt(n))
            if best is None or snr > best[3]:
                best = (float(tc), d, float(dep), float(snr))
    return best


def best_period(ref, events, depth, dur, flat):
    """Longest period, not excluded by the data, that puts the most events on one ephemeris."""
    from duo import alias_table
    D = dur / 24
    times = np.array([ref] + [o["t"] for o in events])
    cands = []
    for o in events:
        dt = abs(o["t"] - ref)
        nmax = int(dt / 1.0)
        P = dt / np.arange(1, max(nmax, 1) + 1)
        P = P[P >= max(1.0, 5 * D)]          # transit duty cycle below 20%
        k = np.round((times[None, :] - ref) / P[:, None])
        on = (np.abs(times[None, :] - ref - k * P[:, None]) < 0.25 * D).sum(1)
        cands += list(zip(on.tolist(), P.tolist()))
    if not cands:
        return None
    top = max(c[0] for c in cands)
    for count in range(top, 1, -1):
        for _, P in sorted({c for c in cands if c[0] == count}, key=lambda c: -c[1])[:30]:
            on_t = [tt for tt in times if abs(((tt - ref) / P) - round((tt - ref) / P)) * P < 0.25 * D]
            rows = alias_table([min(on_t), max(on_t)], 0.6 * depth, dur, flat, p_min=0.999 * P, n_max=None)
            ok = [a for a in rows if a["status"] == "allowed" and abs(a["P"] - P) < 1e-3 * P]
            if ok:
                return {"P": float(P), "n_events_on_it": int(count)}
    return None


def hunt(row, known):
    from dossier import all_sectors, best_lcs
    from duo import alias_table, flatten_sectors, single_period_constraints
    tic = int(row["tic"])
    out = {"tic": tic, "id": row["id"], "epoch": row["epoch"], "depth_ppm": row["depth_ppm"], "dur_h": row["dur_h"],
           "disp": row["disp"]}
    try:
        dur = float(row["dur_h"]) if row["dur_h"] == row["dur_h"] and row["dur_h"] > 0 else 6.0
        depth = float(row["depth_ppm"]) * 1e-6 if row["depth_ppm"] == row["depth_ppm"] and row["depth_ppm"] > 0 else None
        secs = all_sectors(tic)
        lcs = best_lcs(tic, secs)
        out["sectors"] = [lc["meta"]["sector"] for lc in lcs]
        flat = flatten_sectors(lcs, window=max(2.5, 4 * dur / 24))
        if flat is None:
            out["status"] = "no data"
            return out
        t, y, e = flat
        # mask other known planets on this star
        mask = np.ones(t.size, bool)
        for (k_tic, P, T0, d_h) in known:
            if k_tic != tic or not (P == P and T0 == T0):
                continue
            w = (d_h if d_h == d_h and d_h > 0 else 4) / 24
            ph = ((t - T0) / P + 0.5) % 1 - 0.5
            mask &= np.abs(ph * P) > 0.75 * w
        flat = (t[mask], y[mask], e[mask])
        t, y, e = flat
        # the catalogued event, re-measured
        out["in_data"] = bool(np.any(np.abs(t - row["epoch"]) < max(dur / 48, 0.1)))
        rf = refine(flat, row["epoch"], dur)
        out["recovered"] = bool(rf and rf[3] >= 7)
        out["recovered_snr"] = rf[3] if rf else None
        if out["recovered"]:
            ref, dur, d_ref = rf[0], rf[1], rf[2]
            out["ref_epoch"], out["dur_fit_h"], out["depth_fit_ppm"] = ref, dur, d_ref * 1e6
        else:
            ref, d_ref = row["epoch"], depth
        # all other events
        others = []
        for (tc, dd, s) in box_search(flat, dur, exclude=[(ref, dur / 24)]):
            if d_ref and not (0.5 * d_ref <= dd <= 2.0 * d_ref):
                continue
            others.append({"t": tc, "depth_ppm": dd * 1e6, "snr": s})
        out["others"] = others
        out["n_others"] = len(others)
        if d_ref is None:
            out["status"] = "no depth"
            return out
        if others:
            best = max(others, key=lambda o: o["snr"])
            rows = alias_table([ref, best["t"]], 0.6 * d_ref, dur, flat, p_min=1.0)
            ok = [a for a in rows if a["status"] == "allowed"]
            out["pair_dt"] = abs(best["t"] - ref)
            out["aliases_allowed"] = [round(a["P"], 4) for a in ok][:40]
            out["n_aliases_allowed"] = len(ok)
            out["best_alias"] = best_period(ref, others, d_ref, dur, flat)
            out["status"] = "new events"
        else:
            grid, allowed, pmin = single_period_constraints(ref, 0.6 * d_ref, dur, flat, pmin=2.0, pmax=3000.0)
            out["p_min_allowed"] = pmin
            out["status"] = "single"
    except Exception as ex:
        out["status"] = "error"
        out["error"] = repr(ex)[:200]
    return out


def _run(args):
    return hunt(*args)


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    limit = int(sys.argv[2]) if len(sys.argv) > 2 else None
    d, known = targets()
    if limit:
        d = d.head(limit)
    print(len(d), "targets", flush=True)
    res = []
    with cf.ProcessPoolExecutor(n) as ex:
        for i, r in enumerate(ex.map(_run, [(row, known) for row in d.to_dict("records")], chunksize=2)):
            res.append(r)
            if i % 50 == 0:
                print(i, r.get("tic"), r.get("status"), r.get("n_others"), flush=True)
    json.dump(res, open(C / "period_hunt.json", "w"), indent=1, default=float)
    flat = pd.DataFrame([{k: v for k, v in r.items() if k not in ("others",)} for r in res])
    flat.to_csv(C / "period_hunt.csv", index=False)
    print(flat.status.value_counts().to_string())


if __name__ == "__main__":
    main()
