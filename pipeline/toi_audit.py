"""Ephemeris audit of the TOI catalogue: do the transits fall where ExoFOP's TOI table says they do?

For every TOI with a period (ExoFOP TOI table of 12 Oct 2026; dispositions other than FP and FA), each transit
that its listed ephemeris predicts inside a TESS light curve (one per sector: SPOC 2-min > TESS-SPOC > QLP) is
measured like the period review's events (period_review.measure: box depth against a local quadratic baseline,
red-noise error). A predicted transit is "seen" when its SNR >= 3 and its depth is more than 0.4x the catalogued
depth, and "missing" when it is covered (>= 50%) and shallower than a third of the catalogued depth by more than
3 sigma. Only TOIs whose single transits are detectable (catalogued depth >= 3x the typical error) are judged.

A TOI is flagged when at least three predicted transits are missing and they are at least a third of those
covered. For a flagged TOI, alternatives are searched through anchors (the seen transits of the listed
ephemeris; if there are none, the strongest dips of the right depth in the sector of the listed epoch): periods
within max(0.5%, 5 sigma) of the listed one and of 1/3, 1/2, 2 and 3 times it are scored on 30-minute binned,
flattened light curves, and the best are re-measured transit by transit. The best alternative that leaves at most
one in ten covered transits missing is kept if it sees clearly more transits than the listed ephemeris.

  N_PROC=4 python toi_audit.py [max_stars]
Output: results/candidates/toi_audit.json (every judged TOI, with alternatives for the flagged ones);
checkpoint WORK/toi_audit.jsonl
"""
import json
import math
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed

import env  # noqa: F401
import numpy as np
import pandas as pd

from env import CATALOGS, RESULTS

E = CATALOGS / "exofop"
C = RESULTS / "candidates"
OUT = C / "toi_audit.json"
CKPT = CATALOGS.parent / "toi_audit.jsonl"
BJD0 = 2457000.0


def tois():
    t = pd.read_csv(E / "toi_20261012.csv", low_memory=False)
    P = pd.to_numeric(t["Period (days)"], errors="coerce")
    # planet candidates only: confirmed and known planets have published ephemerides of their own
    t = t[(P > 0) & t["TFOPWG Disposition"].isin(["PC", "APC"]) & t["Epoch (BJD)"].notna()
          & t["Duration (hours)"].notna() & t["Depth (ppm)"].notna()]
    return t


def sectors_for(tics):
    """dossier.all_sectors for many stars at once (one filtered read of each product manifest)."""
    out = {int(t): {"spoc": [], "tspoc": [], "qlp": []} for t in tics}
    for p in ("spoc", "tspoc", "qlp"):
        m = pd.read_parquet(CATALOGS / f"{p}_manifest.parquet", columns=["tic", "sector"],
                            filters=[("tic", "in", [int(t) for t in tics])])
        for tic, g in m.groupby("tic"):
            out[int(tic)][p] = sorted(g.sector.astype(int).tolist())
    return out


def predicted(lcs, T0, P, D):
    out = []
    for lc in lcs:
        t = lc["time"]
        lo, hi = float(np.nanmin(t)), float(np.nanmax(t))
        out += [(n, T0 + n * P) for n in range(math.ceil((lo + D / 2 - T0) / P), math.floor((hi - D / 2 - T0) / P) + 1)]
    return sorted(set(out))


def series(lcs, T0, P, D, dref, cap=3000):
    from period_review import lc_for, measure
    ev = []
    pts = predicted(lcs, T0, P, D)
    if len(pts) > cap:   # very short periods over many sectors: a even subsample is enough to judge
        idx = np.linspace(0, len(pts) - 1, cap).astype(int)
        pts = [pts[i] for i in idx]
    for n, tc in pts:
        lc = lc_for(lcs, tc)
        m = measure(lc, tc, D) if lc is not None else None
        if m is None or m["coverage"] < 0.5 or not m["snr"]:
            continue
        err = abs(m["depth"] / m["snr"])
        d = m["depth"]
        ev.append({"n": n, "tc": round(tc, 4), "sector": lc["meta"].get("sector"), "depth_ppm": round(d * 1e6, 1),
                   "err_ppm": round(err * 1e6, 1), "snr": round(m["snr"], 2),
                   "seen": bool(m["snr"] >= 3 and d > 0.4 * dref),
                   "missing": bool(d < dref / 3 and (dref - d) / err > 3)})
    return ev


def counts(ev):
    return {"n_data": len(ev), "n_seen": sum(e["seen"] for e in ev), "n_missing": sum(e["missing"] for e in ev)}


def fold_score(flat, T0, Ps, D):
    """Summed in-transit deficit / error over all predicted transits, for each trial period (vectorised)."""
    t, y, e = flat
    out = np.empty(len(Ps))
    for i, P in enumerate(Ps):
        ph = ((t - T0) / P + 0.5) % 1 - 0.5
        inn = np.abs(ph * P) < 0.4 * D
        if inn.sum() < 2:
            out[i] = -np.inf
            continue
        w = 1 / e[inn] ** 2
        mean = np.sum(w * (1 - y[inn])) / w.sum()
        out[i] = mean * np.sqrt(w.sum())
    return out


def anchors(lcs, ev_listed, T0, D, dref):
    from period_review import measure
    seen = sorted([e for e in ev_listed if e["seen"]], key=lambda e: -e["snr"])
    out = [(e["tc"], "seen transit of the listed ephemeris") for e in seen[:2]]
    if len(seen) >= 3:
        return out
    # few or no transits where the ephemeris says: also look for dips of the right depth, first in the sector of
    # the listed epoch, then in the sectors nearest to it, until three are found
    def gap(l):
        lo, hi = float(np.nanmin(l["time"])), float(np.nanmax(l["time"]))
        return 0.0 if lo <= T0 <= hi else min(abs(lo - T0), abs(hi - T0))
    k = 0
    for lc in sorted(lcs, key=gap)[:6]:
        t = lc["time"]
        found = []
        for tc in np.arange(float(np.nanmin(t)) + D, float(np.nanmax(t)) - D, D / 4):
            m = measure(lc, tc, D)
            if m and m["coverage"] >= 0.7 and m["snr"] >= 6 and 0.4 * dref < m["depth"] < 2.5 * dref:
                found.append((m["snr"], tc))
        found.sort(reverse=True)
        for snr, tc in found:
            if all(abs(tc - o) > 2 * D for o, _ in out):
                out.append((float(tc), f"dip of SNR {snr:.1f} in sector {lc['meta'].get('sector')}"))
                k += 1
            if k == 3:
                return out
    return out


def alternatives(lcs, flat, r, ev_listed, D, dref):
    P0, T0 = float(r["P"]), float(r["T0"])
    Pe = float(r["P_err"]) if r["P_err"] == r["P_err"] else 0.0
    span = max(float(np.nanmax(lc["time"])) for lc in lcs) - min(float(np.nanmin(lc["time"])) for lc in lcs)
    best = None
    tried = []
    for ta, why in anchors(lcs, ev_listed, T0, D, dref):
        for ratio in (1, 0.5, 2, 1 / 3, 3):
            Pc = P0 * ratio
            if Pc < 0.2 or Pc > span:
                continue
            f = min(max(0.005, 5 * Pe * ratio / Pc), 0.05)
            step = max(Pc * (D / 4) / max(span, Pc), 1e-7)
            grid = np.arange(Pc * (1 - f), Pc * (1 + f) + step, step)
            if grid.size > 6000:
                grid = np.linspace(Pc * (1 - f), Pc * (1 + f), 6000)
            sc = fold_score(flat, ta, grid, D)
            top = []
            for i in np.argsort(-sc):
                if all(abs(grid[i] - g) > 3 * step for g in top):
                    top.append(float(grid[i]))
                if len(top) == 2:
                    break
            for P in top:
                ev = series(lcs, ta, P, D, dref)
                c = counts(ev)
                tried.append({"anchor": ta, "why": why, "ratio": ratio, "P": P, **c})
                good = c["n_seen"] >= 3 and c["n_missing"] <= max(1, c["n_data"] // 10)
                if good and (best is None or (c["n_seen"], -c["n_missing"]) > (best["n_seen"], -best["n_missing"])):
                    best = {"anchor": ta, "why": why, "ratio": ratio, "P": P, **c, "events": ev}
    return best, tried


def audit_star(tic, rows, secs=None):
    """All TOIs of one star. rows: list of dicts with toi, P, P_err, T0, T0_err, D, depth, disp."""
    from dossier import all_sectors, best_lcs
    from duo import flatten_sectors
    out = []
    try:
        lcs = [lc for lc in best_lcs(tic, secs if secs is not None else all_sectors(tic)) if lc["time"].size > 100]
    except Exception as ex:
        return [{"tic": tic, "toi": r["toi"], "error": repr(ex)[:200]} for r in rows]
    if not lcs:
        return [{"tic": tic, "toi": r["toi"], "error": "no light curves"} for r in rows]
    flat = None
    for r in rows:
        D = r["D"] / 24
        dref = r["depth"] * 1e-6
        o = {"tic": tic, "toi": r["toi"], "disp": r["disp"], "P": r["P"], "T0": r["T0"], "dur_h": r["D"],
             "depth_ppm": r["depth"], "sectors": sorted({lc["meta"]["sector"] for lc in lcs})}
        try:
            ev = series(lcs, r["T0"], r["P"], D, dref)
            c = counts(ev)
            o.update(c)
            errs = [e["err_ppm"] for e in ev]
            o["detectable"] = bool(errs and r["depth"] >= 3 * float(np.median(errs)))
            o["flag"] = bool(o["detectable"] and c["n_missing"] >= 3 and c["n_missing"] >= c["n_data"] / 3)
            if o["flag"]:
                if flat is None:
                    flat = flatten_sectors(lcs, window=max(2.5, 4 * D))
                best, tried = alternatives(lcs, flat, r, ev, D, dref)
                o["listed_events"] = ev
                o["best"] = best
                o["tried"] = tried
        except Exception as ex:
            o["error"] = repr(ex)[:200]
        out.append(o)
    return out


def main(max_stars=None):
    t = tois()
    rows = {}
    for _, r in t.iterrows():
        rows.setdefault(int(r["TIC ID"]), []).append({
            "toi": f"{float(r['TOI']):.2f}", "P": float(r["Period (days)"]),
            "P_err": float(r["Period (days) err"]) if r["Period (days) err"] == r["Period (days) err"] else float("nan"),
            "T0": float(r["Epoch (BJD)"]) - BJD0,
            "T0_err": float(r["Epoch (BJD) err"]) if r["Epoch (BJD) err"] == r["Epoch (BJD) err"] else float("nan"),
            "D": float(r["Duration (hours)"]), "depth": float(r["Depth (ppm)"]), "disp": str(r["TFOPWG Disposition"])})
    done = set()
    if CKPT.exists():
        for line in open(CKPT):
            try:
                done.add(json.loads(line)["tic"])
            except Exception:
                pass
    todo = [tic for tic in rows if tic not in done]
    if os.environ.get("AUDIT_TICS"):   # test on chosen stars, without the checkpoint
        only = [int(x) for x in os.environ["AUDIT_TICS"].split(",")]
        for tic in only:
            for o in audit_star(tic, rows[tic]):
                b = o.get("best") or {}
                print(json.dumps({k: o.get(k) for k in ("toi", "disp", "P", "n_data", "n_seen", "n_missing", "detectable", "flag")}),
                      "best:", {k: b.get(k) for k in ("P", "ratio", "anchor", "why", "n_data", "n_seen", "n_missing")})
        return
    if max_stars:
        todo = todo[:max_stars]
    print(len(rows), "stars with", len(t), "TOIs;", len(done), "done;", len(todo), "to do", flush=True)
    secs = sectors_for(todo)
    n_proc = int(os.environ.get("N_PROC", "4"))
    with ProcessPoolExecutor(n_proc) as ex, open(CKPT, "a") as ck:
        futs = {ex.submit(audit_star, tic, rows[tic], secs[tic]): tic for tic in todo}
        for i, fu in enumerate(as_completed(futs)):
            tic = futs[fu]
            try:
                res = fu.result()
            except Exception as e:
                res = [{"tic": tic, "toi": r["toi"], "error": repr(e)[:200]} for r in rows[tic]]
            ck.write(json.dumps({"tic": tic, "res": res}, default=float) + "\n")
            ck.flush()
            for o in res:
                if o.get("flag"):
                    b = o.get("best") or {}
                    print(f"TOI-{o['toi']} {o.get('disp')} P={o['P']:.5f} listed seen {o['n_seen']}/{o['n_data']} missing "
                          f"{o['n_missing']} -> best {b.get('P', float('nan')):.5f} (x{b.get('ratio', '-')}) seen "
                          f"{b.get('n_seen', '-')}/{b.get('n_data', '-')} missing {b.get('n_missing', '-')}", flush=True)
            if i % 200 == 0:
                print("progress", i, "/", len(todo), flush=True)
    # merge the checkpoint
    allres = {}
    for line in open(CKPT):
        try:
            d = json.loads(line)
        except Exception:
            continue
        for o in d["res"]:
            allres[f"{o['tic']}|{o['toi']}"] = o
    json.dump(allres, open(OUT, "w"), indent=0, default=float)
    n_j = sum(1 for o in allres.values() if o.get("detectable"))
    n_f = sum(1 for o in allres.values() if o.get("flag"))
    n_b = sum(1 for o in allres.values() if o.get("flag") and o.get("best"))
    print("judged", n_j, "flagged", n_f, "with a better ephemeris", n_b)


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else None)
