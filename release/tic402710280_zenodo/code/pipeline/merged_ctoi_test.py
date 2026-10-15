"""Are the three community TOIs that ExoFOP merged into TOIs really separate signals? (Independent audits, 15 October
2026: equal depths could mean one object with a wrong TOI ephemeris rather than two.)

For each star every event is measured in the same way: the community TOI's own dips, and every transit the TOI's
catalogue ephemeris predicts inside the data. Each event gets a trapezoid fit (toi_audit_timing.time_one) for its
mid-time and depth, at the duration that fits it best (scanned from 0.4 to 2.5 times the catalogue duration).
Then two tests:
  - shape: do the community TOI's events and the TOI's transits share a duration and depth?
  - one object: is there any period, from 0.5 d up to the shortest gap between events, that puts every measured
    event on one linear ephemeris (each within 3 sigma, sigma including 5 min of timing systematics) and is not
    contradicted by data (a predicted transit, at least 60% covered, measured at less than a third of the depth
    and more than 4 sigma below it)? Every such period is listed.

  python merged_ctoi_test.py
Output: results/candidates/merged_ctoi_test.json
"""
import json
import math

import env  # noqa: F401
import numpy as np
import pandas as pd

from env import CATALOGS, RESULTS
from toi_audit_timing import time_one

OUT = RESULTS / "candidates" / "merged_ctoi_test.json"
CASES = [
    # star, community TOI, its event times (BTJD) and duration (h), TOI
    (156716001, "CTOI 156716001.01", [1485.406, 2224.520, 3702.770], 8.9, "7602.01"),
    (369532319, "CTOI 369532319.01", [1755.810, 2831.560], 5.8, "5967.01"),
    (67711670, "CTOI 67711670.01", [1798.985, 3622.800], 6.4, "6033.01"),
]
SYS_MIN = 5.0


def best_duration(lc, tc, D0):
    best = None
    for fD in (0.4, 0.55, 0.7, 0.85, 1.0, 1.2, 1.4, 1.7, 2.0, 2.5):
        m = time_one(lc, tc, D0 * fD)
        if m and not m.get("edge") and (best is None or m["dchi2_flat"] > best[1]["dchi2_flat"]):
            best = (D0 * fD, m)
    if best is None:
        return None
    D, m = best
    return {"tc": m["tc"], "err_min": m["err"] * 1440, "depth_ppm": m["depth"] * 1e6, "depth_err_ppm": m["depth_err"] * 1e6,
            "snr": m["snr"], "dur_h": D * 24, "dchi2": m["dchi2_flat"]}


def lc_at(lcs, t):
    for lc in lcs:
        if lc["time"].min() - 0.05 <= t <= lc["time"].max() + 0.05:
            return lc
    return None


def box_depth(lc, tc, D):
    t, f = lc["time"], lc["flux"]
    w = np.isfinite(f) & (np.abs(t - tc) < D / 2 + 0.4)
    inn = w & (np.abs(t - tc) < 0.4 * D)
    oot = w & (np.abs(t - tc) > 0.6 * D)
    cad = float(np.median(np.diff(t)))
    if inn.sum() < 3 or oot.sum() < 15:
        return None
    c = np.polyfit(t[oot] - tc, f[oot], 2)
    r = f / np.polyval(c, t - tc)
    sig = 1.4826 * np.median(np.abs(np.diff(r[oot]))) / math.sqrt(2)
    return {"coverage": min(1.0, inn.sum() * cad / (0.8 * D)), "depth": float(1 - r[inn].mean()),
            "err": float(sig / math.sqrt(inn.sum()))}


def one_object(lcs, ev, D, depth):
    """Periods that put every event on one ephemeris and are not contradicted by flat data."""
    tt = np.array([e["tc"] for e in ev])
    st = np.hypot(np.array([e["err_min"] for e in ev]), SYS_MIN) / 1440
    o = np.argsort(tt)
    tt, st = tt[o], st[o]
    gaps = np.diff(tt)
    gmin = gaps.min()
    cands = []
    # every period that fits the shortest gap an integer number of times, down to 0.5 d
    for m in range(1, int(gmin / 0.5) + 1):
        P0 = gmin / m
        n = np.round((tt - tt[0]) / P0)
        A = np.vstack([np.ones_like(n), n]).T
        W = np.diag(1 / st ** 2)
        cov = np.linalg.inv(A.T @ W @ A)
        p = cov @ A.T @ W @ tt
        res = (tt - A @ p) / st
        if np.all(np.abs(res) < 3):
            cands.append((float(p[1]), float(p[0]), float(np.sum(res ** 2))))
    out = []
    for P, T0, chi2 in cands:
        miss, seen, covered = 0, 0, 0
        for lc in lcs:
            t = lc["time"]
            for k in range(math.ceil((t.min() - T0) / P), math.floor((t.max() - T0) / P) + 1):
                tc = T0 + k * P
                if np.min(np.abs(tt - tc)) < 0.2:
                    continue
                b = box_depth(lc, tc, D)
                if not b or b["coverage"] < 0.6:
                    continue
                covered += 1
                if b["depth"] < depth / 3 and (depth - b["depth"]) > 4 * b["err"]:
                    miss += 1
                elif b["depth"] > 3 * b["err"]:
                    seen += 1
        out.append({"P": P, "T0": T0, "chi2": chi2, "other_covered": covered, "other_missing": miss, "other_seen": seen})
    return out


def main():
    from fit_new import lcs_for
    toi = pd.read_csv(CATALOGS / "exofop" / "toi_20261012.csv", low_memory=False)
    toi["key"] = toi["TOI"].map(lambda x: f"{float(x):.2f}")
    res = {}
    for tic, name, ctimes, cdur, toi_id in CASES:
        r = toi[toi.key == toi_id].iloc[0]
        P_t, T0_t = float(r["Period (days)"]), float(r["Epoch (BJD)"]) - 2457000
        D_t, dep_t = float(r["Duration (hours)"]) / 24, float(r["Depth (ppm)"])
        lcs = [lc for lc in lcs_for(tic) if lc["time"].size > 100]
        ctoi_ev = []
        for tc in ctimes:
            lc = lc_at(lcs, tc)
            m = best_duration(lc, tc, cdur / 24) if lc else None
            if m:
                ctoi_ev.append({**m, "sector": lc["meta"].get("sector"), "group": "community TOI"})
        toi_ev, toi_flat = [], []
        for lc in lcs:
            t = lc["time"]
            for k in range(math.ceil((t.min() - T0_t) / P_t), math.floor((t.max() - T0_t) / P_t) + 1):
                tc = T0_t + k * P_t
                m = best_duration(lc, tc, D_t)
                if m and m["snr"] >= 3 and abs(m["tc"] - tc) < 0.5 * D_t:
                    toi_ev.append({**m, "sector": lc["meta"].get("sector"), "cycle": k, "group": "TOI"})
                else:
                    b = box_depth(lc, tc, D_t)
                    if b and b["coverage"] >= 0.6:
                        toi_flat.append({"tc": tc, "sector": lc["meta"].get("sector"), "depth_ppm": b["depth"] * 1e6,
                                         "err_ppm": b["err"] * 1e6})

        def summ(ev):
            if not ev:
                return None
            d = np.array([e["dur_h"] for e in ev])
            dp = np.array([e["depth_ppm"] for e in ev])
            de = np.array([e["depth_err_ppm"] for e in ev])
            w = 1 / de ** 2
            return {"n": len(ev), "dur_h_median": float(np.median(d)), "dur_h_range": [float(d.min()), float(d.max())],
                    "depth_ppm": float(np.sum(w * dp) / w.sum()), "depth_err_ppm": float(1 / math.sqrt(w.sum()))}
        allev = ctoi_ev + toi_ev
        D_one = float(np.median([e["dur_h"] for e in allev])) / 24 if allev else D_t
        dep_one = float(np.median([e["depth_ppm"] for e in allev])) * 1e-6 if allev else dep_t * 1e-6
        single = one_object(lcs, allev, D_one, dep_one) if len(allev) >= 2 else []
        ok = [s for s in single if s["other_missing"] == 0]
        res[name] = {"tic": tic, "toi": f"TOI-{toi_id}", "toi_catalogue": {"P": P_t, "T0_btjd": T0_t, "dur_h": D_t * 24, "depth_ppm": dep_t},
                     "ctoi_events": ctoi_ev, "toi_events": toi_ev, "toi_predicted_flat": toi_flat,
                     "ctoi_shape": summ(ctoi_ev), "toi_shape": summ(toi_ev),
                     "single_object_periods_tested": len(single), "single_object_allowed": ok[:20],
                     "single_object_best_rejected": sorted(single, key=lambda s: s["other_missing"])[:5]}
        print(name, "| CTOI", summ(ctoi_ev), "| TOI", summ(toi_ev), "| TOI predicted but flat:", len(toi_flat),
              "| one-object periods fitting all times:", len(single), "not contradicted:", len(ok),
              [round(s["P"], 4) for s in ok[:8]], flush=True)
    json.dump(res, open(OUT, "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
