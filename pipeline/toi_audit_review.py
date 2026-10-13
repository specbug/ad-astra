"""Second stage of the TOI ephemeris audit (toi_audit.py): corrected ephemerides with errors, checked against
what ExoFOP already has.

For every flagged TOI, the alternatives tried by the audit are ranked by the fraction of covered transits they
see (a shorter period predicts more transits, so the count alone would favour it), preferring the listed period
over its multiples when they do equally well. For the chosen ephemeris:
  - each seen transit is timed: the box centre within +-D/2 of the prediction that maximises the box SNR, with an
    error of D / (2 SNR), floored at half a cadence;
  - a weighted straight-line fit to the times gives P and T0 with errors (scaled by sqrt(chi2_red) when above 1);
  - odd and even transits are compared (a period that is half the true one shows two depths);
  - the listed and the corrected ephemerides are compared on 1 December 2026 (BTJD 4100.5): the offset in hours
    and in transit durations, and the listed ephemeris's own error there;
  - ExoFOP's planet-parameter uploads for the star (target pages, cached) are searched for a period within 3
    sigma of the corrected one, i.e. someone who has already posted it.

  N_PROC=4 python toi_audit_review.py
Output: results/candidates/toi_audit_corrections.json and .csv
"""
import json
import math
import os
from concurrent.futures import ProcessPoolExecutor, as_completed

import env  # noqa: F401
import numpy as np
import pandas as pd

from env import CATALOGS, RESULTS

C = RESULTS / "candidates"
CKPT = CATALOGS.parent / "toi_audit.jsonl"
OUT = C / "toi_audit_corrections.json"
REF = 4100.5          # BTJD of 2026-12-01 00:00 UTC (BJD 2461100.5)


def flagged():
    out = {}
    for line in open(CKPT):
        try:
            d = json.loads(line)
        except Exception:
            continue
        for o in d["res"]:
            if o.get("flag"):
                out[f"{o['tic']}|{o['toi']}"] = o
    return out


def choose(o):
    """The shortest period the data do not contradict: its multiples see a subset of the same transits, so a higher
    seen fraction alone would favour them. Shorter periods must see at least 80% of the best fraction (half the
    true period sees about half its transits)."""
    good = [t for t in o.get("tried", []) if t["n_seen"] >= 3 and t["n_missing"] <= max(1, t["n_data"] // 10)]
    if not good:
        return None
    frac = lambda t: t["n_seen"] / max(t["n_data"], 1)
    top = max(frac(t) for t in good)
    ok = [t for t in good if frac(t) >= 0.8 * top]
    return min(ok, key=lambda t: (round(t["P"], 4), -t["n_seen"]))


def time_transit(lcs, tc0, D):
    from period_review import lc_for, measure
    lc = lc_for(lcs, tc0)
    if lc is None:
        return None
    grid = tc0 + np.linspace(-D / 2, D / 2, 41)
    best = None
    for tc in grid:
        m = measure(lc, tc, D)
        if m and m["coverage"] >= 0.5 and (best is None or m["snr"] > best[1]["snr"]):
            best = (tc, m)
    if best is None or best[1]["snr"] < 3:
        return None
    cad = float(np.median(np.diff(lc["time"])))
    return float(best[0]), max(D / (2 * best[1]["snr"]), cad / 2), best[1]["depth"], best[1]["snr"]


def review_one(key, o):
    import toi_audit as ta
    from dossier import all_sectors, best_lcs
    tic = int(o["tic"])
    ch = choose(o)
    res = {"key": key, "tic": tic, "toi": o["toi"], "disp": o.get("disp"), "P_listed": o["P"], "T0_listed": o["T0"],
           "dur_h": o["dur_h"], "depth_ppm": o["depth_ppm"], "listed": {k: o.get(k) for k in ("n_data", "n_seen", "n_missing")}}
    if ch is None:
        res["status"] = "no alternative"
        return res
    D = o["dur_h"] / 24
    dref = o["depth_ppm"] * 1e-6
    lcs = [lc for lc in best_lcs(tic, all_sectors(tic)) if lc["time"].size > 100]
    ev = ta.series(lcs, ch["anchor"], ch["P"], D, dref)
    seen = [e for e in ev if e["seen"]]
    times = []
    for e in seen:
        x = time_transit(lcs, e["tc"], D)
        if x:
            times.append((e["n"], *x))
    res.update(ratio=ch["ratio"], anchor=ch["anchor"], why=ch["why"], chosen={k: ch[k] for k in ("n_data", "n_seen", "n_missing")})
    if len(times) < 3:
        res["status"] = "too few timed transits"
        return res
    n = np.array([t[0] for t in times], float)
    tt = np.array([t[1] for t in times])
    st = np.array([t[2] for t in times])
    n0 = np.round(np.average(n, weights=1 / st ** 2))
    A = np.vstack([np.ones_like(n), n - n0]).T
    W = np.diag(1 / st ** 2)
    cov = np.linalg.inv(A.T @ W @ A)
    p = cov @ A.T @ W @ tt
    chi2 = float(((tt - A @ p) / st) ** 2 @ np.ones_like(tt))
    dof = max(len(tt) - 2, 1)
    scale = math.sqrt(max(chi2 / dof, 1.0))
    T0, P = float(p[0]), float(p[1])
    eT0, eP = float(math.sqrt(cov[0, 0]) * scale), float(math.sqrt(cov[1, 1]) * scale)
    res.update(P=P, P_err=eP, T0=T0, T0_err=eT0, n_timed=len(tt), chi2_red=chi2 / dof,
               depth_seen_ppm=float(np.median([t[3] for t in times]) * 1e6))
    # odd / even
    dep = np.array([t[3] for t in times])
    odd = (n.astype(int) % 2) == 1
    if odd.sum() >= 2 and (~odd).sum() >= 2:
        do, de = dep[odd].mean(), dep[~odd].mean()
        so = dep[odd].std(ddof=1) / math.sqrt(odd.sum()) if odd.sum() > 1 else np.nan
        se = dep[~odd].std(ddof=1) / math.sqrt((~odd).sum()) if (~odd).sum() > 1 else np.nan
        res["odd_even_sigma"] = float(abs(do - de) / math.hypot(so, se)) if so == so and se == se and (so or se) else None
    # where the two ephemerides put the transit on 1 Dec 2026
    k_c = round((REF - T0) / P)
    t_c = T0 + k_c * P
    k_l = round((t_c - o["T0"]) / o["P"])
    t_l = o["T0"] + k_l * o["P"]
    Pe_l = o.get("P_err")
    res.update(offset_h=float((t_l - t_c) * 24), offset_dur=float(abs(t_l - t_c) / D),
               listed_err_h=float(abs(k_l) * Pe_l * 24) if Pe_l else None,
               corrected_err_h=float(math.hypot(eT0, k_c * eP) * 24))
    le = res["listed_err_h"]
    res["sigma_listed"] = float(abs(res["offset_h"]) / le) if le else None
    if res["offset_dur"] < 1:
        res["status"] = "within a transit duration"
    elif le and abs(res["offset_h"]) < 3 * le:
        res["status"] = "refined"          # the listed error already allows the offset; ours is much smaller
    else:
        res["status"] = "corrected"        # the listed ephemeris is wrong by more than its own error
    res["solid"] = bool(res["chosen"]["n_seen"] >= 5 and res["chosen"]["n_seen"] >= 0.7 * res["chosen"]["n_data"]
                        and res["chi2_red"] < 5)
    return res


def exofop_uploads(tic):
    """Periods posted on the ExoFOP target page (cached JSON; fetched once if missing)."""
    from exofop_check import fetch
    try:
        d = fetch(tic)
    except Exception:
        return []
    out = []
    for e in d.get("planet_parameters", []) if isinstance(d, dict) else []:
        try:
            per = float(e.get("per"))
        except (TypeError, ValueError):
            continue
        if per > 0:
            out.append({"name": e.get("name") or e.get("toi"), "per": per, "per_e": e.get("per_e"), "epoch": e.get("epoch"),
                        "user": e.get("puser"), "date": str(e.get("pdate"))[:10]})
    return out


def main():
    import toi_audit as ta
    fl = flagged()
    t = ta.tois()
    perr = {f"{float(x):.2f}": (float(a) if a == a else None) for x, a in zip(t["TOI"], t["Period (days) err"])}
    for o in fl.values():
        o["P_err"] = perr.get(o["toi"])
    res = json.load(open(OUT)) if OUT.exists() else {}
    todo = {k: o for k, o in fl.items() if k not in res}
    print(len(fl), "flagged;", len(todo), "to review", flush=True)
    n_proc = int(os.environ.get("N_PROC", "4"))
    with ProcessPoolExecutor(n_proc) as ex:
        futs = {ex.submit(review_one, k, o): k for k, o in todo.items()}
        for i, fu in enumerate(as_completed(futs)):
            k = futs[fu]
            try:
                r = fu.result()
            except Exception as e:
                r = {"key": k, "status": "error", "error": repr(e)[:200]}
            res[k] = r
            if i % 25 == 0:
                json.dump(res, open(OUT, "w"), indent=1, default=float)
                print(i, "/", len(todo), flush=True)
    # ExoFOP uploads for the corrected ones (one request per star, cached)
    for k, r in res.items():
        if r.get("status") in ("corrected", "refined") and "exofop" not in r:
            ups = exofop_uploads(r["tic"])
            r["exofop"] = ups
            near = [u for u in ups if abs(u["per"] - r["P"]) < max(3 * r["P_err"], 1e-4 * r["P"])]
            r["already_posted"] = [f"{u['user']} {u['date']} P={u['per']}" for u in near]
    json.dump(res, open(OUT, "w"), indent=1, default=float)
    rows = [{k: v for k, v in r.items() if not isinstance(v, (list, dict))} | {"already_posted": "; ".join(r.get("already_posted", []))}
            for r in res.values()]
    pd.DataFrame(rows).to_csv(C / "toi_audit_corrections.csv", index=False)
    st = pd.Series([r.get("status") for r in res.values()]).value_counts()
    print(st.to_string())


if __name__ == "__main__":
    main()
