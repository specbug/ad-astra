"""Alias test redone with tolerances suited to weak, short transits, for the period review's new periods and its
validation cases (targets whose period is already on record).

period_review.py takes the duration from the box that maximises the signal-to-noise near the catalogued transit
and requires every event within half a duration of the alias grid. For short, weak transits both are too strict:
a 1-hour box makes a 50-minute timing error fatal, and the true period can be excluded (CTOI 243417115.01: a period
posted by another user, 34.951 d, fits dips at all four predicted times, but the 1-hour box mis-timed one by 52
minutes). Here the box is at least 1.5 h long, events may sit up to max(half the box, 1 h) from the grid, and
predicted transits are tested at half the measured depth.

  python period_recheck.py
Output: results/candidates/period_recheck.json, keyed "tic|id": the relaxed list of allowed periods.
"""
import concurrent.futures as cf
import json
import re

import env  # noqa: F401

from env import RESULTS

C = RESULTS / "candidates"
OUT = C / "period_recheck.json"


def known(o):
    return re.findall(r"P=([0-9.]+)", (o.get("nasa") or "") + (o.get("exofop_periods") or ""))


def recheck(o):
    from duo import alias_table, flatten_sectors
    from fit_new import lcs_for
    try:
        lcs = lcs_for(int(o["tic"]))
        D_h = max(float(o["dur_h"]), 1.5)
        good = sorted([o["ref_epoch"]] + [e["t"] for e in o["events"] if e.get("ok")])
        flat = flatten_sectors(lcs, window=max(2.5, 4 * D_h / 24))
        rows = alias_table(good, 0.5 * o["depth_ref_ppm"] * 1e-6, D_h, flat, p_min=1.0,
                           tol=max(0.5 * D_h / 24, 1 / 24))
        al = [a["P"] for a in rows if a["status"] == "allowed"]
        return f"{o['tic']}|{o['id']}", {"aliases": [round(p, 5) for p in al], "n_aliases": len(al), "dur_box_h": D_h}
    except Exception as ex:
        return f"{o['tic']}|{o['id']}", {"error": repr(ex)[:200]}


def main():
    rv = json.load(open(C / "period_review.json"))
    todo = [o for o in rv if "error" not in o and o.get("n_credible", 0) >= 1 and not str(o.get("id", "")).startswith("single")
            and (known(o) or (o.get("n_aliases") and o["n_aliases"] <= 3))]
    res = json.load(open(OUT)) if OUT.exists() else {}
    todo = [o for o in todo if f"{o['tic']}|{o['id']}" not in res]
    print(len(todo), "to recheck", flush=True)
    with cf.ProcessPoolExecutor(4) as ex:
        for i, (k, v) in enumerate(ex.map(recheck, todo)):
            res[k] = v
            if i % 20 == 0:
                print(i, k, v.get("n_aliases"), flush=True)
                json.dump(res, open(OUT, "w"), indent=1)
    json.dump(res, open(OUT, "w"), indent=1)
    print("done", len(res))


if __name__ == "__main__":
    main()
