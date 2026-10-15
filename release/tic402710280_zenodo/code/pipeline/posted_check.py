"""Periods posted to ExoFOP that the period hunt does not reproduce: which one do the data support?

For the TOIs and community TOIs where someone has posted a period for the same candidate (ExoFOP planet
parameters) and that period is not among the hunt's allowed periods (period_review.json, relaxed with
period_recheck.json), every transit that the posted period predicts in the TESS data is measured, and so is
every transit predicted by the hunt's own best periods. The ephemeris runs through the catalogued transit;
the period is refined within 0.3% of each trial value (the posted values are rounded, and a small period error
grows over hundreds of cycles) by maximising the summed signal of the predicted transits.

Each predicted transit is measured like the review's events (period_review.measure: box depth against a local
quadratic baseline, red-noise error). A period is contradicted when transits it predicts are covered by data and
show no dip (depth below a third of the catalogued depth at more than 3 sigma from it), and supported when the
predicted transits show dips of the right depth. For a period that is half of another, the odd and even transits
are compared.

  python posted_check.py
Output: results/candidates/posted_check.json
"""
import json
import math
import re

import env  # noqa: F401
import numpy as np

from env import RESULTS
from fit_new import lcs_for
from period_review import lc_for, measure

C = RESULTS / "candidates"
OUT = C / "posted_check.json"


def cases():
    rv = json.load(open(C / "period_review.json"))
    rc = json.load(open(C / "period_recheck.json"))
    out = []
    for o in rv:
        if "error" in o or str(o.get("id", "")).startswith("single") or o.get("n_credible", 0) < 1:
            continue
        x = rc.get(f"{o['tic']}|{o['id']}")
        if x and "aliases" in x:
            o["aliases"], o["n_aliases"] = x["aliases"], x["n_aliases"]
        al = o.get("aliases") or []
        al = json.loads(al) if isinstance(al, str) else al
        m = re.search(r"\.(\d\d)$", str(o["id"]))
        if not m:
            continue
        posted = [(name, float(p)) for name, p in re.findall(r"(\S+) P=([0-9.]+)", o.get("exofop_periods") or "")
                  if name.endswith("." + m.group(1)) and (str(o["tic"]) in name or "TOI" in name.upper())]
        if posted and not any(abs(a / p - 1) < 2e-3 for a in al for _, p in posted):
            out.append((o, posted, al))
    return out


def predicted(lcs, T0, P, D):
    """Transit times T0 + nP that fall inside the data (with half a transit of margin)."""
    tcs = []
    for lc in lcs:
        t = lc["time"]
        lo, hi = float(np.nanmin(t)), float(np.nanmax(t))
        n0, n1 = math.ceil((lo + D / 2 - T0) / P), math.floor((hi - D / 2 - T0) / P)
        tcs += [(n, T0 + n * P) for n in range(n0, n1 + 1)]
    return sorted(set(tcs))


def score(lcs, T0, P, D, depth_ref):
    ev = []
    for n, tc in predicted(lcs, T0, P, D):
        lc = lc_for(lcs, tc)
        m = measure(lc, tc, D) if lc is not None else None
        if m is None or m["coverage"] < 0.5:
            continue
        err = m["depth"] / m["snr"] if m["snr"] else np.nan
        ev.append({"n": n, "tc": round(tc, 4), "sector": lc["meta"].get("sector"), "depth_ppm": m["depth"] * 1e6,
                   "err_ppm": abs(err) * 1e6 if err == err else None, "snr": m["snr"], "coverage": m["coverage"]})
    z = [e["depth_ppm"] / e["err_ppm"] for e in ev if e["err_ppm"]]
    total = float(np.sum(z) / math.sqrt(len(z))) if z else 0.0
    dref = depth_ref * 1e6
    for e in ev:
        if e["err_ppm"]:
            e["dip"] = e["snr"] >= 3 and e["depth_ppm"] > 0.4 * dref
            e["no_dip"] = e["depth_ppm"] < dref / 3 and (dref - e["depth_ppm"]) / e["err_ppm"] > 3
    return total, ev


def refine(lcs, T0, P0, D, depth_ref, frac=3e-3):
    """Best period within P0 (1 +- frac): the step keeps the phase drift over the data below D/4."""
    span = max(float(np.nanmax(lc["time"])) for lc in lcs) - T0
    span = max(span, T0 - min(float(np.nanmin(lc["time"])) for lc in lcs))
    step = max(P0 * (D / 4) / max(span, P0), 1e-6)
    grid = np.arange(P0 * (1 - frac), P0 * (1 + frac) + step, step)
    if grid.size > 4000:
        grid = np.linspace(P0 * (1 - frac), P0 * (1 + frac), 4000)
    best = (-np.inf, P0)
    for P in grid:
        tot, _ = score(lcs, T0, P, D, depth_ref)
        if tot > best[0]:
            best = (tot, P)
    return best[1]


def summarise(ev):
    with_data = [e for e in ev if e.get("err_ppm")]
    return {"n_data": len(with_data), "n_dip": sum(1 for e in with_data if e["dip"]),
            "n_no_dip": sum(1 for e in with_data if e["no_dip"])}


def odd_even(ev):
    def wm(x):
        x = [e for e in x if e.get("err_ppm")]
        if not x:
            return None, None
        w = np.array([1 / e["err_ppm"] ** 2 for e in x])
        d = np.array([e["depth_ppm"] for e in x])
        return float((w * d).sum() / w.sum()), float(1 / math.sqrt(w.sum()))
    o = wm([e for e in ev if e["n"] % 2])
    e = wm([e for e in ev if not e["n"] % 2])
    if o[0] is None or e[0] is None:
        return None
    return {"odd_ppm": o[0], "odd_err": o[1], "even_ppm": e[0], "even_err": e[1],
            "diff_sigma": abs(o[0] - e[0]) / math.hypot(o[1], e[1])}


def main():
    res = {}
    for o, posted, al in cases():
        key = f"{o['tic']}|{o['id']}"
        D = float(o["dur_h"]) / 24
        T0 = float(o["ref_epoch"])
        dref = float(o["depth_ref_ppm"]) * 1e-6
        lcs = [lc for lc in lcs_for(int(o["tic"])) if lc["time"].size > 100]
        trials = [("posted " + name, p) for name, p in posted] + [("hunt", a) for a in al[:3]]
        # the halves and doubles of the posted period: a factor of two is the commonest error
        for name, p in posted:
            trials += [("half of posted", p / 2), ("double of posted", 2 * p)]
        out = {"tic": int(o["tic"]), "id": o["id"], "depth_ref_ppm": dref * 1e6, "dur_h": o["dur_h"], "ref_epoch": T0,
               "posted": posted, "hunt_aliases": al[:6], "hunt_n_aliases": o.get("n_aliases"),
               "hunt_events": [(e["sector"], round(e["t"], 3), round(e.get("snr", 0), 1)) for e in o["events"] if e.get("ok")],
               "trials": []}
        for label, P0 in trials:
            P = refine(lcs, T0, P0, D, dref)
            tot, ev = score(lcs, T0, P, D, dref)
            tr = {"label": label, "P0": P0, "P": P, "total_snr": tot, **summarise(ev), "odd_even": odd_even(ev),
                  "events": ev}
            out["trials"].append(tr)
            print(key, f"{label:>22s} P0={P0:9.4f} P={P:9.4f} snr={tot:6.1f} data={tr['n_data']:3d} "
                  f"dip={tr['n_dip']:3d} no_dip={tr['n_no_dip']:3d}", flush=True)
        res[key] = out
        json.dump(res, open(OUT, "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
