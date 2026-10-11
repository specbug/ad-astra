"""Period constraints for the Track B single/duo candidates (results/candidates/trackb_singles.json).

For each star with verdict single or repeating in trackb_final.csv:
  - every TESS sector is loaded (SPOC 2-min > TESS-SPOC > QLP), flattened and binned to 30 min;
  - one event: the shortest period that the other data allow (duo.single_period_constraints);
    two or more events: the allowed aliases dt/n (duo.alias_table) and the shortest of them;
  - the period a circular orbit with b = 0 would need to give the observed duration, from the TIC density
    (an upper bound on the period for non-grazing, circular orbits);
  - a search for a secondary eclipse at phase 0.5 for repeating events with a period (EB test).

  python vet_singles.py [n_threads]
"""
import json
import math
import sys
import concurrent.futures as cf

import env  # noqa: F401
import numpy as np
import pandas as pd

import ticdb
from dossier import all_sectors, best_lcs
from duo import alias_table, flatten_sectors, single_period_constraints
from env import RESULTS

G = 6.674e-11
RHO_SUN = 1408.0
OUT = RESULTS / "candidates"


def p_from_duration(t14_h, rho_solar, k=0.1):
    """Period (d) giving T14 for a central (b=0) circular transit around a star of density rho."""
    # T14 = P/pi * (1+k) / (a/R*),  a/R* = (G rho P^2 / 3pi)^(1/3)  =>  P = (pi T14 / (1+k))^3 * G rho / (3 pi)
    t = t14_h * 3600
    return (math.pi * t / (1 + k)) ** 3 * G * rho_solar * RHO_SUN / (3 * math.pi) / 86400


def one(r):
    tic = int(r["tic"])
    out = {"tic": tic}
    try:
        times = sorted({round(float(x), 4) for x in str(r.get("event_times") or r["t0s"]).split(";") if x})
        depth, dur = float(r["depth"]), float(r["t14_h"])
        st = ticdb.query([tic], ["mass", "rad", "rho"]).iloc[0]
        rho = st["rho"] if st["rho"] == st["rho"] else (st["mass"] / st["rad"] ** 3 if st["mass"] == st["mass"] else 1.0)
        out["rho_solar"] = float(rho)
        out["p_circ_b0"] = p_from_duration(dur, rho, math.sqrt(max(depth, 0)))
        lcs = best_lcs(tic, all_sectors(tic))
        flat = flatten_sectors(lcs)
        out["sectors"] = [lc["meta"]["sector"] for lc in lcs]
        if flat is None:
            return out
        if len(times) == 1:
            grid, allowed, pmin = single_period_constraints(times[0], depth, dur, flat, pmin=2.0, pmax=5000.0)
            out["p_min_allowed"] = pmin
            out["frac_allowed_below_p_circ"] = float(allowed[grid < out["p_circ_b0"]].mean()) if (grid < out["p_circ_b0"]).any() else None
        else:
            rows = alias_table(times, depth, dur, flat, p_min=1.0)
            ok = [a for a in rows if a["status"] == "allowed"]
            out["aliases_allowed"] = [round(a["P"], 4) for a in ok][:200]
            out["n_aliases"] = len(rows)
            out["p_min_allowed"] = min(a["P"] for a in ok) if ok else None
            # secondary eclipse at phase 0.5 of the shortest allowed period (EB test)
            if ok:
                P = out["p_min_allowed"]
                t, y, e = flat
                D = dur / 24
                ph = ((t - times[0]) / P + 0.25) % 1 - 0.25
                sec = np.abs(ph - 0.5) * P < D / 2
                if sec.sum() >= 3:
                    d2 = 1 - np.mean(y[sec])
                    e2 = np.median(e[sec]) / math.sqrt(sec.sum())
                    out["secondary_depth"] = float(d2)
                    out["secondary_sigma"] = float(d2 / e2)
    except Exception as ex:
        out["error"] = repr(ex)[:200]
    return out


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    import report
    tb = report.trackb_shortlist(pd.read_csv(OUT / "trackb_final.csv"))
    res = {}
    with cf.ProcessPoolExecutor(n) as ex:  # processes: the work holds the GIL
        for r in ex.map(one, tb.to_dict("records")):
            res[str(r["tic"])] = r
            print(r["tic"], {k: (round(v, 2) if isinstance(v, float) else v) for k, v in r.items()
                             if k in ("p_min_allowed", "p_circ_b0", "secondary_sigma", "error")}, flush=True)
    json.dump(res, open(OUT / "trackb_singles.json", "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
