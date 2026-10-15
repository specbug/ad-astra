"""Track B final vetting: star-level checks for every planet-like single-dip star that is not a known TOI.

For each candidate star:
  1. all-sector dip search (starcheck): deeper eclipses elsewhere => eclipsing binary (reject);
     consistent repeats => duo/periodic candidate, with a linear-ephemeris fit when >= 3 events
  2. SPOC TCEs from DV reports (2-min targets): did NASA's pipeline already detect it?
  3. TOIs within 150" with a matching period (signal assigned to a neighbour)
Writes results/candidates/trackb_final.csv (+ .json with event lists).

  python finalize_b.py [n_proc]
"""
import json
import sys
import concurrent.futures as cf

import env  # noqa: F401
import numpy as np
import pandas as pd

import known
from starcheck import check_star
from env import RESULTS

OUT = RESULTS / "candidates"


def ephemeris(times, dur_h, pmin=1.0, tol_frac=0.25):
    """Longest period consistent with all event times.

    Candidate periods are pair separations divided by small integers; each is refined by a linear
    least-squares fit t = t0 + k*P, and accepted if every event lies within tol_frac*duration of the
    fitted ephemeris. Shorter aliases (P/m) also fit by construction; those are tested against the
    data separately (duo.alias_table). Returns (P, rms residual in minutes) or None.
    """
    t = np.unique(np.round(np.sort(np.asarray(times, float)), 4))
    if t.size < 2:
        return None
    D = dur_h / 24
    cands = []
    for i in range(t.size):
        for j in range(i + 1, t.size):
            dt = t[j] - t[i]
            for n in range(1, int(dt / pmin) + 1):
                cands.append(dt / n)
    best = None
    for P in sorted(set(np.round(cands, 6)), reverse=True):
        k = np.round((t - t[0]) / P)
        if len(np.unique(k)) < t.size:
            continue  # two events assigned to the same epoch
        A = np.vstack([np.ones_like(k), k]).T
        (t0, Pf), *_ = np.linalg.lstsq(A, t, rcond=None)
        res = t - (t0 + k * Pf)
        if np.max(np.abs(res)) < tol_frac * D:
            best = (float(Pf), float(np.sqrt(np.mean(res ** 2)) * 1440))
            break
    return best


def one(row):
    tic = int(row["tic"])
    t0s = [float(x) for x in str(row["t0s"]).split(";") if x]
    out = {"tic": tic}
    try:
        res, _ = check_star(tic, float(row["depth"]), float(row["t14_h"]), t0s)
        out.update({k: res[k] for k in ("n_sectors", "n_deeper", "n_similar")})
        out["sectors_all"] = ";".join(str(s) for s in res["sectors"])
        times = sorted(set([round(t, 3) for t in t0s] + [round(e["t0"], 3) for e in res["similar"]]))
        out["event_times"] = ";".join(f"{t:.3f}" for t in times)
        eph = ephemeris(times, float(row["t14_h"])) if len(times) >= 2 else None
        out["period_guess"] = eph[0] if eph else np.nan
        out["eph_resid_min"] = eph[1] if eph else np.nan
        out["deeper_examples"] = json.dumps(res["deeper"][:5])
    except Exception as e:
        out["error"] = repr(e)[:200]
    try:
        tc = known.spoc_tces(tic)
        out["spoc_tces"] = json.dumps([(x["run"], round(x["period"], 4), round(x.get("depth_ppm", 0)), round(x["mes"], 1))
                                       for x in tc])
        out["n_spoc_tce"] = len(tc)
    except Exception:
        out["n_spoc_tce"] = -1
    if np.isfinite(row.get("ra", np.nan)):
        p = out.get("period_guess")
        near = known.nearby_tois(row["ra"], row["dec"], 150, p if p and np.isfinite(p) else None)
        out["nearby_toi_match"] = json.dumps(near[:3], default=str)
        out["n_nearby_toi_match"] = len(near)
    return out


def main():
    n_proc = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    c = pd.read_csv(OUT / "trackb_candidates.csv")
    c = c[~c.is_toi].copy()
    done = pd.read_csv(OUT / "trackb_final.csv") if (OUT / "trackb_final.csv").exists() else pd.DataFrame(columns=["tic"])
    todo = c[~c.tic.isin(done.tic)]
    print(f"{len(c)} non-TOI planet-like stars; {len(todo)} to finalize")
    # work in chunks and write after each one, so an interrupted run resumes where it stopped
    allr = done
    recs = todo.to_dict("records")
    chunk = 150
    for i0 in range(0, len(recs), chunk):
        # processes, not threads: the per-star work (numba dip search, pandas) holds the GIL
        with cf.ProcessPoolExecutor(n_proc) as ex:
            rows = list(ex.map(one, recs[i0:i0 + chunk]))
        new = c.merge(pd.DataFrame(rows), on="tic", how="inner") if rows else pd.DataFrame()
        allr = pd.concat([allr, new], ignore_index=True) if len(allr) else new
        if len(allr):
            allr["verdict"] = np.select(
                [allr.n_deeper > 0, allr.n_nearby_toi_match.fillna(0) > 0, allr.n_similar.fillna(0) >= 1],
                ["EB (deeper eclipses)", "TOI on neighbour", "repeating"], "single")
            allr.to_csv(OUT / "trackb_final.csv", index=False, float_format="%.6g")
        print(f"  finalized {min(i0 + chunk, len(recs))}/{len(recs)}", flush=True)
    if len(allr):
        print(allr.verdict.value_counts().to_dict())
    return allr


if __name__ == "__main__":
    main()
