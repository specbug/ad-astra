"""Combine per-sector survey outputs (results/survey/sNNNN/) into cross-sector candidate lists.

  python aggregate.py [stage]      stage in {b, c, d, all}

Track B: stage-2 vet every planet-sized single dip on a dwarf (real and flipped), group the planet-like
         ones by star across sectors, attach known-TOI matches, write results/candidates/trackb_*.
Track C: EBs whose template-subtracted residuals show another periodic eclipse signal; keep signals
         that are not harmonics of the main period and, when the star has several sectors, recur at the
         same period; write results/candidates/trackc_*.
Track D: stage-2 fit of anomalous dips; write results/candidates/trackd_*.
"""
import glob
import json
import sys

import env  # noqa: F401
import numpy as np
import pandas as pd

import known
import stage2
import ticdb
from env import RESULTS

OUT = RESULTS / "candidates"


def load(name):
    files = sorted(glob.glob(str(RESULTS / "survey" / "s*" / name)))
    dfs = []
    for f in files:
        d = pd.read_parquet(f)
        if len(d):
            d["product"] = "spoc" if "_spoc" in f else "qlp"
            dfs.append(d)
    return pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()


def track_b(n_proc=None):
    ev = load("trackb_events.parquet")
    if ev.empty:
        print("no Track B events yet")
        return None
    sel = ev[ev.planet_sized & ev.dwarf & (ev.snr >= 8)].copy()
    cache = OUT / "trackb_stage2.parquet"
    done = pd.read_parquet(cache) if cache.exists() else pd.DataFrame(columns=["tic", "sector", "t0", "invert"])
    key = ["tic", "sector", "t0", "invert"]
    todo = sel.merge(done[key], on=key, how="left", indicator=True)
    todo = todo[todo._merge == "left_only"].drop(columns="_merge")
    print(f"Track B: {len(ev)} events, {len(sel)} planet-sized dwarf events (snr>=8), {len(todo)} new to vet")
    if len(todo):
        res = stage2.classify_shape(stage2.run(todo, n_proc))
        done = pd.concat([done, res], ignore_index=True) if len(done) else res
        OUT.mkdir(parents=True, exist_ok=True)
        done.to_parquet(cache, index=False)
    done = stage2.classify_shape(done)
    n_lc = {}
    for f in glob.glob(str(RESULTS / "survey" / "s*" / "summary.json")):
        s = json.load(open(f))
        for k in ("False", "True"):
            n_lc[k] = n_lc.get(k, 0) + s["rates"][k]["n_lc"]
    for inv in (False, True):
        x = done[done.invert == inv]
        n = n_lc.get(str(inv), 1)
        print(f"  invert={inv}: {n} LCs, vetted {len(x)}, planet-like {int(x.planet_like.sum())} "
              f"({1e3 * x.planet_like.sum() / max(n, 1):.3f} per 1000 LCs)")
    pl = done[done.planet_like & ~done.invert].copy()
    # group by star
    rows = []
    for tic, g in pl.groupby("tic"):
        g = g.sort_values("t0")
        rows.append({"tic": int(tic), "n_events": len(g), "sectors": ";".join(str(int(s)) for s in g.sector),
                     "t0s": ";".join(f"{t:.4f}" for t in g.fit_t0), "depth": g.fit_depth.median(),
                     "t14_h": g.fit_t14_h.median(), "ingress": g.fit_ingress_frac.median(),
                     "snr_max": g.snr_here.max(), "rp_rearth": g.fit_rp_rearth.median(), "Tmag": g.Tmag.iloc[0],
                     "Teff": g.Teff.iloc[0], "rad": g.rad.iloc[0], "contratio": g.contratio.iloc[0]})
    cand = pd.DataFrame(rows)
    if len(cand):
        st = ticdb.query(cand.tic.values, ["ra", "dec", "d", "logg", "mass"])
        cand = cand.merge(st, on="tic", how="left")
        cand = known.flag_known(cand)
        cand = cand.sort_values(["is_toi", "snr_max"], ascending=[True, False])
    OUT.mkdir(parents=True, exist_ok=True)
    cand.to_csv(OUT / "trackb_candidates.csv", index=False, float_format="%.6g")
    print(f"  planet-like stars: {len(cand)}, of which known TOIs: {int(cand.is_toi.sum()) if len(cand) else 0}")
    return cand


def track_c():
    eb = load("trackc_eb.parquet")
    if eb.empty or "P2" not in eb:
        print("no Track C results yet")
        return None
    e = eb[eb.P2.notna()].copy()
    r = e.P2 / e.P1
    harm = np.zeros(len(e), bool)
    for a in (1, 2, 3, 4, 5, 6):
        for b in (1, 2, 3, 4, 5, 6):
            harm |= np.abs(r - a / b) < 0.02 * a / b
    e["harmonic_any"] = harm
    e["near_P1"] = np.abs(np.log(r)) < np.log(1.5)  # residual of an imperfect EB template, usually
    strong = e[(e.snr2 >= 15) & ~e.harmonic_any & (e.n_tr2 >= 2)].copy()
    # eclipse-likeness (box depth / sinusoid peak-to-trough) is missing from older worker outputs:
    # recompute it for the strong second signals by re-running the EB analysis
    if "P2_box_over_sine" not in strong or strong.P2_box_over_sine.isna().any():
        import concurrent.futures as cf
        from run_sector import _eb_one
        need = strong[strong.get("P2_box_over_sine", pd.Series(np.nan, index=strong.index)).isna()]
        rows = [{"product": r["product"], "tic": r.tic, "sector": r.sector} for _, r in need.iterrows()]
        with cf.ProcessPoolExecutor(4) as ex:
            res = list(ex.map(_eb_one, rows, chunksize=2))
        bos = {(x["tic"], x["sector"]): (x.get("P2_box_over_sine"), x.get("P2"), x.get("snr2")) for x in res if x}
        strong["P2_box_over_sine"] = [bos.get((r.tic, r.sector), (np.nan,))[0] if np.isnan(r.get("P2_box_over_sine", np.nan))
                                      else r.P2_box_over_sine for _, r in strong.iterrows()]
    strong = strong[strong.P2_box_over_sine > 3]
    rows = []
    for tic, g in strong.groupby("tic"):
        g = g.sort_values("sector")
        # recurrence: same P2 (2%) in >= 2 sectors
        p2 = g.P2.values
        rec = max(int(np.sum(np.abs(p2 / p - 1) < 0.02)) for p in p2)
        rows.append({"tic": int(tic), "n_sectors": len(g), "sectors": ";".join(str(int(s)) for s in g.sector),
                     "P1": g.P1.median(), "depth1": g.depth1.median(), "P2": g.P2.iloc[np.argmax(g.snr2.values)],
                     "depth2": g.depth2.median(), "snr2_max": g.snr2.max(), "P2_recurs": rec,
                     "near_P1": bool(g.near_P1.all()), "P3": g.P3.iloc[np.argmax(g.snr2.values)] if "P3" in g else np.nan,
                     "snr3": g.snr3.max() if "snr3" in g else np.nan,
                     "box_over_sine": g.P2_box_over_sine.median()})
    cand = pd.DataFrame(rows)
    if len(cand):
        st = ticdb.query(cand.tic.values, ["ra", "dec", "Tmag", "Teff", "rad", "d"])
        cand = known.flag_known(cand.merge(st, on="tic", how="left"))
        cand = cand.sort_values(["P2_recurs", "snr2_max"], ascending=False)
    OUT.mkdir(parents=True, exist_ok=True)
    cand.to_csv(OUT / "trackc_candidates.csv", index=False, float_format="%.6g")
    print(f"Track C: {len(eb)} EB analyses, {len(e)} with residual signal, {len(strong)} strong non-harmonic, "
          f"{len(cand)} stars; recurring in >=2 sectors: {int((cand.P2_recurs >= 2).sum()) if len(cand) else 0}")
    return cand


def track_d(n_proc=None):
    d = load("trackd_events.parquet")
    if d.empty:
        print("no Track D events yet")
        return None
    d["invert"] = False
    sel = d[d.d_reason.isin(["long_deep", "asymmetric"]) | (d.d_reason == "aperiodic_multi")]
    res = stage2.run(sel, n_proc)
    OUT.mkdir(parents=True, exist_ok=True)
    res.to_parquet(OUT / "trackd_stage2.parquet", index=False)
    print(f"Track D: {len(d)} events, vetted {len(res)}; by reason {res.d_reason.value_counts().to_dict()}")
    return res


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    if which in ("b", "all"):
        track_b()
    if which in ("c", "all"):
        track_c()
    if which in ("d", "all"):
        track_d()
