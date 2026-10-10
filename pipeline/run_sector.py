"""Run the full survey for one QLP sector and write compact, committable results.

  python run_sector.py SECTOR [n_proc]

Steps: job list (QLP target list of the sector, intersected with the survey sample) -> survey.py
(dip catalogue, with 10% of light curves also searched flipped) -> survey_reduce (systematics, blends,
quality flags) -> classify (Track B/C/D) -> ebsearch on Track C light curves.
Outputs go to results/survey/sNNNN/ (small parquet files); bulky intermediates stay in WORK and are
deleted at the end unless KEEP=1.
"""
import json
import os
import shutil
import subprocess
import sys
import time
import concurrent.futures as cf

import env  # noqa: F401
import numpy as np
import pandas as pd

import s3io
import survey_reduce
import ticdb
import classify
from ebsearch import analyse_eb
from env import CATALOGS, RESULTS, SURVEY

HERE = os.path.dirname(os.path.abspath(__file__))
EVENT_COLS = ["tic", "sector", "camera", "ccd", "product", "invert", "t0", "dur_h", "depth", "snr", "local_snr",
              "noise", "coverage", "n_durations", "depth_core", "depth_wing", "depth_first", "depth_second",
              "min_flux", "edge_dist_d", "depth_sap_flux", "depth_det_flux", "depth_flux_sml", "depth_flux_lag",
              "z_bkg", "z_cx", "z_cy", "sys_ccd", "sys_sc", "blend_n", "blend_src", "blend_is_deepest", "clean",
              "n_clean", "snr_max"]


B_COLS = ["tic", "sector", "camera", "ccd", "invert", "t0", "dur_h", "depth", "snr", "local_snr", "noise", "coverage",
          "n_durations", "depth_core", "depth_wing", "depth_first", "depth_second", "edge_dist_d", "depth_sap_flux",
          "depth_det_flux", "depth_flux_sml", "depth_flux_lag", "z_bkg", "z_cx", "z_cy", "blend_n", "n_clean",
          "snr_max", "Tmag", "Teff", "logg", "rad", "contratio", "rp_rearth", "p_dur_d", "shape_vu", "dwarf",
          "planet_sized", "asym", "d_reason", "periodic"]


def slim(df, cols=B_COLS):
    df = df[[c for c in cols if c in df]].copy()
    for col in df.select_dtypes("float64").columns:
        if col != "t0":
            df[col] = df[col].astype("float32")
    return df


def sample_ids():
    return pd.read_parquet(env.REPO / "results" / "samples" / "survey_sample.parquet", columns=["tic"]).tic.values


PRODUCT = os.environ.get("PRODUCT", "qlp")


def build_jobs(sector, ids):
    if PRODUCT == "spoc":  # every 2-minute target of the sector (they were pre-selected; keep them all)
        m = pd.read_parquet(CATALOGS / "spoc_manifest.parquet", columns=["tic", "sector"])
        tic = m[m.sector == sector].tic.astype("int64")
        return pd.DataFrame({"product": "spoc", "tic": tic.values, "sector": sector})
    df = s3io.read_target_list("qlp", sector)
    tic = pd.to_numeric(df.tic, errors="coerce").dropna().astype("int64")
    tic = tic[tic.isin(pd.Index(ids))]
    return pd.DataFrame({"product": "qlp", "tic": tic.values, "sector": sector})


def _eb_one(r):
    try:
        lc = s3io.load(r["product"], int(r["tic"]), int(r["sector"]))
        if lc is None:
            return None
        res = analyse_eb(lc)
        ex = res.pop("extra_events", [])
        res.update({"tic": int(r["tic"]), "sector": int(r["sector"]), "n_extra": len(ex), "extra": json.dumps(ex)})
        return res
    except Exception as e:
        return {"tic": int(r["tic"]), "sector": int(r["sector"]), "error": repr(e)[:200]}


def ebs(rows):
    with cf.ProcessPoolExecutor(os.cpu_count()) as ex:
        return [x for x in ex.map(_eb_one, rows, chunksize=4) if x is not None]


def main():
    sector = int(sys.argv[1])
    n_proc = int(sys.argv[2]) if len(sys.argv) > 2 else os.cpu_count()
    out = RESULTS / "survey" / (f"s{sector:04d}" + ("" if PRODUCT == "qlp" else f"_{PRODUCT}"))
    if (out / "summary.json").exists():
        print("already done", sector)
        return
    work = SURVEY / (f"s{sector:04d}" + ("" if PRODUCT == "qlp" else f"_{PRODUCT}"))
    work.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    jobs = build_jobs(sector, sample_ids())
    if os.environ.get("LIMIT"):  # testing: a random subset of the sector
        jobs = jobs.sample(min(len(jobs), int(os.environ["LIMIT"])), random_state=0)
    jobs.to_parquet(work / "jobs.parquet", index=False)
    print(f"s{sector}: {len(jobs)} light curves", flush=True)
    env_vars = dict(os.environ, INVERT_FRAC=os.environ.get("INVERT_FRAC", "0.1"), CHUNK="500")
    subprocess.run([sys.executable, os.path.join(HERE, "survey.py"), str(work / "jobs.parquet"), str(work),
                    str(n_proc), "16"], check=True, env=env_vars)
    t1 = time.time()
    # the TIC master table may still be building in the background (worker.sh); wait up to 3 h
    for _ in range(360):
        if ticdb.READY.exists():
            break
        time.sleep(30)
    ev, st, ccd_hot, sc_hot = survey_reduce.reduce_dirs([str(work)])
    ok_frac = float((st.status == "ok").mean()) if len(st) else 0.0
    if ok_frac < 0.5:  # e.g. a product re-released under a new file name: fail loudly, never write empty results
        raise RuntimeError(f"only {ok_frac:.1%} of light curves loaded: {st.status.value_counts().to_dict()}")
    tic = ticdb.stars(ev.tic.unique())
    ev = classify.lc_context(ev)
    ev = classify.attach_star(ev, tic)
    out.mkdir(parents=True, exist_ok=True)
    # compact table of strong clean dips (real and flipped) for later cross-sector work
    keep = ev[ev.clean & (ev.snr >= 10)]
    keep = keep[[c for c in EVENT_COLS if c in keep]].copy()
    for col in keep.select_dtypes("float64").columns:
        if col not in ("t0",):
            keep[col] = keep[col].astype("float32")
    keep.to_parquet(work.parent / f"events_strong_s{sector:04d}.parquet", index=False, compression="zstd")  # local only
    b = classify.track_b(ev)
    slim(b).to_parquet(out / "trackb_events.parquet", index=False, compression="zstd")
    c = classify.track_c_targets(ev)
    d = classify.track_d(ev)
    slim(d).to_parquet(out / "trackd_events.parquet", index=False, compression="zstd")
    eb = pd.DataFrame(ebs(c.to_dict("records"))) if len(c) else pd.DataFrame()
    if "n_extra" in eb:
        eb.loc[eb.n_extra == 0, "extra"] = ""
    slim(eb, list(eb.columns)).to_parquet(out / "trackc_eb.parquet", index=False, compression="zstd")
    hot = pd.concat([ccd_hot, sc_hot])
    hot[hot.hot].to_parquet(out / "hot_bins.parquet", index=False)
    # noise vs magnitude summary (per-light-curve stats are too bulky to commit)
    lst = st[(~st.invert.astype(bool)) & (st.status == "ok")].copy()
    lst["tbin"] = (lst.tmag * 2).round() / 2
    noise = lst.groupby("tbin").agg(n=("tic", "size"), **{f"{k}_med": (k, "median") for k in
                                    ("sig_1h", "sig_3h", "sig_9h", "rms_10m") if k in lst})
    noise.reset_index().to_csv(out / "noise_vs_tmag.csv", index=False, float_format="%.6g")
    rates = {}
    for inv in (False, True):
        e = ev[(ev.invert == inv) & ev.clean]
        n_lc = int(((st.invert == inv) & (st.status == "ok")).sum())
        rates[str(inv)] = {"n_lc": n_lc, **{f"snr{x}": int((e.snr >= x).sum()) for x in (7, 8, 10, 12, 15, 20, 30)},
                           "by_dur_snr8": {f"{d:g}": int(((e.dur_h == d) & (e.snr >= 8)).sum())
                                           for d in sorted(e.dur_h.unique())}}
    summary = {"sector": sector, "n_jobs": int(len(jobs)), "status": st.status.value_counts().to_dict(),
               "n_events": int(len(ev)), "n_clean": int(ev.clean.sum()), "n_trackb": int(len(b)),
               "n_trackc": int(len(c)), "n_trackd": int(len(d)), "rates": rates,
               "survey_s": round(t1 - t0), "total_s": round(time.time() - t0)}
    (out / "summary.json").write_text(json.dumps(summary, indent=1, default=int))
    print(json.dumps(summary, default=int), flush=True)
    if os.environ.get("KEEP") != "1":
        shutil.rmtree(work, ignore_errors=True)


if __name__ == "__main__":
    main()
