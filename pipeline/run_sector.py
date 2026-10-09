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
import classify
from ebsearch import analyse_eb
from env import CATALOGS, RESULTS, SURVEY

HERE = os.path.dirname(os.path.abspath(__file__))
EVENT_COLS = ["tic", "sector", "camera", "ccd", "product", "invert", "t0", "dur_h", "depth", "snr", "local_snr",
              "noise", "coverage", "n_durations", "depth_core", "depth_wing", "depth_first", "depth_second",
              "min_flux", "edge_dist_d", "depth_sap_flux", "depth_det_flux", "depth_flux_sml", "depth_flux_lag",
              "z_bkg", "z_cx", "z_cy", "sys_ccd", "sys_sc", "blend_n", "blend_src", "blend_is_deepest", "clean",
              "n_clean", "snr_max"]


def sample_ids():
    return pd.read_parquet(env.REPO / "results" / "samples" / "survey_sample.parquet", columns=["tic"]).tic.values


def build_jobs(sector, ids):
    df = s3io.read_target_list("qlp", sector)
    tic = pd.to_numeric(df.tic, errors="coerce").dropna().astype("int64")
    tic = tic[tic.isin(pd.Index(ids))]
    return pd.DataFrame({"product": "qlp", "tic": tic.values, "sector": sector})


def ebs(rows):
    def one(r):
        try:
            lc = s3io.load(r["product"], int(r["tic"]), int(r["sector"]))
            if lc is None:
                return None
            res = analyse_eb(lc)
            ex = res.pop("extra_events", [])
            res.update({"tic": int(r["tic"]), "sector": int(r["sector"]), "n_extra": len(ex),
                        "extra": json.dumps(ex)})
            return res
        except Exception as e:
            return {"tic": int(r["tic"]), "sector": int(r["sector"]), "error": repr(e)[:200]}
    with cf.ProcessPoolExecutor(os.cpu_count()) as ex:
        return [x for x in ex.map(one, rows, chunksize=4) if x is not None]


def main():
    sector = int(sys.argv[1])
    n_proc = int(sys.argv[2]) if len(sys.argv) > 2 else os.cpu_count()
    out = RESULTS / "survey" / f"s{sector:04d}"
    if (out / "summary.json").exists():
        print("already done", sector)
        return
    work = SURVEY / f"s{sector:04d}"
    work.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    jobs = build_jobs(sector, sample_ids())
    jobs.to_parquet(work / "jobs.parquet", index=False)
    print(f"s{sector}: {len(jobs)} light curves", flush=True)
    env_vars = dict(os.environ, INVERT_FRAC=os.environ.get("INVERT_FRAC", "0.1"), CHUNK="500")
    subprocess.run([sys.executable, os.path.join(HERE, "survey.py"), str(work / "jobs.parquet"), str(work),
                    str(n_proc), "16"], check=True, env=env_vars)
    t1 = time.time()
    ev, st, ccd_hot, sc_hot = survey_reduce.reduce_dirs([str(work)])
    tic = pd.read_parquet(CATALOGS / "tic_master.parquet")
    tic = tic[tic.tic.isin(ev.tic.unique())]
    ev = classify.lc_context(ev)
    ev = classify.attach_star(ev, tic)
    out.mkdir(parents=True, exist_ok=True)
    # compact event table: clean dips with snr >= 7 (real and flipped), for cross-sector work and FAP
    keep = ev[ev.clean & (ev.snr >= 7)]
    keep[[c for c in EVENT_COLS if c in keep]].to_parquet(out / "events_clean.parquet", index=False)
    b = classify.track_b(ev)
    b.to_parquet(out / "trackb_events.parquet", index=False)
    c = classify.track_c_targets(ev)
    d = classify.track_d(ev)
    d[[c_ for c_ in EVENT_COLS + ["asym", "d_reason", "Tmag", "Teff", "rad"] if c_ in d]].to_parquet(
        out / "trackd_events.parquet", index=False)
    eb = pd.DataFrame(ebs(c.to_dict("records"))) if len(c) else pd.DataFrame()
    eb.to_parquet(out / "trackc_eb.parquet", index=False)
    hot = pd.concat([ccd_hot, sc_hot])
    hot[hot.hot].to_parquet(out / "hot_bins.parquet", index=False)
    # noise and status per light curve (compact) for completeness/rates
    lst = st[~st.invert.astype(bool)]
    lcols = ["tic", "status", "camera", "ccd", "tmag", "n_used", "rms_10m", "sig_1h", "sig_3h", "sig_9h", "n_events"]
    lst[[c_ for c_ in lcols if c_ in lst]].astype({"tic": "int64"}).to_parquet(out / "lcstats.parquet", index=False)
    rates = {}
    for inv in (False, True):
        e = ev[(ev.invert == inv) & ev.clean]
        n_lc = int(((st.invert == inv) & (st.status == "ok")).sum())
        rates[str(inv)] = {"n_lc": n_lc, **{f"snr{x}": int((e.snr >= x).sum()) for x in (7, 8, 10, 12, 15, 20, 30)}}
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
