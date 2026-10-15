"""Survey pass: download light curves from S3, flatten them and catalogue every significant dip.

One job = (product, tic, sector). Output per chunk: an events table (one row per dip) and an lcstats
table (one row per light curve: noise levels, coverage, number of dips) used later to recognise
systematics shared by many stars and to calibrate rates.

Usage: python survey.py JOBS.parquet OUTDIR [n_proc] [n_threads]
"""
import os
import sys
import time
import traceback
import concurrent.futures as cf
import multiprocessing as mp

import env  # noqa: F401  (sets thread env vars before numpy import)
import numpy as np
import pandas as pd

import s3io
from lcprep import bin_lc, detrend_fast, clip_upward, segments, robust_std
from dips import find_dips

BIN_MIN = 10.0
WIN_SHORT = 0.75
WIN_LONG = 2.5
INVERT_FRAC = float(os.environ.get("INVERT_FRAC", "0"))  # fraction of light curves also searched flipped


def process_lc(lc, snr_min=6.0, invert=False):
    """Flatten one light curve and catalogue its dips. invert=True runs the identical pipeline on the
    flux flipped about its median (dips become bumps), giving a false-alarm baseline."""
    meta = lc["meta"]
    t, f = lc["time"], lc["flux"].copy()
    flux_like = ("sap_flux", "det_flux", "flux_sml", "flux_lag")
    if invert:
        f = 2 * np.nanmedian(f) - f
        lc = dict(lc)
        for k in flux_like:
            if k in lc:
                lc[k] = 2 * np.nanmedian(lc[k]) - lc[k]
    stats = {k: meta[k] for k in ("tic", "sector", "camera", "ccd", "product", "tmag", "crowdsap")}
    stats["invert"] = bool(invert)
    stats["n_raw"] = int(t.size)
    if t.size < 100:
        stats["status"] = "short"
        return stats, []
    order = np.argsort(t)
    t, f = t[order], f[order]
    aux_raw = {k: lc[k][order] for k in ("bkg", "cx", "cy") + flux_like if k in lc}
    cad = meta["cadence_s"] / 60.0 if np.isfinite(meta["cadence_s"]) else 10.0
    if cad < BIN_MIN * 0.9:
        t, f, cnt, aux = bin_lc(t, f, BIN_MIN, aux_raw)
        cad_d = BIN_MIN / 1440
    else:
        aux = aux_raw
        cad_d = cad / 1440
    # normalise each contiguous segment, drop segments too short to detrend
    keep = np.zeros(t.size, bool)
    for a, b in segments(t, gap=0.5):
        if t[b - 1] - t[a] < 0.75:
            continue
        f[a:b] = f[a:b] / np.nanmedian(f[a:b])
        keep[a:b] = True
    t, f = t[keep], f[keep]
    aux = {k: v[keep] for k, v in aux.items()}
    if t.size < 100:
        stats["status"] = "short"
        return stats, []
    try:
        ys, _ = detrend_fast(t, f, WIN_SHORT, knot_step=WIN_SHORT / 36)
        yl, _ = detrend_fast(t, f, WIN_LONG, knot_step=WIN_LONG / 36)
    except Exception:
        stats["status"] = "detrend_fail"
        return stats, []
    ok = np.isfinite(ys) & np.isfinite(yl) & clip_upward(ys, 4.0)
    t, ys, yl = t[ok], ys[ok], yl[ok]
    aux = {k: v[ok] for k, v in aux.items()}
    for k in flux_like:
        if k in aux:
            aux[k] = aux[k] / np.nanmedian(aux[k])
    segs = segments(t, gap=0.5)
    edges = [t[a] for a, b in segs] + [t[b - 1] for a, b in segs]
    events, sig = find_dips(t, ys, yl, cad_d, snr_min=snr_min, aux=aux, seg_edges=edges)
    stats.update({"status": "ok", "n_used": int(t.size), "span_d": float(t[-1] - t[0]),
                  "n_seg": len(segs), "rms_10m": float(robust_std(ys)),
                  "n_events": len(events), "max_snr": max([e["snr"] for e in events], default=0.0)})
    for D_h, sd in sig.items():
        stats[f"sig_{D_h:g}h"] = sd
    for e in events:
        e.update({"tic": meta["tic"], "sector": meta["sector"], "camera": meta["camera"],
                  "ccd": meta["ccd"], "product": meta["product"], "invert": bool(invert)})
    return stats, events


def _fetch_parse(job):
    product, tic, sector = job
    try:
        key = s3io.key_for(product, tic, sector)
        raw = s3io.fetch(key) if key else None
        if raw is None:
            return job, None, "missing"
        return job, s3io.parse_lc(raw, product), None
    except Exception as ex:  # network or FITS problems: record and move on
        return job, None, f"error:{type(ex).__name__}"


def run_chunk(jobs, n_threads=16):
    stats_rows, event_rows = [], []
    with cf.ThreadPoolExecutor(n_threads) as ex:
        for job, lc, err in ex.map(_fetch_parse, jobs):
            if lc is None:
                stats_rows.append({"product": job[0], "tic": job[1], "sector": job[2], "status": err})
                continue
            modes = [False]
            if INVERT_FRAC > 0 and (job[1] * 2654435761 % 1000) < INVERT_FRAC * 1000:
                modes.append(True)
            for inv in modes:
                try:
                    st, ev = process_lc(lc, invert=inv)
                except Exception:
                    st, ev = {"product": job[0], "tic": job[1], "sector": job[2], "invert": inv,
                              "status": "process_fail:" + traceback.format_exc(limit=1)[-120:]}, []
                stats_rows.append(st)
                event_rows.extend(ev)
    return stats_rows, event_rows


def _worker(args):
    idx, jobs, outdir, n_threads = args
    out_ev = os.path.join(outdir, f"events_{idx:05d}.parquet")
    out_st = os.path.join(outdir, f"lcstats_{idx:05d}.parquet")
    if os.path.exists(out_st):
        return idx, 0, 0.0
    t0 = time.time()
    st, ev = run_chunk(jobs, n_threads)
    pd.DataFrame(ev).to_parquet(out_ev + ".tmp", index=False)
    pd.DataFrame(st).to_parquet(out_st + ".tmp", index=False)
    os.replace(out_ev + ".tmp", out_ev)
    os.replace(out_st + ".tmp", out_st)
    return idx, len(jobs), time.time() - t0


def main():
    jobs_path, outdir = sys.argv[1], sys.argv[2]
    n_proc = int(sys.argv[3]) if len(sys.argv) > 3 else os.cpu_count()
    n_threads = int(sys.argv[4]) if len(sys.argv) > 4 else 12
    chunk = int(os.environ.get("CHUNK", "400"))
    os.makedirs(outdir, exist_ok=True)
    jobs = pd.read_parquet(jobs_path)
    tuples = list(zip(jobs["product"], jobs["tic"].astype(int), jobs["sector"].astype(int)))
    chunks = [(i, tuples[i * chunk:(i + 1) * chunk], outdir, n_threads)
              for i in range((len(tuples) + chunk - 1) // chunk)]
    t0 = time.time()
    done = 0
    with mp.get_context("fork").Pool(n_proc) as pool:
        for idx, n, dt in pool.imap_unordered(_worker, chunks):
            done += n
            el = time.time() - t0
            print(f"chunk {idx} ({n} lcs, {dt:.0f}s); total {done}/{len(tuples)} in {el:.0f}s "
                  f"= {done / max(el, 1):.1f} lc/s", flush=True)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
