"""Track A: multi-sector BLS search of every target in a list, real and flipped (false-alarm baseline).

  python run_track_a.py TARGETS.csv OUTDIR [n_proc] [start] [stop]

TARGETS.csv needs columns tic, rad, mass (from track_a_targets.py). For each star every available
sector is loaded from S3 (SPOC 2-min preferred, then TESS-SPOC, then QLP), flattened and stitched, and an
iterative BLS (up to 3 signals) is run on the light curve and on its upside-down copy. One JSON file per
star in OUTDIR/stars/.
"""
import json
import os
import sys
import time
import multiprocessing as mp

import env  # noqa: F401
import numpy as np
import pandas as pd

from env import CATALOGS
from periodic import load_star, prepare_star, bls_search

_MAN = {}


def manifests():
    if not _MAN:
        for p in ("qlp", "tspoc", "spoc"):
            m = pd.read_parquet(CATALOGS / f"{p}_manifest.parquet", columns=["tic", "sector"])
            _MAN[p] = m.groupby("tic").sector.apply(list).to_dict()
    return _MAN


def one(row, outdir):
    tic = int(row["tic"])
    out = os.path.join(outdir, "stars", f"{tic}.json")
    if os.path.exists(out):
        return tic, "cached", 0.0
    t0 = time.time()
    if "sec_qlp" in row:
        sec = {p: [int(float(x)) for x in str(row.get(f"sec_{p}") or "").split(";") if x not in ("", "nan")]
               for p in ("spoc", "tspoc", "qlp")}
    else:
        man = manifests()
        sec = {p: man[p].get(tic, []) for p in ("spoc", "tspoc", "qlp")}
    res = {"tic": tic, "sectors": sec}
    try:
        lcs = load_star(tic, sec)
        lc = prepare_star(lcs)
        if lc is None or lc["time"].size < 500:
            res["status"] = "no_data"
        else:
            rad = float(row.get("rad") or 0.6)
            mass = float(row.get("mass") or rad)
            rho = mass / rad ** 3 if rad > 0 and mass > 0 else 1.0
            res.update({"status": "ok", "n_points": int(lc["time"].size), "bad_sectors": lc["bad_sectors"],
                        "noise_1h": lc["noise_1h"], "rho": rho,
                        "products": sorted({l["meta"]["product"] for l in lcs}),
                        "signals": bls_search(lc, rho, max_signals=2)})
            if tic % 2 == 0:  # flipped-light-curve false-alarm baseline on half of the stars
                res["signals_inverted"] = bls_search(lc, rho, max_signals=1, invert=True)
    except Exception as e:
        res["status"] = f"error: {type(e).__name__}: {str(e)[:200]}"
    res["runtime_s"] = round(time.time() - t0, 1)
    with open(out + ".tmp", "w") as f:
        json.dump(res, f)
    os.replace(out + ".tmp", out)
    return tic, res["status"], res["runtime_s"]


def _star(args):
    return one(*args)


def main():
    targets, outdir = sys.argv[1], sys.argv[2]
    n_proc = int(sys.argv[3]) if len(sys.argv) > 3 else os.cpu_count()
    df = pd.read_csv(targets, dtype={"sec_spoc": str, "sec_tspoc": str, "sec_qlp": str})
    start = int(sys.argv[4]) if len(sys.argv) > 4 else 0
    stop = int(sys.argv[5]) if len(sys.argv) > 5 else len(df)
    df = df.iloc[start:stop]
    os.makedirs(os.path.join(outdir, "stars"), exist_ok=True)
    if "sec_qlp" not in df:
        manifests()  # load once in the parent so forked workers share it
    t0 = time.time()
    with mp.get_context("fork").Pool(n_proc) as pool:
        for i, (tic, status, dt) in enumerate(pool.imap_unordered(_star, [(r, outdir) for r in df.to_dict("records")])):
            print(f"{i + 1}/{len(df)} TIC {tic} {status} {dt}s (elapsed {time.time() - t0:.0f}s)", flush=True)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
