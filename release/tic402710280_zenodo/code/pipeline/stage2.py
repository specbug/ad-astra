"""Stage 2 for single-dip candidates (Track B, and Track D long/asymmetric dips).

  python stage2.py OUT.parquet EVENTS.parquet [EVENTS.parquet ...]

Reads candidate events (survey output tables), re-fetches each light curve from S3 and runs
vet_events.vet_one (also on flipped-light-curve events, so the cuts can be calibrated), then applies
the shape/uniqueness classification below.
"""
import os
import sys
import concurrent.futures as cf

import env  # noqa: F401
import numpy as np
import pandas as pd

import s3io
from vet_events import vet_one

R_EARTH_PER_RSUN = 109.1


def _one(r):
    try:
        lc = s3io.load(r.get("product", "qlp") or "qlp", int(r["tic"]), int(r["sector"]))
        if lc is None:
            return {"status": "missing"}
        if r.get("invert"):
            med = np.nanmedian(lc["flux"])
            lc["flux"] = 2 * med - lc["flux"]
            if "sap_flux" in lc:
                lc["sap_flux"] = 2 * np.nanmedian(lc["sap_flux"]) - lc["sap_flux"]
        out = vet_one(lc, float(r["t0"]), float(r["dur_h"]))
        out["status"] = "ok"
        return out
    except Exception as e:
        return {"status": f"error:{type(e).__name__}:{str(e)[:80]}"}


def run(events, n_proc=None):
    rows = events.to_dict("records")
    with cf.ProcessPoolExecutor(n_proc or os.cpu_count()) as ex:
        res = list(ex.map(_one, rows, chunksize=8))
    return pd.concat([events.reset_index(drop=True), pd.DataFrame(res)], axis=1)


def classify_shape(df):
    """Planet-like single transits: flat-bottomed, symmetric, unique, present in raw flux, planet-sized."""
    d = df.copy()
    rs = d.rad.fillna(1.0)
    d["fit_rp_rearth"] = np.sqrt(d.fit_depth.clip(lower=0) * (1 + d.contratio.fillna(0))) * rs * R_EARTH_PER_RSUN
    snr = d.snr_here
    d["c_unique_bump"] = d.max_bump_snr < 0.5 * snr
    d["c_unique_dip"] = (d.max_other_dip_snr < 0.5 * snr) | (d.n_clean >= 2)  # a 2nd similar dip may be a repeat
    d["c_shape_u"] = d.fit_ingress_frac <= 0.5
    d["c_symmetric"] = d.dchi2_asym < np.maximum(10, 0.02 * d.dchi2_trap)
    d["c_sap"] = (d.sap_fit_depth / d.fit_depth).between(0.6, 1.6)
    d["c_signif"] = d.dchi2_trap > 100
    d["c_planet"] = d.fit_rp_rearth < 2.2 * 11.2
    d["c_fitok"] = d.redchi2_trap < 2.0
    crit = [c for c in d.columns if c.startswith("c_")]
    d["n_pass"] = d[crit].sum(axis=1)
    d["planet_like"] = d[crit].all(axis=1)
    return d


if __name__ == "__main__":
    out = sys.argv[1]
    ev = pd.concat([pd.read_parquet(p) for p in sys.argv[2:]], ignore_index=True)
    res = classify_shape(run(ev))
    res.to_parquet(out, index=False)
    for inv in (False, True):
        x = res[res.invert == inv]
        print(f"invert={inv}: {len(x)} events, planet_like {int(x.planet_like.sum())}; "
              f"criteria pass rates: " + ", ".join(f"{c}={x[c].mean():.2f}" for c in x.columns if c.startswith('c_')))
