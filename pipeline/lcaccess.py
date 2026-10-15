"""Light curves of one star from every source: the public SPOC 2-minute, TESS-SPOC and QLP products on the MAST S3
bucket (best one per sector) and our own photometry of the newest full-frame images (ffi_new.py).
"""
import env  # noqa: F401

import s3io
import ticdb
from env import CATALOGS
from ffi_new import new_lcs


def all_sectors(tic):
    import pandas as pd
    out = {}
    for p in ("spoc", "tspoc", "qlp"):
        m = pd.read_parquet(CATALOGS / f"{p}_manifest.parquet", columns=["tic", "sector"],
                            filters=[("tic", "==", int(tic))])
        out[p] = sorted(m.sector.astype(int).tolist())
    return out


def best_lcs(tic, secs):
    """One light curve per sector: SPOC 2-min > TESS-SPOC > QLP."""
    got, lcs = set(), []
    for p in ("spoc", "tspoc", "qlp"):
        for s in secs.get(p, []):
            if s in got:
                continue
            try:
                lc = s3io.load(p, tic, s)
            except Exception:
                lc = None
            if lc is not None and lc["time"].size > 100:
                lcs.append(lc)
                got.add(s)
    return sorted(lcs, key=lambda l: l["meta"]["sector"])


def lcs_for(tic):
    st = ticdb.query([tic], ["ra", "dec"]).iloc[0]
    lcs = best_lcs(tic, all_sectors(tic))
    have = {lc["meta"]["sector"] for lc in lcs}
    return lcs + [lc for lc in new_lcs(tic, float(st.ra), float(st.dec), build_missing=False, target_time=True) if lc["meta"]["sector"] not in have]
