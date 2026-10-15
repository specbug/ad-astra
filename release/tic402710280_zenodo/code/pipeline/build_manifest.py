"""Build the TIC master table and the per-sector light-curve manifests.

  catalogs/tic_master.parquet      every TIC star with Tmag < 13.5 or d < 200 pc (selected columns)
  catalogs/qlp_manifest.parquet    (tic, sector) for every QLP light curve on S3 of stars in tic_master
  catalogs/tspoc_manifest.parquet  same for TESS-SPOC
"""
import glob
import sys
import concurrent.futures as cf

import env  # noqa: F401
import numpy as np
import pandas as pd
import pyarrow.parquet as pq

import s3io
from env import CATALOGS


def master():
    out = CATALOGS / "tic_master.parquet"
    if out.exists():
        return pd.read_parquet(out)
    parts = sorted(glob.glob(str(CATALOGS / "tic_v2" / "*.parquet")))
    df = pd.concat([pq.read_table(p).to_pandas() for p in parts], ignore_index=True)
    df = df.rename(columns={"ID": "tic"})
    df["tic"] = df["tic"].astype("int64")
    df = df.drop_duplicates("tic")
    df.to_parquet(str(out) + ".tmp", index=False)
    import os
    os.replace(str(out) + ".tmp", out)  # atomic: run_sector polls for this file
    return df


def manifest(product, sectors, ids):
    out = CATALOGS / f"{product}_manifest.parquet"
    if out.exists():
        return pd.read_parquet(out)
    idset = pd.Index(np.unique(ids))

    def one(s):
        df = s3io.read_target_list(product, s)
        if df is None:
            return None
        tic = pd.to_numeric(df["tic"], errors="coerce").dropna().astype("int64")
        tic = tic[tic.isin(idset)]
        return pd.DataFrame({"tic": tic.values, "sector": np.int16(s)})

    with cf.ThreadPoolExecutor(8) as ex:
        res = [r for r in ex.map(one, sectors) if r is not None]
    m = pd.concat(res, ignore_index=True).drop_duplicates()
    m.to_parquet(out, index=False)
    return m


if __name__ == "__main__":
    import ticdb
    ids = ticdb.query(columns=[])["tic"].values
    print("TIC extract rows:", len(ids))
    q = manifest("qlp", range(1, 99), ids)
    print("qlp manifest rows:", len(q), "stars:", q.tic.nunique())
    t = manifest("tspoc", range(1, 82), ids)
    print("tspoc manifest rows:", len(t), "stars:", t.tic.nunique())
