"""Memory-light access to the local TIC extract (CATALOGS/tic_v2/*.parquet, one file per HATS partition).

The extract has ~22.5 million rows; loading it whole into pandas needs >15 GB, so every query here
projects columns and filters inside Arrow, one file at a time.
"""
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.dataset as ds

from env import CATALOGS

DIR = CATALOGS / "tic_v2"
READY = DIR / "_SUCCESS"


def dataset():
    return ds.dataset(str(DIR), format="parquet")  # files starting with '_' or '.' are ignored


def query(ids=None, columns=None, expr=None):
    """Rows for the given TIC IDs and/or Arrow filter expression, as pandas with column 'tic'."""
    f = None
    if ids is not None:
        f = pc.field("ID").isin(pa.array(np.unique(np.asarray(ids, dtype="int64"))))
    if expr is not None:
        f = expr if f is None else (f & expr)
    cols = None
    if columns is not None:
        cols = ["ID"] + [c for c in columns if c not in ("ID", "tic")]
    t = dataset().to_table(columns=cols, filter=f)
    df = t.to_pandas()
    return df.rename(columns={"ID": "tic"}).astype({"tic": "int64"})


STAR_COLS = ["ra", "dec", "Tmag", "Teff", "logg", "rad", "mass", "rho", "lumclass", "d", "contratio"]


def stars(ids):
    return query(ids, STAR_COLS)
