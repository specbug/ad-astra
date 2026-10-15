"""List SPOC 2-min light curves on S3 for every sector -> (tic, sector, has_fast)."""
import re, concurrent.futures as cf

import env, s3io
import pandas as pd
from env import CATALOGS
pat = re.compile(r"-s(\d{4})-(\d{16})-\d{4}-(s_lc|a_fast-lc)\.fits$")

def one_prefix(pfx):
    rows = []
    for k, _ in s3io.list_keys(pfx):
        m = pat.search(k)
        if m:
            rows.append((int(m.group(2)), int(m.group(1)), m.group(3) == "a_fast-lc"))
    return rows

def sector(s):
    # split listing by the first TIC-path level to parallelise
    tops = [p for p, _ in s3io.list_keys(f"tess/public/tid/s{s:04d}/", delimiter="/")]
    subs = []
    for t in tops:
        subs += [p for p, _ in s3io.list_keys(t, delimiter="/")]
    rows = []
    with cf.ThreadPoolExecutor(16) as ex:
        for r in ex.map(one_prefix, subs or tops):
            rows += r
    df = pd.DataFrame(rows, columns=["tic", "sector", "fast"])
    return df[df.sector == s]

out = []
with cf.ThreadPoolExecutor(4) as ex:
    for df in ex.map(sector, range(1, 108)):
        out.append(df); print(df.sector.iloc[0] if len(df) else None, len(df), flush=True)
m = pd.concat(out, ignore_index=True)
m = m.groupby(["tic", "sector"], as_index=False).fast.max()
m.to_parquet(CATALOGS / "spoc_manifest.parquet", index=False)
print("DONE", len(m), m.tic.nunique())
