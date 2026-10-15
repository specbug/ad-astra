"""Extract TIC rows with Tmag<13.5 or d<200 pc from the HATS TIC on S3 (column-projected range reads)."""
import sys, os, io, csv, time, urllib.request, concurrent.futures as cf

import env, s3parquet
import pyarrow.parquet as pq, pyarrow.compute as pc

OUT = str(env.CATALOGS / "tic_v2")
COLS = ["ID", "ra", "dec", "Tmag", "Teff", "logg", "rad", "mass", "rho", "lumclass", "d", "plx",
        "contratio", "numcont", "GAIA", "HIP", "objType", "disposition", "wdflag", "priority",
        "MH", "Vmag", "Kmag", "gaiabp", "gaiarp", "GAIAmag"]
os.makedirs(OUT, exist_ok=True)
parts = list(csv.reader(io.StringIO(urllib.request.urlopen(
    "https://stpubdata.s3.amazonaws.com/tess/public/hats/tic/tic/partition_info.csv").read().decode())))[1:]

def do(p):
    o, n = int(p[0]), int(p[1])
    out = f"{OUT}/o{o}_n{n}.parquet"
    if os.path.exists(out): return 0
    key = f"tess/public/hats/tic/tic/dataset/Norder={o}/Dir={(n // 10000) * 10000}/Npix={n}.parquet"
    for attempt in range(4):
        try:
            tab = s3parquet.read_columns(key, COLS, threads=12)
            m = pc.or_(pc.fill_null(pc.less(tab["Tmag"], 13.5), False), pc.fill_null(pc.less(tab["d"], 200.0), False))
            tab = tab.filter(m)
            pq.write_table(tab, out + ".tmp", compression="zstd"); os.replace(out + ".tmp", out)
            return tab.num_rows
        except Exception as e:
            err = e; time.sleep(5)
    print("FAILED", key, repr(err)[:200], flush=True); return -1

t0 = time.time(); tot = 0; done = 0; failed = 0
with cf.ThreadPoolExecutor(int(sys.argv[1]) if len(sys.argv) > 1 else 6) as ex:
    for r in ex.map(do, parts):
        done += 1; tot += max(r, 0); failed += (r < 0)
        if done % 25 == 0: print(f"{done}/{len(parts)} partitions, {tot} rows kept, {time.time()-t0:.0f}s", flush=True)
print("DONE", done, tot, f"{time.time()-t0:.0f}s", "failed", failed, flush=True)
if failed == 0:
    open(os.path.join(OUT, "_SUCCESS"), "w").write("ok\n")
else:
    sys.exit(1)
