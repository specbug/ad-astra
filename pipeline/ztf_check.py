"""Which star eclipses? ZTF light curves (1" pixels) for the Track C doubly eclipsing systems.

TESS pixels are 21" wide, so a TESS light curve blends every star within about an arcminute; the pixel tests
(eb_pixels.py) locate each eclipse set to a few arcsec. ZTF resolves the neighbours directly.

Data: ZTF DR23 light curves, read from IRSA's public HATS/Parquet copy on AWS (ipac-irsa-ztf), so no IRSA
service is needed. Only the position columns of a partition are read to find the objects, then the light curves
of just those rows. ZTF covers dec > -30; stars brighter than about 12.5 mag saturate.

For each system, the target and every Gaia star within 120" bright enough to make either eclipse depth:
  - ZTF objects within 1.5" (one per field, CCD and filter); points with catflags = 0, mag fainter than 12.5,
    each object normalised to its median flux;
  - for each eclipse set (period P, duration from TESS), a box search over P(1 +- 5e-4) with free phase
    (the TESS ephemeris drifts over ZTF's 2018-2024 span), and the same search in 300 random windows of
    periods for the false-alarm probability;
  - the depth this star would need if it were the eclipser: TESS depth x (target flux / star flux).
A star whose light curve shows the set at the right period, significantly and at roughly the needed depth, is the
eclipser; a neighbour with good data and no dip where it would need a deep one is ruled out.

  python ztf_check.py [tic ...]       (default: the 30 Track C systems of eb_pixels.json north of dec -28)
Output: results/candidates/ztf_check.json
"""
import json
import math
import sys

import env  # noqa: F401
import numpy as np
import pandas as pd

import pixels
import ticdb
from env import RESULTS

C = RESULTS / "candidates"
ROOT = "ipac-irsa-ztf/contributed/dr23/lc/hats/ztf_dr23_lc-hats"
_fs = None
_parts = None
_files = {}


def fs():
    global _fs
    if _fs is None:
        import s3fs
        _fs = s3fs.S3FileSystem(anon=True)
    return _fs


def partition_dir(ra, dec):
    global _parts
    import healpy as hp
    if _parts is None:
        p = pd.read_csv(fs().open(ROOT + "/partition_info.csv"))
        _parts = {(int(o), int(n)) for o, n in zip(p.Norder, p.Npix)}
    pix = int(hp.ang2pix(2 ** 9, ra, dec, lonlat=True, nest=True))
    for o in range(9, -1, -1):
        n = pix >> (2 * (9 - o))
        if (o, n) in _parts:
            return f"{ROOT}/dataset/Norder={o}/Dir={(n // 10000) * 10000}/Npix={n}"
    return None


class PartFile:
    """One Parquet file of the light-curve catalogue: positions read once, light curves by row group."""

    def __init__(self, path):
        import pyarrow.parquet as pq
        self.pf = pq.ParquetFile(fs().open(path, block_size=8 << 20))
        t = self.pf.read(columns=["objra", "objdec"])
        self.ra = t.column("objra").to_numpy().astype(float)
        self.dec = t.column("objdec").to_numpy().astype(float)
        md = self.pf.metadata
        self.starts = np.cumsum([0] + [md.row_group(i).num_rows for i in range(md.num_row_groups)])
        self.groups = {}

    def objects(self, ra, dec, r_arcsec):
        d = np.hypot((self.ra - ra) * math.cos(math.radians(dec)), self.dec - dec) * 3600
        out = []
        rows = np.where(d < r_arcsec)[0]
        groups = sorted({int(np.searchsorted(self.starts, i, side="right") - 1) for i in rows})
        for g in groups:   # read a row group, keep only the wanted rows (row groups are ~100 MB in memory)
            tb = self.pf.read_row_group(g, columns=["objectid", "filterid", "objra", "objdec", "lightcurve"])
            out += tb.take([int(i - self.starts[g]) for i in rows if self.starts[g] <= i < self.starts[g + 1]]).to_pylist()
            del tb
        return out


def ztf_objects(ra, dec, r_arcsec=1.5):
    d = partition_dir(ra, dec)
    if d is None:
        return []
    out = []
    for path in sorted(fs().ls(d)):
        if not path.endswith(".parquet"):
            continue
        if path not in _files:
            _files[path] = PartFile(path)
        out += _files[path].objects(ra, dec, r_arcsec)
    return out


def light_curve(objs, bright_limit=12.5):
    """Clean, per-object normalised flux; time in BTJD (HMJD + 2400000.5 - 2457000)."""
    T, F, E, info = [], [], [], []
    for o in objs:
        lc = o["lightcurve"]
        t = np.array(lc["hmjd"], float)
        m = np.array(lc["mag"], float)
        e = np.array(lc["magerr"], float)
        fl = np.array(lc["catflags"], int)
        ok = (fl == 0) & np.isfinite(m) & (m > bright_limit)
        info.append({"objectid": int(o["objectid"]), "filter": int(o["filterid"]), "n": int(t.size), "n_good": int(ok.sum()),
                     "median_mag": float(np.median(m[ok])) if ok.any() else None})
        if ok.sum() < 30:
            continue
        f = 10 ** (-0.4 * (m[ok] - np.median(m[ok])))
        s = 1.4826 * np.median(np.abs(f - 1))
        keep = f > 1 - 50 * s  # keep eclipses, drop wild outliers below
        keep &= f < 1 + 5 * s
        T.append(t[ok][keep] + 2400000.5 - 2457000.0)
        F.append(f[keep])
        E.append(np.maximum(e[ok][keep] * 0.921 * f[keep], s * 0.5))
    if not T:
        return None, info
    t, f, e = (np.concatenate(x) for x in (T, F, E))
    o = np.argsort(t)
    return (t[o], f[o], e[o]), info


def box_scan(t, f, e, P0, D, rel=5e-4, n_null=300, seed=2):
    """Best box depth and SNR over P0(1 +- rel) with free phase; FAP from random period windows."""
    from astropy.timeseries import BoxLeastSquares
    span = t.max() - t.min()
    dP = 0.25 * D * P0 / span
    grid = np.arange(P0 * (1 - rel), P0 * (1 + rel), max(dP, 1e-7))
    bls = BoxLeastSquares(t, f, dy=e)
    r = bls.power(grid, D, objective="snr")
    i = int(np.argmax(r.depth_snr))
    best = {"P": float(grid[i]), "depth": float(r.depth[i]), "snr": float(r.depth_snr[i]), "t0": float(r.transit_time[i])}
    rng = np.random.default_rng(seed)
    null = []
    for c in np.exp(rng.uniform(np.log(0.5), np.log(15.0), n_null)):
        if any(abs(c * k / P0 - 1) < 2 * rel for k in (0.5, 1, 2, 1 / 3, 3)):
            continue
        g = np.arange(c * (1 - rel), c * (1 + rel), max(0.25 * D * c / span, 1e-7))
        null.append(float(np.max(bls.power(g, D, objective="snr").depth_snr)))
    null = np.array(null)
    best["fap"] = float((null >= best["snr"]).mean())
    best["null_p99"] = float(np.percentile(null, 99))
    # sensitivity: 3-sigma upper limit on a dip at the period (phase free), from the noise of the best fit
    best["depth_err"] = float(r.depth_err[i])
    return best


def check(tic, rec):
    st = ticdb.query([tic], ["ra", "dec", "Tmag"]).iloc[0]
    ra, dec, tmag = float(st.ra), float(st.dec), float(st.Tmag)
    sets = {}
    for j in ("1", "2"):
        per = [s["sets"][j] for s in rec["sectors"] if "sets" in s and j in s["sets"]]
        if per:
            sets[j] = {"P": float(np.median([x["P"] for x in per])), "dur_h": float(np.median([x["dur_h"] for x in per])),
                       "depth": float(np.median([x["depth"] for x in per]))}
    if not sets:
        return {"status": "no eclipse parameters"}
    g = pixels.gaia_cone(ra, dec, 120.0)
    dmin = min(s["depth"] for s in sets.values())
    lim = tmag - 2.5 * math.log10(max(dmin, 1e-4))
    srcs = [{"name": "target", "ra": ra, "dec": dec, "tmag": tmag, "sep": 0.0}]
    for _, r in g.iterrows():
        if r.sep_arcsec > 3 and r.tmag_est < lim:
            srcs.append({"name": f"G{r.phot_g_mean_mag:.1f} at {r.sep_arcsec:.0f}\"", "ra": float(r.ra), "dec": float(r.dec),
                         "tmag": float(r.tmag_est), "sep": float(r.sep_arcsec)})
    out = {"status": "checked", "sets": sets, "sources": []}
    _files.clear()   # positions of the partitions of earlier systems are not needed any more
    for s in srcs:
        lc, info = light_curve(ztf_objects(s["ra"], s["dec"]))
        row = dict(s, ztf_objects=info, n_points=int(lc[0].size) if lc else 0)
        if lc is not None:
            t, f, e = lc
            row["rms"] = float(1.4826 * np.median(np.abs(f - 1)))
            for j, p in sets.items():
                need = p["depth"] * 10 ** (0.4 * (s["tmag"] - tmag))
                b = box_scan(t, f, e, p["P"], p["dur_h"] / 24)
                b["needed_depth"] = float(need)
                b["call"] = ("eclipses" if b["snr"] >= 7 and b["fap"] < 0.01 and b["depth"] > 0.2 * min(need, 1) else
                             "ruled out" if need < 1 and b["depth"] + 3 * b["depth_err"] < 0.3 * need else "inconclusive")
                row[f"set{j}"] = b
        out["sources"].append(row)
        print("  ", tic, s["name"], row["n_points"], {j: (row[f"set{j}"]["call"], round(row[f"set{j}"]["snr"], 1),
                                                          round(row[f"set{j}"]["depth"], 4), round(row[f"set{j}"]["needed_depth"], 4))
                                                      for j in sets if f"set{j}" in row}, flush=True)
    # per set: who eclipses
    for j in sets:
        ecl = [r["name"] for r in out["sources"] if r.get(f"set{j}", {}).get("call") == "eclipses"]
        ruled = [r["name"] for r in out["sources"] if r.get(f"set{j}", {}).get("call") == "ruled out"]
        out[f"set{j}_eclipser"] = ecl
        out[f"set{j}_ruled_out"] = ruled
    return out


def main(tics=None):
    rec = json.load(open(C / "eb_pixels.json"))
    pos = pd.read_csv("/home/user/work/trackc30_coords.csv") if tics is None else None
    todo = tics or [int(t) for t in rec if float(pos[pos.tic == int(t)].dec.iloc[0]) > -28]
    p = C / "ztf_check.json"
    res = json.load(open(p)) if p.exists() else {}
    for tic in todo:
        try:
            res[str(tic)] = check(int(tic), rec[str(tic)])
        except Exception as ex:
            res[str(tic)] = {"status": "error", "error": repr(ex)[:200]}
        r = res[str(tic)]
        print(tic, r["status"], {j: (r.get(f"set{j}_eclipser"), len(r.get(f"set{j}_ruled_out", []))) for j in ("1", "2")}, flush=True)
        json.dump(res, open(p, "w"), indent=1, default=float)


if __name__ == "__main__":
    main([int(a) for a in sys.argv[1:]] or None)
