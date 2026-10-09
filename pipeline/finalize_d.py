"""Track D final vetting: anomalous dips (asymmetric / very long) that survive stage 2.

Selection after stage 2 (results/candidates/trackd_stage2.parquet):
  - contamination ratio < 0.5 (otherwise the dip may belong to a brighter neighbour)
  - dip SNR >= 10 and the strongest same-duration bump < half of it
  - asymmetric: dchi2(asymmetric vs symmetric trapezoid) > max(30, 10% of the transit dchi2) with T14 > 2 h
    (shorter "asymmetries" are 10-min binning artefacts of fast eclipses), or
  - long: T14 >= 15 h and depth >= 0.5%
Then every sector of the star is searched: stars with strictly periodic, repeating dips of the same depth
(eclipsing binaries) are dropped unless the selected dip itself is asymmetric.

  python finalize_d.py [n_proc]
"""
import json
import sys
import concurrent.futures as cf

import env  # noqa: F401
import numpy as np
import pandas as pd

import ticdb
from starcheck import check_star
from finalize_b import ephemeris
from env import RESULTS

OUT = RESULTS / "candidates"


def select(d):
    d = d[d.status == "ok"].copy()
    st = ticdb.query(d.tic.unique(), ["contratio", "ra", "dec"])
    d = d.drop(columns=[c for c in ("contratio",) if c in d]).merge(st, on="tic", how="left")
    base = (d.contratio.fillna(0) < 0.5) & (d.snr_here >= 10) & (d.max_bump_snr < 0.5 * d.snr_here)
    asym = (d.dchi2_asym > np.maximum(30, 0.1 * d.dchi2_trap)) & (d.fit_t14_h > 2)
    long_ = (d.fit_t14_h >= 15) & (d.fit_depth >= 0.005)
    d["d_kind"] = np.select([asym, long_], ["asymmetric", "long"], "")
    return d[base & (d.d_kind != "")]


def one(r):
    out = {"tic": int(r["tic"]), "sector": int(r["sector"]), "t0": float(r["t0"])}
    try:
        res, _ = check_star(int(r["tic"]), float(r["fit_depth"]), float(r["fit_t14_h"]), [float(r["t0"])])
        out.update({k: res[k] for k in ("n_sectors", "n_deeper", "n_similar")})
        times = [float(r["t0"])] + [e["t0"] for e in res["similar"]]
        eph = ephemeris(times, float(r["fit_t14_h"])) if len(times) >= 3 else None
        out["periodic_all"] = bool(eph is not None and eph[1] < 30)
        out["similar"] = json.dumps(res["similar"][:8])
    except Exception as e:
        out["error"] = repr(e)[:200]
    return out


def main():
    n_proc = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    d = select(pd.read_parquet(OUT / "trackd_stage2.parquet"))
    print("Track D selected:", len(d), d.d_kind.value_counts().to_dict())
    with cf.ThreadPoolExecutor(n_proc) as ex:
        rows = list(ex.map(one, d.to_dict("records")))
    f = d.merge(pd.DataFrame(rows), on=["tic", "sector", "t0"], how="left")
    keep = (f.d_kind == "asymmetric") | ~f.periodic_all.fillna(False).astype(bool)
    f["verdict"] = np.where(keep & (f.n_deeper.fillna(0) == 0), "anomaly candidate",
                            np.where(f.n_deeper.fillna(0) > 0, "EB (deeper dips)", "periodic (EB/planet)"))
    f.to_csv(OUT / "trackd_final.csv", index=False, float_format="%.6g")
    print(f.verdict.value_counts().to_dict())
    return f


if __name__ == "__main__":
    main()
