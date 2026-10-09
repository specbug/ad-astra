"""Define the survey sample for Tracks B/C/D and count light curves per sector.

Sample: TIC stars (objType STAR, not duplicates/artifacts) that are dwarfs or subgiants
(radius <= 2.5 Rsun, or logg >= 3.5 when the radius is missing) with Tmag <= TMAX, plus every
M dwarf (Teff < 3900 K, radius < 0.65 Rsun) down to Tmag 13.5, where even the faint stars show giant-planet
transits several percent deep.

Writes results/samples/survey_sample.parquet (tic only, sorted) and prints per-sector counts.
"""
import sys

import env  # noqa: F401
import numpy as np
import pandas as pd

from env import CATALOGS, RESULTS

TMAX = float(sys.argv[1]) if len(sys.argv) > 1 else 12.0

import pyarrow.compute as pc
import ticdb
tic = ticdb.query(columns=["Tmag", "Teff", "logg", "rad", "objType", "disposition", "lumclass"],
                  expr=(pc.field("Tmag") <= 13.5))
base = (tic.objType == "STAR") & ~tic.disposition.isin(["DUPLICATE", "ARTIFACT", "SPLIT"])
dwarfish = (tic.rad <= 2.5) | (tic.rad.isna() & (tic.logg >= 3.5))
main = base & dwarfish & (tic.Tmag <= TMAX)
mdw = base & (tic.Teff < 3900) & (tic.rad < 0.65) & (tic.Tmag <= 13.5)
s = tic[main | mdw]
print("lumclass of all T<=TMAX stars:", tic[base & (tic.Tmag <= TMAX)].lumclass.value_counts(dropna=False).to_dict())
print(f"sample: {len(s)} stars ({int((main).sum())} dwarfs/subgiants T<={TMAX}, {int(mdw.sum())} M dwarfs T<=13.5)")
out = RESULTS / "samples"
out.mkdir(parents=True, exist_ok=True)
pd.DataFrame({"tic": np.sort(s.tic.values.astype("int64"))}).to_parquet(out / "survey_sample.parquet", index=False,
                                                                       compression="zstd")
q = pd.read_parquet(CATALOGS / "qlp_manifest.parquet")
q = q[q.tic.isin(s.tic)]
per = q.groupby("sector").size()
print("QLP light curves in sample:", len(q))
print(per.to_string())
per.rename("n_lc").to_csv(out / "survey_sample_per_sector.csv")
