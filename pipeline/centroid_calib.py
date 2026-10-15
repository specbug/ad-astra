"""Calibration of the transit centroid test (centroid_vet.vet) on confirmed planets.

Sample: TOIs that are confirmed or known planets (TFOPWG disposition CP or KP in the ExoFOP TOI table of
12 Oct 2026), T = 8-11, depth 700-10000 ppm, period 2-40 d, with SPOC 2-minute data in at least one sector.
For each, up to four transits from up to two of its 2-minute sectors (predicted from the TOI ephemeris) and
two from a full-frame-only sector when there is one, so that both kinds of pixels are tested.
The fraction called "on target" is the test's true-positive rate on real planets; the fraction put "on" a
neighbour is its false-alarm rate for the report's "off-target" verdicts.

  python centroid_calib.py [n]       (default 40 planets)
Output: results/candidates/centroid_calib.json
"""
import json
import sys

import env  # noqa: F401
import numpy as np
import pandas as pd

from centroid_vet import vet
from dossier import all_sectors
from env import CATALOGS, RESULTS

OUT = RESULTS / "candidates" / "centroid_calib.json"
BTJD0 = 2457000.0


def sample(n, seed=3):
    t = pd.read_csv(CATALOGS / "exofop" / "toi_20261012.csv", low_memory=False)
    P = pd.to_numeric(t["Period (days)"], errors="coerce")
    sel = t[t["TFOPWG Disposition"].isin(["CP", "KP"]) & t["TESS Mag"].between(8, 11) &
            t["Depth (ppm)"].between(700, 10000) & P.between(2, 40)]
    sel = sel.drop_duplicates("TIC ID").sample(frac=1, random_state=seed)
    return sel.head(int(n * 2.5))


def spans():
    return {int(k): v for k, v in json.load(open(CATALOGS / "sector_spans.json")).items()}


def events_for(r, secs, sp):
    P = float(r["Period (days)"])
    T0 = float(r["Epoch (BJD)"]) - BTJD0
    out = []
    spoc = [s for s in secs.get("spoc", []) if s in sp][-2:]
    ffi_only = [s for s in sorted(set(secs.get("qlp", [])) | set(secs.get("tspoc", []))) if s not in secs.get("spoc", []) and s in sp][-1:]
    for s, kmax in [(s, 2) for s in spoc] + [(s, 2) for s in ffi_only]:
        lo, hi = sp[s]
        ks = np.arange(np.ceil((lo + 0.5 - T0) / P), np.floor((hi - 0.5 - T0) / P) + 1)
        for k in ks[:kmax]:
            out.append((s, float(T0 + k * P)))
    return out


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 40
    res = json.load(open(OUT)) if OUT.exists() else {}
    sp = spans()
    for _, r in sample(n).iterrows():
        if len(res) >= n:
            break
        tic = int(r["TIC ID"])
        if str(tic) in res:
            continue
        secs = all_sectors(tic)
        if not secs.get("spoc"):
            continue
        ev = events_for(r, secs, sp)
        if len(ev) < 2:
            continue
        try:
            o = vet(tic, ev, float(r["Duration (hours)"]), float(r["Depth (ppm)"]) * 1e-6)
        except Exception as ex:
            o = {"verdict": "error", "error": repr(ex)[:200]}
        o.update(toi=str(r["TOI"]), name=str(r.get("Planet Name") or ""), tmag=float(r["TESS Mag"]))
        res[str(tic)] = o
        print(tic, r["TOI"], o.get("verdict"), round(o.get("offset_arcsec", -1), 1), flush=True)
        json.dump(res, open(OUT, "w"), indent=1, default=float)
    v = pd.Series([o.get("verdict", "") for o in res.values()])
    print(v.str.slice(0, 12).value_counts().to_string())


if __name__ == "__main__":
    main()
