"""Pixel test (difference-image centroids, centroid_vet.vet) of the transits behind every new period that the
period review (period_review.py) narrows to three or fewer choices and that is on no list.

Input: results/candidates/period_review.json. Output: results/candidates/period_vet.json, keyed "tic|id".

  python period_vet.py
"""
import glob
import json
import re

import env  # noqa: F401
import numpy as np

from centroid_vet import vet
from env import CATALOGS, RESULTS

C = RESULTS / "candidates"
OUT = C / "period_vet.json"


def sector_spans():
    sp = {int(k): v for k, v in json.load(open(CATALOGS / "sector_spans.json")).items()}
    for f in glob.glob(str(CATALOGS.parent / "cache" / "ffi_lc" / "*_s*.npz")):
        s = int(re.search(r"_s(\d+)\.npz$", f).group(1))
        if s in sp and s < 99:
            continue
        t = np.load(f)["time"]
        lo, hi = float(t.min()), float(t.max())
        sp[s] = [min(lo, sp.get(s, [lo, hi])[0]), max(hi, sp.get(s, [lo, hi])[1])]
    return sp


def sector_of(t, spans):
    for s, (lo, hi) in spans.items():
        if lo - 0.5 <= t <= hi + 0.5:
            return s
    return None


def new_periods(rv):
    def known(o):
        return re.findall(r"P=([0-9.]+)", (o.get("nasa") or "") + (o.get("exofop_periods") or ""))
    return [o for o in rv if "error" not in o and o.get("n_credible", 0) >= 1 and not known(o)
            and o.get("n_aliases") and o["n_aliases"] <= 3 and not str(o.get("id", "")).startswith("single")]


def main():
    rv = json.load(open(C / "period_review.json"))
    spans = sector_spans()
    res = json.load(open(OUT)) if OUT.exists() else {}
    for o in new_periods(rv):
        key = f"{o['tic']}|{o['id']}"
        if key in res:
            continue
        ev = [(sector_of(o["ref_epoch"], spans), o["ref_epoch"])] + [(e["sector"], e["t"]) for e in o["events"] if e.get("ok")]
        ev = [(s, t) for s, t in ev if s is not None]
        try:
            r = vet(int(o["tic"]), ev, float(o["dur_h"]), max(float(o["depth_ref_ppm"]) * 1e-6, 1e-4))
        except Exception as ex:
            r = {"verdict": "error", "error": repr(ex)[:200]}
        res[key] = r
        print(key, r.get("verdict"), round(r.get("offset_arcsec", -1), 1), flush=True)
        json.dump(res, open(OUT, "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
