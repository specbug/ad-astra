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


def main(shard=0, n_shards=1):
    """Shards (python period_vet.py i n) split the stars between processes; each writes period_vet_i.json, and
    a run without arguments merges them into period_vet.json."""
    rv = json.load(open(C / "period_review.json"))
    spans = sector_spans()
    res = json.load(open(OUT)) if OUT.exists() else {}
    shards = sorted(C.glob("period_vet_*.json"))
    for p in shards:
        res.update(json.load(open(p)))
    out = OUT if n_shards == 1 else C / f"period_vet_{shard}.json"
    if n_shards == 1 and shards:
        json.dump(res, open(OUT, "w"), indent=1, default=float)
        for p in shards:
            p.unlink()
    mine = json.load(open(out)) if (n_shards > 1 and out.exists()) else res
    for i, o in enumerate(new_periods(rv)):
        key = f"{o['tic']}|{o['id']}"
        if key in res or key in mine or i % n_shards != shard:
            continue
        ev = [(sector_of(o["ref_epoch"], spans), o["ref_epoch"], 1e9)] + \
             [(e["sector"], e["t"], e.get("snr", 0)) for e in o["events"] if e.get("ok")]
        ev = [x for x in ev if x[0] is not None]
        # the four sectors with the strongest transits (one image cube per sector is the cost)
        best = {}
        for s, t, snr in ev:
            best[s] = max(best.get(s, 0), snr)
        keep = set(sorted(best, key=lambda s: -best[s])[:4])
        ev = [(s, t) for s, t, _ in ev if s in keep]
        try:
            r = vet(int(o["tic"]), ev, float(o["dur_h"]), max(float(o["depth_ref_ppm"]) * 1e-6, 1e-4))
        except Exception as ex:
            r = {"verdict": "error", "error": repr(ex)[:200]}
        mine[key] = r
        print(key, r.get("verdict"), round(r.get("offset_arcsec", -1), 1), flush=True)
        json.dump(mine, open(out, "w"), indent=1, default=float)


if __name__ == "__main__":
    import sys
    main(*(int(a) for a in sys.argv[1:3])) if len(sys.argv) > 2 else main()
