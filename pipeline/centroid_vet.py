"""Where does each transit come from? Difference-image centroids for the transit candidates.

The standard pixel test (as in the TESS data validation reports): for every transit, the out-of-transit image
(a linear fit to the pixels just outside the transit) minus the in-transit image shows only the light that
disappeared; a Gaussian PSF with free position, width and amplitude is fitted to it (centroids._centroid). The
offsets are measured from where the target really is on each image: the out-of-transit image is fitted with
Gaussians at the Gaia positions moved to the observing date (centroids._calibrate), which removes errors of the
image's coordinate solution and the target's proper motion. They are combined over all transits of a star
(3-sigma clipped, error from their scatter, with a 2" floor or the scatter between sectors) and compared with the
target and with every Gaia star close enough and bright enough to make the dip (an eclipse of at most 100%).

This replaces, for the report's claims, the earlier source test (vet_deep._source_test: the difference image
fitted as a shifted copy of the whole out-of-transit image), which a calibration on published quadruple stars
(eb_pixels.py --kostov) and ZTF light curves (ztf_check.py) showed to be unreliable in crowded fields.

Pixels: SPOC 2-min target pixel files where they exist, else 15x15 full-frame cutouts (astrocut).
Input: the events of results/candidates/deep_vet.json, plus any given on the command line.

  python centroid_vet.py                      (every star in deep_vet.json)
  python centroid_vet.py TIC:dur_h:depth:sector@tc,sector@tc,...
  python centroid_vet.py --recombine          (re-derive verdicts from the stored per-transit centroids)
Output: results/candidates/centroid_vet.json
"""
import io
import json
import sys

import env  # noqa: F401
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

import pixels
import s3io
import ticdb
from centroids import _centroid, centroid_verdict, combine_centroids, gaia_field, neighbour_offsets
from lcaccess import all_sectors
from env import RESULTS

C = RESULTS / "candidates"
OUT = C / "centroid_vet.json"


def cube(tic, sector, ra, dec, secs):
    if sector in secs.get("spoc", []):
        tm = s3io.spoc_templates()[str(sector)]
        key = f"tess/public/tid/s{sector:04d}/{s3io.tic_dir(tic)}/{tm[0]}-s{sector:04d}-{tic:016d}-{tm[1]}-s_tp.fits"
        raw = s3io.fetch(key)
        if raw is not None:
            h = fits.open(io.BytesIO(raw))
            kind = "2-min pixels"
        else:
            h = None
    else:
        h = None
    if h is None:
        from tess_stars2px import tess_stars2px_function_entry as t2p
        r = t2p([0], [ra], [dec])
        cc = [(c, d) for s, c, d in zip(r[3], r[4], r[5]) if int(s) == int(sector)]
        if not cc:
            return None
        path = pixels.cutout(ra, dec, sector, cc[0][0], cc[0][1], size=15)
        h = fits.open(path)
        kind = "full-frame cutout"
    d = h[1].data
    t, f, q = d["TIME"], d["FLUX"], d["QUALITY"]
    g = np.isfinite(t) & (q == 0) & np.isfinite(f).all(axis=(1, 2))
    return np.asarray(t[g], float), np.asarray(f[g], float), WCS(h[2].header), kind


def vet(tic, events, dur_h, depth):
    st = ticdb.query([tic], ["ra", "dec", "Tmag"]).iloc[0]
    ra, dec, tmag = float(st.ra), float(st.dec), float(st.Tmag)
    secs = all_sectors(tic)
    gf = gaia_field(ra, dec, tmag)
    D = dur_h / 24
    cents, per = [], []
    by_sector = {}
    for s, tc in events:
        by_sector.setdefault(int(s), []).append(float(tc))
    for s, tcs in sorted(by_sector.items()):
        try:
            cb = cube(tic, s, ra, dec, secs)
        except Exception as ex:
            per.append({"sector": s, "error": repr(ex)[:120]})
            continue
        if cb is None:
            per.append({"sector": s, "error": "no pixels"})
            continue
        t, F, w, kind = cb
        for tc in tcs:
            c = _centroid(t, F, w, tc, D, ra, dec, min_snr=4.0, gf=gf)
            row = {"sector": s, "tc": tc, "kind": kind}
            if c:
                row.update(c)
                if "dra" in c:
                    cents.append(dict(c, sector=s))
            per.append(row)
    cand = neighbour_offsets(gf, depth, 150.0)
    cen = combine_centroids(cents)
    return {"tic": tic, "dur_h": dur_h, "depth": depth, "per_event": per, "centroid": cen, "n_candidates": len(cand),
            **centroid_verdict(cen, cand)}


def recombine():
    """Re-derive every star's combined centroid and verdict from the stored per-transit centroids (no pixel
    data needed), e.g. after a change to centroids.combine_centroids."""
    res = json.load(open(OUT))
    for tic, r in res.items():
        if "per_event" not in r:
            continue
        st = ticdb.query([int(tic)], ["ra", "dec", "Tmag"]).iloc[0]
        gf = gaia_field(float(st.ra), float(st.dec), float(st.Tmag))
        cen = combine_centroids([e for e in r["per_event"] if "dra" in e])
        r.update(centroid=cen, n_candidates=len(neighbour_offsets(gf, r["depth"], 150.0)),
                 **centroid_verdict(cen, neighbour_offsets(gf, r["depth"], 150.0)))
    json.dump(res, open(OUT, "w"), indent=1, default=float)
    return res


def from_deep_vet():
    d = json.load(open(C / "deep_vet.json"))
    out = []
    for tic, o in d.items():
        ev = [(e["sector"], e["tc"]) for e in o.get("events", []) if e.get("sector")]
        dur = None
        for e in o.get("events", []):
            dur = dur or e.get("dur_h")
        dep = o.get("depth_used") or np.median([e["depth"] for e in o["events"] if e.get("depth")] or [1e-3])
        out.append((int(tic), ev, dur, float(dep)))
    return out


def main(args):
    todo = []
    if args:
        for a in args:
            tic, dur, dep, evs = a.split(":")
            todo.append((int(tic), [(int(x.split("@")[0]), float(x.split("@")[1])) for x in evs.split(",")], float(dur), float(dep)))
    else:
        # durations: from the Track B table when the deep-vet record has none
        import pandas as pd
        tb = pd.read_csv(C / "trackb_final.csv")
        durs = dict(zip(tb.tic.astype(int), tb.t14_h))
        for tic, ev, dur, dep in from_deep_vet():
            todo.append((tic, ev, dur or durs.get(tic) or 4.0, dep))
    res = json.load(open(OUT)) if OUT.exists() else {}
    for tic, ev, dur, dep in todo:
        if str(tic) in res and not args:
            continue
        try:
            r = vet(tic, ev, dur, dep)
        except Exception as ex:
            r = {"tic": tic, "verdict": "error", "error": repr(ex)[:200]}
        res[str(tic)] = r
        c = r.get("centroid") or {}
        print(tic, r.get("verdict"), round(r.get("offset_arcsec", -1), 1), "±", round(c.get("err", -1), 1),
              "n", c.get("n_eclipses"), r.get("closest_neighbours", [])[:1], flush=True)
        json.dump(res, open(OUT, "w"), indent=1, default=float)


if __name__ == "__main__":
    if sys.argv[1:] == ["--recombine"]:
        recombine()
    else:
        main(sys.argv[1:])
