"""Which star makes each eclipse set? Pixel test for the Track C doubly eclipsing candidates.

For each star, up to three sectors where both eclipse signals were found:
  - the eclipse analysis is re-run (ebsearch.analyse_eb) to get both ephemerides in that sector;
  - for every eclipse of each set, the difference image (out-of-eclipse minus in-eclipse) is fitted with the
    out-of-eclipse image shifted to the target and to each Gaia star within 120" that is bright enough to make
    that eclipse depth (vet_deep._source_test); delta-chi2 is summed over the eclipses of the set.
Both sets on the target: the two binaries are in one unresolved system (2+2 quadruple or triple candidate).
A set that favours a neighbour: a blend of two unrelated binaries.
Second, independent test (the standard one): for every eclipse, a Gaussian PSF with free position is fitted to the
difference image; the offsets from the target are combined over all eclipses of a set (3-sigma clipped, error from
their scatter, with a 2" floor) and compared with the target and with every Gaia star able to make the depth.
Output: results/candidates/eb_pixels.json

  python eb_pixels.py [tic ...]      (default: the 30 stars of the report's Track C table)
"""
import json
import math
import sys

import env  # noqa: F401
import numpy as np
import pandas as pd

import pixels
import s3io
import ticdb
from ebsearch import analyse_eb
from env import RESULTS
from vet_deep import _source_test
from centroids import _centroid, centroid_verdict, combine_centroids, gaia_field, neighbour_offsets, pixel_cube

C = RESULTS / "candidates"


def test_star(row, max_sectors=3):
    tic = int(row["tic"])
    st = ticdb.query([tic], ["ra", "dec", "Tmag"]).iloc[0]
    ra, dec, tmag = float(st.ra), float(st.dec), float(st.Tmag)
    g = pixels.gaia_cone(ra, dec, 120.0)
    gf = gaia_field(ra, dec, tmag)
    out = {"tic": tic, "P1_cat": float(row["P1"]), "P2_cat": float(row["P2"]), "sectors": []}
    secs = [int(s) for s in str(row["sectors"]).split(";") if s]
    have = {p: set(v) for p, v in __import__("dossier").all_sectors(tic).items()}
    done = 0
    cents, depths = {1: [], 2: []}, {1: [], 2: []}
    for s in secs[::-1]:  # latest sectors first: best data
        if done >= max_sectors:
            break
        product = "spoc" if s in have.get("spoc", ()) else "qlp" if s in have.get("qlp", ()) else "tspoc" if s in have.get("tspoc", ()) else None
        if product is None:
            continue
        try:
            lc = s3io.load(product if product != "spoc" else "spoc", tic, s)
            if lc is None and product == "spoc":
                product, lc = "qlp", s3io.load("qlp", tic, s)
            if lc is None:
                continue
            a = analyse_eb(lc)
            if "P2" not in a:
                continue
            # the ephemerides must match the catalogued periods (allowing a factor of 2)
            def match(p, ref):
                return any(abs(p * k / ref - 1) < 0.01 for k in (0.5, 1, 2))
            if not (match(a["P1"], row["P1"]) and match(a["P2"], row["P2"])):
                out["sectors"].append({"sector": s, "skipped": f"periods {a['P1']:.4f}/{a['P2']:.4f} differ"})
                continue
            t_cube, f_cube, w, kind = pixel_cube(tic, s, product, lc, ra, dec)
            sec_res = {"sector": s, "kind": kind, "sets": {}}
            for j in (1, 2):
                P, t0, D, dep = a[f"P{j}"], a[f"t0_{j}"], a[f"dur{j}_h"] / 24, a[f"depth{j}"]
                depths[j].append(dep)
                lim = tmag - 2.5 * math.log10(max(dep, 1e-4))  # neighbours able to make this depth at <= 100% eclipse
                srcs = [{"name": f"G{r.phot_g_mean_mag:.1f} at {r.sep_arcsec:.0f}\"", "ra": float(r.ra), "dec": float(r.dec),
                         "sep": float(r.sep_arcsec)} for _, r in g.iterrows() if r.sep_arcsec > 3 and r.tmag_est < lim]
                k0 = int(np.ceil((t_cube.min() - t0) / P)); k1 = int(np.floor((t_cube.max() - t0) / P))
                comb, n_ev = {}, 0
                for k in range(k0, k1 + 1):
                    tc = t0 + k * P
                    cen = _centroid(t_cube, f_cube, w, tc, D, ra, dec, gf=gf)
                    if cen and "dra" in cen:
                        cents[j].append(dict(cen, sector=s, tc=float(tc)))
                    r = _source_test(t_cube, f_cube, w, tc, D, ra, dec, srcs)
                    if not r or r["snr_diff"] < 5:
                        continue
                    n_ev += 1
                    for src in r["sources"]:
                        c = comb.setdefault(src["name"], {"dchi2": 0.0, "sep": src["sep"]})
                        c["dchi2"] += src["dchi2"]
                worst = min(comb.items(), key=lambda kv: kv[1]["dchi2"]) if comb else None
                verdict = ("no usable eclipses" if n_ev == 0 else
                           "on target" if not worst or worst[1]["dchi2"] > 9 else
                           f"favours {worst[0]}" if worst[1]["dchi2"] < -9 else f"ambiguous ({worst[0]})")
                sec_res["sets"][j] = {"P": P, "t0": t0, "dur_h": D * 24, "depth": dep, "n_eclipses": n_ev,
                                      "n_sources": len(srcs), "best_alternative": worst[0] if worst else None,
                                      "best_dchi2": worst[1]["dchi2"] if worst else None, "verdict": verdict,
                                      "dchi2_by_source": {k: v["dchi2"] for k, v in comb.items()}}
            out["sectors"].append(sec_res)
            done += 1
        except Exception as ex:
            out["sectors"].append({"sector": s, "error": repr(ex)[:200]})
    # combine over sectors: delta-chi2 per neighbour summed over every sector where it was tested
    # (sources are matched by position, rounded to 1", since names carry rounded magnitudes)
    v = [x["sets"] for x in out["sectors"] if "sets" in x]
    summary = {}
    for j in (1, 2):
        tot = {}
        for sset in v:
            for name, dc in sset.get(j, {}).get("dchi2_by_source", {}).items():
                tot[name] = tot.get(name, 0.0) + dc
        n_sec = sum(1 for sset in v if j in sset and sset[j]["n_eclipses"] > 0)
        worst = min(tot.items(), key=lambda kv: kv[1]) if tot else None
        summary[j] = {"n_sectors": n_sec, "best_alternative": worst[0] if worst else None,
                      "combined_dchi2": worst[1] if worst else None,
                      "verdict": ("untested" if n_sec == 0 else "on target" if not worst or worst[1] > 9 else
                                  f"favours {worst[0]}" if worst[1] < -9 else "ambiguous"),
                      "per_sector": [sset[j]["verdict"] for sset in v if j in sset]}
    out["sets"] = summary
    # independent second test: difference-image centroids, all eclipses of a set over all sectors
    cen_res = {}
    for j in (1, 2):
        dep = float(np.median(depths[j])) if depths[j] else None
        cand = neighbour_offsets(gf, dep, 120.0) if dep else []
        cen = combine_centroids(cents[j])
        cen_res[j] = dict(centroid=cen, **centroid_verdict(cen, cand), n_candidates=len(cand),
                          per_eclipse=[{k: (round(v, 3) if isinstance(v, float) else v) for k, v in x.items()} for x in cents[j]])
    out["centroid"] = cen_res
    c1, c2 = cen_res[1]["verdict"], cen_res[2]["verdict"]
    out["centroid_verdict"] = ("both sets on target" if c1 == c2 == "on target" else
                               "a set on a neighbour" if c1.startswith("on G") or c2.startswith("on G") else
                               "unresolved or untested")
    s1, s2 = summary[1]["verdict"], summary[2]["verdict"]
    if "untested" in (s1, s2):
        out["verdict"] = "untested"
    elif s1.startswith("favours") or s2.startswith("favours"):
        out["verdict"] = "blend: an eclipse set comes from a neighbour"
    elif s1 == s2 == "on target":
        split = any(x.startswith("favours") for j in (1, 2) for x in summary[j]["per_sector"])
        out["verdict"] = "both eclipse sets on the target" + (" (one sector dissents)" if split else "")
    else:
        out["verdict"] = "inconclusive"
    return out


def recombine(path):
    """Re-derive the centroid verdicts of a results file from its stored per-eclipse centroids (no pixel data
    needed), e.g. after a change to combine_centroids."""
    res = json.load(open(path))
    for tic, o in res.items():
        if "centroid" not in o:
            continue
        st = ticdb.query([int(tic)], ["ra", "dec", "Tmag"]).iloc[0]
        ra, dec, tmag = float(st.ra), float(st.dec), float(st.Tmag)
        gf = gaia_field(ra, dec, tmag)
        for j, v in o["centroid"].items():
            deps = [s["sets"][j]["depth"] for s in o["sectors"] if "sets" in s and j in s["sets"]]
            cand = neighbour_offsets(gf, float(np.median(deps)), 120.0) if deps else []
            cen = combine_centroids(v["per_eclipse"])
            v.update(centroid=cen, **centroid_verdict(cen, cand), n_candidates=len(cand))
        c1, c2 = o["centroid"]["1"]["verdict"], o["centroid"]["2"]["verdict"]
        o["centroid_verdict"] = ("both sets on target" if c1 == c2 == "on target" else
                                 "a set on a neighbour" if c1.startswith("on G") or c2.startswith("on G") else
                                 "unresolved or untested")
    json.dump(res, open(path, "w"), indent=1, default=float)
    return res


def kostov_rows(max_period=10.0):
    """Calibration sample: Kostov et al. (2022) quadruple candidates, whose two eclipse sets the authors traced to
    one star with their own pixel analysis (VizieR J/ApJS/259/66, saved by hand to the work catalogue folder)."""
    from dossier import all_sectors
    from env import CATALOGS
    k = pd.read_csv(CATALOGS / "kostov2022_quadruples.csv")
    rows = []
    for tic, g in k.groupby("TIC"):
        if len(g) != 2:
            continue
        p1, p2 = sorted(g.Per.astype(float))
        if not (p2 < max_period):
            continue
        secs = sorted({int(s) for v in all_sectors(int(tic)).values() for s in v})
        rows.append({"tic": int(tic), "P1": p1, "P2": p2, "sectors": ";".join(map(str, secs))})
    return pd.DataFrame(rows)


def main(tics=None, calib=False):
    if calib:
        x, out = kostov_rows(), C / "eb_pixels_kostov.json"
    else:
        tc = pd.read_csv(C / "trackc_candidates.csv")
        x, out = tc[tc.P2_recurs >= 2].sort_values(["P2_recurs", "snr2_max"], ascending=False).head(30), C / "eb_pixels.json"
    if tics:
        x = x[x.tic.isin(tics)]
    res = json.load(open(out)) if out.exists() else {}
    for _, r in x.iterrows():
        if str(int(r.tic)) in res and tics is None:
            continue
        try:
            o = test_star(r)
        except Exception as ex:
            o = {"tic": int(r.tic), "verdict": "error", "error": repr(ex)[:200], "sectors": []}
        res[str(int(r.tic))] = o
        print(int(r.tic), o["verdict"], "| centroid:", o.get("centroid_verdict"),
              {j: (v["verdict"], round(v.get("offset_arcsec", -1), 1), round(v.get("offset_sigma", -1), 1)) for j, v in o.get("centroid", {}).items()},
              [(s.get("sector"), {j: v["verdict"] for j, v in s.get("sets", {}).items()}) for s in o["sectors"]], flush=True)
        json.dump(res, open(out, "w"), indent=1, default=float)


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if a != "--kostov"]
    main([int(a) for a in args] or None, calib="--kostov" in sys.argv)
