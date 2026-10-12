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
import io
import json
import math
import sys

import env  # noqa: F401
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS

import pixels
import s3io
import ticdb
from ebsearch import analyse_eb
from env import RESULTS
from vet_deep import _source_test

C = RESULTS / "candidates"


def pixel_cube(tic, sector, product, lc, ra, dec):
    if product == "spoc":
        tm = s3io.spoc_templates()[str(sector)]
        key = f"tess/public/tid/s{sector:04d}/{s3io.tic_dir(tic)}/{tm[0]}-s{sector:04d}-{tic:016d}-{tm[1]}-s_tp.fits"
        h = fits.open(io.BytesIO(s3io.fetch(key)))
        kind = "2-min pixels"
    else:
        path = pixels.cutout(ra, dec, sector, lc["meta"]["camera"], lc["meta"]["ccd"], size=15)
        h = fits.open(path)
        kind = "full-frame cutout"
    d = h[1].data
    t, f, q = d["TIME"], d["FLUX"], d["QUALITY"]
    g = np.isfinite(t) & (q == 0) & np.isfinite(f).all(axis=(1, 2))
    return np.asarray(t[g], float), np.asarray(f[g], float), WCS(h[2].header), kind


def _diff_image(t, F, tc, D):
    """Out-of-eclipse (linear fit around the eclipse) minus in-eclipse image, with per-pixel errors."""
    inn = np.abs(t - tc) < 0.35 * D
    oot = (np.abs(t - tc) > 0.7 * D) & (np.abs(t - tc) < 0.7 * D + max(0.4, 2 * D))
    if inn.sum() < 2 or oot.sum() < 6:
        return None
    A = np.vstack([np.ones(oot.sum()), t[oot] - tc]).T
    ny, nx = F.shape[1:]
    coef, *_ = np.linalg.lstsq(A, F[oot].reshape(oot.sum(), -1), rcond=None)
    base = coef[0].reshape(ny, nx)
    diff = base - F[inn].mean(axis=0)
    err = ((F[oot].reshape(oot.sum(), -1) - A @ coef).std(axis=0) / math.sqrt(inn.sum())).reshape(ny, nx)
    return diff, err, base


GAIA_COLS = ("source_id", "ra", "dec", "phot_g_mean_mag", "phot_bp_mean_mag", "phot_rp_mean_mag", "pmra", "pmdec")


def gaia_field(ra, dec, tmag, radius=250.0):
    """Gaia DR3 stars around a target, with proper motions and estimated T magnitudes; the target (the star
    within 6" whose estimated T is closest to the TIC value) is flagged."""
    g = pixels.gaia_cone(ra, dec, radius, cols=GAIA_COLS)
    g["is_target"] = False
    near = g[g.sep_arcsec < 6]
    if len(near):
        g.loc[(near.tmag_est - tmag).abs().idxmin(), "is_target"] = True
    return g


def _calibrate(base, w, gf, t_mid):
    """Where the stars really are on this image: the out-of-transit image fitted with Gaussians at the Gaia
    positions (moved to the observing date with their proper motions), fluxes in Gaia's proportions, one common
    shift, width, scale and background. Returns the target's pixel position and the shift."""
    from scipy.optimize import least_squares
    tg = gf[gf.is_target]
    if tg.empty:
        return None
    dt = 2000.0 + (t_mid + 2457000.0 - 2451545.0) / 365.25 - 2016.0
    cosd = math.cos(math.radians(float(tg.dec.iloc[0])))
    ra_s = gf.ra + dt * gf.pmra.fillna(0) / 3.6e6 / cosd
    de_s = gf.dec + dt * gf.pmdec.fillna(0) / 3.6e6
    xs, ys = (np.asarray(v, float) for v in w.all_world2pix(ra_s.values, de_s.values, 0))
    ny, nx = base.shape
    t_tm = float(tg.tmag_est.iloc[0])
    use = (xs > -3) & (xs < nx + 2) & (ys > -3) & (ys < ny + 2) & (gf.tmag_est.values < t_tm + 6)
    i_t = int(np.where(gf.is_target.values)[0][0])
    if not use[i_t]:
        return None
    xs, ys = xs[use], ys[use]
    fl = 10 ** (-0.4 * (gf.tmag_est.values[use] - t_tm))
    yy, xx = np.mgrid[0:ny, 0:nx]
    good = np.isfinite(base)
    err = np.sqrt(np.abs(base)) + 1.0

    def model(p):
        a, dx, dy, s, c = p
        img = np.zeros_like(base)
        for x0, y0, f0 in zip(xs, ys, fl):
            img += f0 * np.exp(-((xx - x0 - dx) ** 2 + (yy - y0 - dy) ** 2) / (2 * s * s))
        return a * img + c

    a0 = float(np.nanmax(base))
    try:
        r = least_squares(lambda p: ((model(p) - base) / err)[good], [a0, 0.0, 0.0, 0.9, float(np.nanmedian(base))],
                          x_scale="jac", bounds=([0, -2, -2, 0.4, -np.inf], [np.inf, 2, 2, 3, np.inf]))
    except Exception:
        return None
    xt0, yt0 = (float(v) for v in w.all_world2pix(ra_s.values[i_t], de_s.values[i_t], 0))
    return {"x": xt0 + r.x[1], "y": yt0 + r.x[2], "shift_px": [float(r.x[1]), float(r.x[2])], "psf_sigma_px": float(r.x[3]),
            "ra": float(ra_s.values[i_t]), "dec": float(de_s.values[i_t])}


def _centroid(t, F, w, tc, D, ra0, dec0, min_snr=5.0, gf=None):
    """Where on the sky the light that disappears in one eclipse comes from: a Gaussian PSF (free position, width,
    amplitude and constant) fitted to the difference image. Returns the offset from the target in arcsec.
    With gf (gaia_field), the offset is measured from where the target really is on this image (_calibrate), which
    removes errors of the image's coordinate solution and the target's proper motion; otherwise from the
    catalogue position through the image's coordinate solution."""
    from astropy.wcs.utils import proj_plane_pixel_scales
    from scipy.optimize import least_squares
    r = _diff_image(t, F, tc, D)
    if r is None:
        return None
    diff, err, base = r
    m = np.isfinite(diff) & np.isfinite(err) & (err > 0)
    snr = float(np.sqrt(max(np.sum((diff[m] / err[m]) ** 2) - m.sum(), 0)))
    if snr < min_snr:
        return {"snr": snr}
    ny, nx = diff.shape
    yy, xx = np.mgrid[0:ny, 0:nx]
    xt, yt = (float(v) for v in w.all_world2pix(ra0, dec0, 0))

    def resid(p):
        a, x, y, s, c = p
        return ((a * np.exp(-((xx - x) ** 2 + (yy - y) ** 2) / (2 * s * s)) + c - diff) / err)[m]
    sm = np.where(m, diff, -np.inf)
    iy, ix = np.unravel_index(int(np.argmax(sm)), sm.shape)
    best = None
    for x0, y0 in ((xt, yt), (float(ix), float(iy))):
        x0, y0 = min(max(x0, 0), nx - 1), min(max(y0, 0), ny - 1)
        a0 = max(float(diff[int(round(y0)), int(round(x0))]), float(np.nanmax(sm)) * 0.5, 1e-6)
        try:
            res = least_squares(resid, [a0, x0, y0, 0.9, 0.0], x_scale="jac",
                                bounds=([0, -0.5, -0.5, 0.4, -np.inf], [np.inf, nx - 0.5, ny - 0.5, 2.5, np.inf]))
        except Exception:
            continue
        if best is None or res.cost < best.cost:
            best = res
    if best is None:
        return {"snr": snr}
    dof = max(int(m.sum()) - 5, 1)
    chi2 = 2 * best.cost
    try:
        cov = np.linalg.inv(best.jac.T @ best.jac) * max(chi2 / dof, 1.0)
        sx, sy = math.sqrt(cov[1, 1]), math.sqrt(cov[2, 2])
    except np.linalg.LinAlgError:
        return {"snr": snr}
    x, y = best.x[1], best.x[2]
    ra, dec = (float(v) for v in w.all_pix2world(x, y, 0))
    scale = float(np.mean(proj_plane_pixel_scales(w))) * 3600
    out = {"snr": snr, "dra": (ra - ra0) * math.cos(math.radians(dec0)) * 3600, "ddec": (dec - dec0) * 3600,
           "err": math.hypot(sx, sy) / math.sqrt(2) * scale, "psf_sigma_px": float(best.x[3]),
           "redchi2": chi2 / dof, "edge": bool(x < 1 or y < 1 or x > nx - 2 or y > ny - 2), "calibrated": False}
    if gf is not None:
        cal = _calibrate(base, w, gf, float(tc))
        if cal is not None:
            rt, dt_ = (float(v) for v in w.all_pix2world(cal["x"], cal["y"], 0))
            out.update(dra=(ra - rt) * math.cos(math.radians(dt_)) * 3600, ddec=(dec - dt_) * 3600, calibrated=True,
                       cal_shift_px=cal["shift_px"], cal_psf_sigma_px=cal["psf_sigma_px"])
    return out


def neighbour_offsets(gf, depth, max_sep=150.0):
    """Gaia stars able to make a dip of this depth (an eclipse of at most 100%), as offsets from the target's
    Gaia position (same epoch, so offsets are not blurred by the target's proper motion)."""
    tg = gf[gf.is_target]
    if tg.empty:
        return []
    t = tg.iloc[0]
    lim = float(t.tmag_est) - 2.5 * math.log10(max(depth, 1e-5))
    out = []
    for _, r in gf.iterrows():
        if r.is_target or r.tmag_est >= lim:
            continue
        dra = (float(r.ra) - float(t.ra)) * math.cos(math.radians(float(t.dec))) * 3600
        ddec = (float(r.dec) - float(t.dec)) * 3600
        sep = math.hypot(dra, ddec)
        if 1.0 < sep < max_sep:
            out.append({"name": f"G{r.phot_g_mean_mag:.1f} at {sep:.0f}\"", "dra": dra, "ddec": ddec})
    return out


def combine_centroids(cents, floor=2.0, clip=3.5):
    """Robust weighted mean of per-eclipse centroid offsets. Each offset gets its fitted error plus a 2"
    systematic floor (about 0.1 pixel); offsets further than 3.5 sigma (2-D) from the mean are dropped in turn;
    the error of the mean is scaled up by sqrt(reduced chi-square) when the offsets scatter more than their
    errors say. Eclipses whose difference image a Gaussian fits very badly are dropped first, with an allowance
    that grows with signal-to-noise (a deep eclipse of a bright star never fits a Gaussian exactly)."""
    c = [x for x in cents if x and "dra" in x and not x["edge"] and x["redchi2"] < max(50.0, 0.002 * x["snr"] ** 2)]
    if not c:
        return None
    X = np.array([[x["dra"], x["ddec"], math.hypot(x["err"], floor)] for x in c])
    keep = np.ones(len(X), bool)
    while True:
        w = 1 / X[keep, 2] ** 2
        mx, my = np.sum(w * X[keep, 0]) / w.sum(), np.sum(w * X[keep, 1]) / w.sum()
        d = np.hypot(X[:, 0] - mx, X[:, 1] - my) / X[:, 2]
        worst = int(np.argmax(np.where(keep, d, -1)))
        if keep.sum() <= 2 or d[worst] < clip:
            break
        keep[worst] = False
    n = int(keep.sum())
    err = 1 / math.sqrt(w.sum())
    if n >= 2:
        chi2 = float(np.sum(d[keep] ** 2))
        red = chi2 / (2 * (n - 1))
        err *= math.sqrt(max(red, 1.0))
    else:
        red = None
    secs = sorted({x.get("sector") for x in c})
    return {"dra": float(mx), "ddec": float(my), "err": float(err), "n_eclipses": n, "n_rejected": int(len(X) - n),
            "rejected_sectors": sorted({c[i].get("sector") for i in np.where(~keep)[0]}), "reduced_chi2": red,
            "n_sectors": len(secs)}


def centroid_verdict(cen, srcs):
    """srcs: Gaia neighbours able to make the depth, with offsets (dra, ddec) from the target in arcsec."""
    if cen is None:
        return {"verdict": "no usable centroids"}
    d_t = math.hypot(cen["dra"], cen["ddec"]) / cen["err"]
    near = sorted(((math.hypot(cen["dra"] - s["dra"], cen["ddec"] - s["ddec"]) / cen["err"], s["name"]) for s in srcs))
    ok_n = [nm for d, nm in near if d < 3]
    v = ("on target" if d_t < 3 and not ok_n else
         f"target or {ok_n[0]} (unresolved)" if d_t < 3 else
         f"on {ok_n[0]}" if ok_n else "offset, no catalogued source fits")
    return {"verdict": v, "offset_arcsec": float(math.hypot(cen["dra"], cen["ddec"])), "offset_sigma": float(d_t),
            "closest_neighbours": [{"name": nm, "sigma": float(d)} for d, nm in near[:3]]}


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
