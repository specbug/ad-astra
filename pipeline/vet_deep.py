"""Deep vetting of individual transit candidates (results/candidates/deep_vet.json).

For each star and its list of transit times:
  star      TIC parameters; Gaia DR3 RUWE, non-single-star flag and RV scatter of the target
  known     TOIs by TIC or within 10', SPOC TCEs from the DV reports
  events    depth of every listed event in its sector's best light curve (SPOC 2-min > TESS-SPOC > QLP)
  period    for >= 2 events: the longest common period (finalize_b.ephemeris) and each alias P/n tested
            against all sectors (predicted transits inside data must show the dip); for 1 event: the
            shortest period the other sectors allow (duo.single_period_constraints)
  pixels    per event, a source test on the pixels (2-min target pixel file, else a full-frame cutout):
            the difference image (out of transit minus in transit) is fitted as a scaled copy of the
            out-of-transit image shifted to each candidate source (the target, and every Gaia star within
            150" bright enough to cause the dip with an eclipse of at most 50%); dchi2 > 0 favours the target.

  python vet_deep.py TIC:t1;t2;... [TIC:...]      or      python vet_deep.py --shortlist [Tmax]
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
from scipy.ndimage import shift as ndshift

import known
import pixels
import s3io
import ticdb
from dossier import all_sectors
from duo import flatten_sectors, single_period_constraints
from env import RESULTS
from finalize_b import ephemeris
from lcprep import bin_lc, detrend_fast, robust_std

OUT = RESULTS / "candidates" / "deep_vet.json"


def sector_of(t, lcs):
    for lc in lcs:
        if lc["time"].min() - 0.5 <= t <= lc["time"].max() + 0.5:
            return lc
    return None


def best_lcs_all(tic):
    secs = all_sectors(tic)
    got, out = set(), []
    for p in ("spoc", "tspoc", "qlp"):
        for s in secs.get(p, []):
            if s in got:
                continue
            try:
                lc = s3io.load(p, tic, s)
            except Exception:
                lc = None
            if lc is not None and lc["time"].size > 100:
                out.append(lc)
                got.add(s)
    return sorted(out, key=lambda l: l["meta"]["sector"]), secs


def flat10(lcs, window=1.5):
    T, Y, S = [], [], []
    for lc in lcs:
        t, f = lc["time"], lc["flux"]
        ok = np.isfinite(t) & np.isfinite(f)
        tb, fb, _, _ = bin_lc(t[ok], f[ok], 10.0)
        y, _ = detrend_fast(tb, fb, window, knot_step=window / 36)
        g = np.isfinite(y)
        T.append(tb[g]); Y.append(y[g]); S.append(np.full(g.sum(), lc["meta"]["sector"]))
    o = np.argsort(np.concatenate(T))
    return np.concatenate(T)[o], np.concatenate(Y)[o], np.concatenate(S)[o]


def box_depth(t, y, tc, D):
    inn = np.abs(t - tc) < 0.35 * D
    if inn.sum() < 3:
        return None
    oot = (np.abs(t - tc) > 0.75 * D) & (np.abs(t - tc) < max(0.8, 3 * D))
    return float(np.nanmedian(y[oot]) - np.nanmean(y[inn])), float(robust_std(y[oot]) / math.sqrt(inn.sum())), int(inn.sum())


def alias_test(t, y, sec, times, depth, D, nmax=12):
    """Aliases P/n of the longest common period of `times` (>= 2 events)."""
    eph = ephemeris(times, D * 24, pmin=0.5)
    if eph is None:
        return None
    P1 = eph[0]
    rows = []
    times = np.asarray(times)
    for n in range(1, nmax + 1):
        P = P1 / n
        ep = np.round((times - times[0]) / P)
        A = np.vstack([np.ones_like(ep), ep]).T
        (t0, Pf), *_ = np.linalg.lstsq(A, times, rcond=None)
        rms = float(np.std(times - (t0 + Pf * ep)) * 1440)
        k = np.arange(np.floor((t.min() - t0) / Pf) - 1, np.ceil((t.max() - t0) / Pf) + 2)
        worst = (0.0, None, None, None)
        ncov = 0
        for p in t0 + k * Pf:
            if np.min(np.abs(p - times)) < 0.5 * D:
                continue
            r = box_depth(t, y, p, D)
            if r is None:
                continue
            ncov += 1
            z = (depth - r[0]) / r[1]
            if z > worst[0]:
                worst = (float(z), float(p), r[0], int(sec[np.argmin(np.abs(t - p))]))
        rows.append({"n": n, "P": float(Pf), "rms_min": rms, "n_covered": ncov, "max_reject_sigma": worst[0],
                     "at": worst[1], "depth_there": worst[2], "sector": worst[3],
                     "status": "excluded" if worst[0] > 5 else "allowed"})
    return {"P_longest": P1, "rms_min": eph[1], "aliases": rows}


def _source_test(cube_t, cube_f, wcs, tc, D, ra0, dec0, srcs):
    inn = np.abs(cube_t - tc) < 0.35 * D
    oot = (np.abs(cube_t - tc) > 0.7 * D) & (np.abs(cube_t - tc) < 0.7 * D + max(0.4, 2 * D))
    if inn.sum() < 2 or oot.sum() < 6:
        return None
    A = np.vstack([np.ones(oot.sum()), cube_t[oot] - tc]).T
    ny, nx = cube_f.shape[1:]
    coef, *_ = np.linalg.lstsq(A, cube_f[oot].reshape(oot.sum(), -1), rcond=None)
    base = coef[0].reshape(ny, nx)
    diff = base - cube_f[inn].mean(axis=0)
    err = ((cube_f[oot].reshape(oot.sum(), -1) - A @ coef).std(axis=0) / math.sqrt(inn.sum())).reshape(ny, nx)
    m = np.isfinite(diff) & np.isfinite(base) & (err > 0)
    xt, yt = [float(v) for v in wcs.all_world2pix(ra0, dec0, 0)]

    def chi2(dx, dy):
        tm = ndshift(np.where(np.isfinite(base), base, 0), (dy, dx), order=1, mode="constant", cval=0.0)
        a = np.sum(tm[m] * diff[m] / err[m] ** 2) / max(np.sum(tm[m] ** 2 / err[m] ** 2), 1e-30)
        return float(np.sum(((diff[m] - a * tm[m]) / err[m]) ** 2))
    c0 = chi2(0, 0)
    out = {"chi2_target": c0, "snr_diff": float(np.sqrt(max(np.sum((diff[m] / err[m]) ** 2) - m.sum(), 0))),
           "sources": []}
    for s in srcs:
        x, y = [float(v) for v in wcs.all_world2pix(s["ra"], s["dec"], 0)]
        dx, dy = x - xt, y - yt
        if abs(dx) > nx - 1 or abs(dy) > ny - 1:
            continue
        out["sources"].append({**s, "dchi2": chi2(dx, dy) - c0, "dx": dx, "dy": dy})
    return out


def pixel_tests(tic, events, D, ra0, dec0, depth, secs):
    g = pixels.gaia_cone(ra0, dec0, 150.0)
    tgt_t = ticdb.query([tic], ["Tmag"]).Tmag.iloc[0]
    # sources able to cause the dip with an eclipse <= 50%
    lim = tgt_t - 2.5 * math.log10(2 * depth)
    srcs = [{"name": f"G{r.phot_g_mean_mag:.1f} at {r.sep_arcsec:.0f}\"", "ra": float(r.ra), "dec": float(r.dec),
             "tmag": float(r.tmag_est), "sep": float(r.sep_arcsec), "needed_eclipse": float(depth * 10 ** (0.4 * (r.tmag_est - tgt_t)))}
            for _, r in g.iterrows() if r.sep_arcsec > 3 and r.tmag_est < lim]
    res = []
    for sector, tc in events:
        r = None
        try:
            if sector in secs.get("spoc", []):
                tm = s3io.spoc_templates()[str(sector)]
                key = f"tess/public/tid/s{sector:04d}/{s3io.tic_dir(tic)}/{tm[0]}-s{sector:04d}-{tic:016d}-{tm[1]}-s_tp.fits"
                h = fits.open(io.BytesIO(s3io.fetch(key)))
                d = h[1].data
                w = WCS(h[2].header)
                t, f, q = d["TIME"], d["FLUX"], d["QUALITY"]
                kind = "2-min pixels"
            else:
                lc = s3io.load("qlp", tic, sector) or s3io.load("tspoc", tic, sector)
                cam, ccd = lc["meta"]["camera"], lc["meta"]["ccd"]
                path = pixels.cutout(ra0, dec0, sector, cam, ccd, size=15)
                h = fits.open(path)
                d = h[1].data
                w = WCS(h[2].header)
                t, f, q = d["TIME"], d["FLUX"], d["QUALITY"]
                kind = "full-frame cutout"
            g_ = np.isfinite(t) & (q == 0) & np.isfinite(f).all(axis=(1, 2))
            r = _source_test(np.asarray(t[g_], float), np.asarray(f[g_], float), w, tc, D, ra0, dec0, srcs)
            if r:
                r.update({"sector": int(sector), "tc": float(tc), "kind": kind})
        except Exception as ex:
            r = {"sector": int(sector), "tc": float(tc), "error": repr(ex)[:200]}
        if r:
            res.append(r)
    # combine dchi2 per source over events
    comb = {}
    for r in res:
        for s in r.get("sources", []):
            c = comb.setdefault(s["name"], {"dchi2": 0.0, "n": 0, "needed_eclipse": s["needed_eclipse"], "sep": s["sep"]})
            c["dchi2"] += s["dchi2"]
            c["n"] += 1
    return {"events": res, "combined": comb, "n_sources_tested": len(srcs)}


def vet(tic, times, dur_h=None, depth=None):
    tic = int(tic)
    times = sorted(float(x) for x in times)
    st = ticdb.query([tic], ["ra", "dec", "Tmag", "Teff", "logg", "rad", "mass", "rho", "d", "contratio", "GAIA", "HIP",
                             "Vmag"]).iloc[0].to_dict()
    out = {"tic": tic, "star": {k: (None if (isinstance(v, float) and not math.isfinite(v)) else v) for k, v in st.items()}}
    g = pixels.gaia_cone(st["ra"], st["dec"], 10.0, cols=("source_id", "ra", "dec", "phot_g_mean_mag", "phot_bp_mean_mag",
                                                          "phot_rp_mean_mag", "ruwe", "non_single_star",
                                                          "radial_velocity_error", "rv_nb_transits"))
    if len(g):
        r0 = g.iloc[0]
        out["gaia"] = {"source_id": str(int(r0.source_id)), "G": float(r0.phot_g_mean_mag), "ruwe": float(r0.ruwe),
                       "nss": int(r0.non_single_star),
                       "rv_err": None if r0.radial_velocity_error != r0.radial_velocity_error else float(r0.radial_velocity_error),
                       "rv_n": None if r0.rv_nb_transits != r0.rv_nb_transits else int(r0.rv_nb_transits),
                       "n_within_10as": int(len(g)) - 1}
    out["known"] = {"toi": known.lookup(tic, st["ra"], st["dec"], 120),
                    "toi_within_10arcmin": known.nearby_tois(st["ra"], st["dec"], 600)}
    try:
        out["known"]["tces"] = known.spoc_tces(tic)
    except Exception as ex:
        out["known"]["tces_error"] = repr(ex)[:100]
    lcs, secs = best_lcs_all(tic)
    out["sectors"] = secs
    t, y, sec = flat10(lcs)
    D = (dur_h or 4.0) / 24
    evs = []
    for tc in times:
        r = box_depth(t, y, tc, D)
        lc = sector_of(tc, lcs)
        evs.append({"tc": tc, "sector": None if lc is None else int(lc["meta"]["sector"]),
                    "product": None if lc is None else lc["meta"]["product"],
                    "depth": None if r is None else r[0], "err": None if r is None else r[1]})
    out["events"] = evs
    dep = depth or np.nanmedian([e["depth"] for e in evs if e["depth"] is not None])
    out["depth_used"] = float(dep)
    if len(times) >= 2:
        out["period"] = alias_test(t, y, sec, times, dep, D)
    else:
        flat = flatten_sectors(lcs)
        grid, allowed, pmin = single_period_constraints(times[0], dep, D * 24, flat, pmin=2.0, pmax=3000.0)
        out["period"] = {"single": True, "p_min_allowed": pmin}
    ev_s = [(e["sector"], e["tc"]) for e in evs if e["sector"] is not None]
    out["pixels"] = pixel_tests(tic, ev_s, D, st["ra"], st["dec"], dep, secs)
    return out


def main():
    res = json.load(open(OUT)) if OUT.exists() else {}
    if sys.argv[1] == "--shortlist":
        import report
        tmax = float(sys.argv[2]) if len(sys.argv) > 2 else 10.5
        tb = pd.read_csv(RESULTS / "candidates" / "trackb_final.csv")
        s = report.trackb_shortlist(tb)
        s = s[s.Tmag < tmax].sort_values("Tmag")
        jobs = [(int(r.tic), [float(x) for x in str(r.event_times).split(";") if x], float(r.t14_h), float(r.depth))
                for _, r in s.iterrows()]
    else:
        jobs = []
        for a in sys.argv[1:]:
            tic, ts = a.split(":")
            jobs.append((int(tic), [float(x) for x in ts.split(";")], None, None))
    for tic, times, dur, dep in jobs:
        if str(tic) in res and "error" not in res[str(tic)]:
            continue
        try:
            r = vet(tic, times, dur, dep)
        except Exception as ex:
            r = {"tic": tic, "error": repr(ex)[:300]}
        res[str(tic)] = r
        json.dump(res, open(OUT, "w"), indent=1, default=float)
        p = r.get("period") or {}
        comb = (r.get("pixels") or {}).get("combined", {})
        worst = min(comb.values(), key=lambda c: c["dchi2"])["dchi2"] if comb else None
        allowed = [a["P"] for a in p.get("aliases", []) if a["status"] == "allowed"] if "aliases" in p else None
        print(tic, "| gaia", r.get("gaia", {}).get("ruwe"), "| TOI", len(r.get("known", {}).get("toi", [])),
              "| TCE", len(r.get("known", {}).get("tces", [])), "| allowed P", allowed[:4] if allowed else p.get("p_min_allowed"),
              "| min dchi2 over sources", worst, "| err", r.get("error"), flush=True)


if __name__ == "__main__":
    main()
