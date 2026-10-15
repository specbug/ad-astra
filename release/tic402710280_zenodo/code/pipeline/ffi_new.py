"""Light curves from the TESS full-frame images of sectors that have no public FFI light curves yet.

Public FFI light curves end at sector 98 (QLP) and near sector 80 (TESS-SPOC), but the calibrated full-frame
images of sectors 99-108 are on AWS as cubes. For a star on silicon in those sectors and without SPOC 2-min data
there, this makes a light curve from a 21x21-pixel cutout (astrocut):
  - background per cadence: median of the faintest 40% of pixels in the time-median image, away from the target;
  - aperture: the pixels within r of the target (r = 1.0 ... 4.0 px) whose light curve has the lowest 1-hour
    scatter once a 1-day running median is removed;
  - cadences flagged by SPOC (QUALITY != 0), and those whose background jumps far above its running median
    (scattered-light glints), are dropped;
  - times: the cube's TIME column is barycentric-corrected for a reference position on the CCD, not for the star
    (TIMECORR); it is recomputed for the star's position (target_times), as the SPOC and QLP light curves do. The
    geocentre stands in for the spacecraft (TESS stays within 0.003 AU of Earth: error < 1.5 s);
  - output in the s3io.parse_lc format (time, flux, flux_err, sap_flux, bkg, cx, cy, meta), product "ffi".
Cache: cache/ffi_lc/{tic}_s{sector}.npz

  python ffi_new.py targets.csv [n_proc]   (columns tic, ra, dec, Tmag)  -> builds every missing light curve
"""
import concurrent.futures as cf
import json
import os
import sys

import env  # noqa: F401
import numpy as np
import pandas as pd

from env import CATALOGS

NEW_SECTORS = range(99, 109)
CACHE = CATALOGS.parent / "cache" / "ffi_lc"


def coverage(tics, ras, decs):
    """{tic: [(sector, camera, ccd), ...]} for the new sectors, from tess-point."""
    from tess_stars2px import tess_stars2px_function_entry as t2p
    r = t2p(list(range(len(tics))), list(ras), list(decs))
    out = {int(t): [] for t in tics}
    for i, s, cam, ccd in zip(r[0], r[3], r[4], r[5]):
        if int(s) in NEW_SECTORS:
            out[int(tics[int(i)])].append((int(s), int(cam), int(ccd)))
    return out


def spoc_sectors(tic):
    m = pd.read_parquet(CATALOGS / "spoc_manifest.parquet", columns=["tic", "sector"], filters=[("tic", "==", int(tic))])
    return set(m.sector.astype(int))


def _running_median(t, y, w):
    out = np.empty_like(y)
    lo = np.searchsorted(t, t - w / 2)
    hi = np.searchsorted(t, t + w / 2)
    for i in range(t.size):
        out[i] = np.median(y[lo[i]:hi[i]])
    return out


def target_times(t, timecorr, ra, dec):
    """BJD_TDB - 2457000 for the star at (ra, dec), from cube times corrected for the CCD reference position."""
    import astropy.units as u
    from astropy.coordinates import EarthLocation, SkyCoord
    from astropy.time import Time
    tobs = Time(np.asarray(t) - np.asarray(timecorr) + 2457000, format="jd", scale="tdb",
                location=EarthLocation.from_geocentric(0, 0, 0, unit="m"))
    ltt = tobs.light_travel_time(SkyCoord(ra, dec, unit="deg"), "barycentric").to(u.d).value
    return np.asarray(t) - np.asarray(timecorr) + ltt


def retime_cached(lc, tic, ra, dec, sector, camera, ccd):
    """Recompute the times of a light curve cached before target_times existed: TIMECORR from a 3x3 cutout of
    the same cube, matched by cadence time; the cache is rewritten with meta time_ref = "target"."""
    from astropy.io import fits
    import pixels
    path = pixels.cutout(ra, dec, sector, camera, ccd, size=3)
    try:
        with fits.open(path) as h:
            tc, corr = np.asarray(h[1].data["TIME"], float), np.asarray(h[1].data["TIMECORR"], float)
    finally:
        try:
            os.remove(path)
        except OSError:
            pass
    ok = np.isfinite(tc) & np.isfinite(corr)
    tc, corr = tc[ok], corr[ok]
    j = np.clip(np.searchsorted(tc, lc["time"]), 1, tc.size - 1)
    j = np.where(np.abs(tc[j - 1] - lc["time"]) < np.abs(tc[j] - lc["time"]), j - 1, j)
    assert np.max(np.abs(tc[j] - lc["time"])) < 1e-6, "cadences do not match the cube"
    new = target_times(lc["time"], corr[j], ra, dec)
    lc["meta"]["time_ref"] = "target"
    lc["meta"]["time_shift_s"] = float(np.median(new - lc["time"]) * 86400)
    lc["time"] = new
    p = CACHE / f"{int(tic)}_s{sector}.npz"
    np.savez_compressed(p, **{k: v for k, v in lc.items() if k != "meta"}, meta=json.dumps(lc["meta"]))
    return lc


def build(tic, ra, dec, sector, camera, ccd, size=21, target_time=False):
    """Make (or read from cache) the light curve of one star in one new sector. New light curves always carry
    target-corrected times; a cached one made before that is retimed only with target_time=True."""
    CACHE.mkdir(parents=True, exist_ok=True)
    p = CACHE / f"{int(tic)}_s{sector}.npz"
    if p.exists():
        z = np.load(p, allow_pickle=True)
        lc = {k: z[k] for k in z.files if k != "meta"} | {"meta": json.loads(str(z["meta"]))}
        if target_time and lc["meta"].get("time_ref") != "target":
            lc = retime_cached(lc, tic, ra, dec, sector, camera, ccd)
        return lc
    from astropy.io import fits
    from astropy.wcs import WCS
    import pixels
    path = pixels.cutout(ra, dec, sector, camera, ccd, size=size)
    try:
        with fits.open(path) as h:
            d = h[1].data
            t = np.asarray(d["TIME"], float)
            tcorr = np.asarray(d["TIMECORR"], float)
            F = np.asarray(d["FLUX"], float)
            E = np.asarray(d["FLUX_ERR"], float)
            q = np.asarray(d["QUALITY"], int)
            w = WCS(h[2].header)
    finally:
        try:
            os.remove(path)   # cutouts are large; the light curve is what is kept
        except OSError:
            pass
    good = np.isfinite(t) & np.isfinite(tcorr) & (q == 0) & np.isfinite(F).all(axis=(1, 2))
    t, F, E = target_times(t[good], tcorr[good], ra, dec), F[good], E[good]
    x0, y0 = (float(v) for v in w.all_world2pix([[ra, dec]], 0)[0])
    ny, nx = F.shape[1:]
    yy, xx = np.mgrid[0:ny, 0:nx]
    rr = np.hypot(xx - x0, yy - y0)
    med = np.median(F, axis=0)
    # background: faint pixels at least 4 px from the target
    far = rr > 4
    thr = np.percentile(med[far], 40)
    bmask = far & (med <= thr)
    bkg = np.median(F[:, bmask], axis=1)
    Fs = F - bkg[:, None, None]
    # scattered-light glints: background far above its 1-day running median
    rmb = _running_median(t, bkg, 1.0)
    mad = 1.4826 * np.median(np.abs(bkg - rmb)) + 1e-9
    ok = (bkg - rmb) < 8 * mad
    t, Fs, E, bkg = t[ok], Fs[ok], E[ok], bkg[ok]
    best = None
    for r in (1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0):
        ap = rr <= r
        if ap.sum() < 1:
            continue
        f = Fs[:, ap].sum(axis=1)
        if not np.isfinite(f).all() or np.median(f) <= 0:
            continue
        fn = f / np.median(f)
        resid = fn - _running_median(t, fn, 1.0)
        # 1-hour scatter
        nb = max(int(round((1 / 24) / np.median(np.diff(t)))), 1)
        m = resid[: resid.size // nb * nb].reshape(-1, nb).mean(axis=1)
        s = 1.4826 * np.median(np.abs(m - np.median(m)))
        if best is None or s < best[0]:
            best = (s, r, ap, f)
    if best is None:
        return None
    s1h, r, ap, f = best
    fe = np.sqrt((E[:, ap] ** 2).sum(axis=1))
    wts = np.clip(Fs[:, ap], 0, None)
    cx = (wts * xx[ap]).sum(axis=1) / np.maximum(wts.sum(axis=1), 1e-9)
    cy = (wts * yy[ap]).sum(axis=1) / np.maximum(wts.sum(axis=1), 1e-9)
    medf = np.median(f)
    out = {"time": t, "flux": f / medf, "flux_err": fe / medf, "sap_flux": f / medf, "bkg": bkg, "cx": cx, "cy": cy,
           "meta": {"tic": int(tic), "sector": int(sector), "camera": int(camera), "ccd": int(ccd), "product": "ffi",
                    "cadence_s": float(np.median(np.diff(t)) * 86400), "aperture_radius_px": r,
                    "n_aperture_px": int(ap.sum()), "scatter_1h_ppm": float(s1h * 1e6), "time_ref": "target"}}
    np.savez_compressed(p, **{k: v for k, v in out.items() if k != "meta"}, meta=json.dumps(out["meta"]))
    return out


def new_lcs(tic, ra, dec, cov=None, build_missing=True, target_time=False):
    """Our light curves of the new sectors for one star (only sectors without SPOC 2-min data)."""
    cov = cov if cov is not None else coverage([tic], [ra], [dec])[int(tic)]
    have = spoc_sectors(tic)
    out = []
    for (s, cam, ccd) in cov:
        if s in have:
            continue
        p = CACHE / f"{int(tic)}_s{s}.npz"
        if not p.exists() and not build_missing:
            continue
        try:
            lc = build(tic, ra, dec, s, cam, ccd, target_time=target_time)
        except Exception as ex:
            print("ffi", tic, s, "failed:", repr(ex)[:120], flush=True)
            lc = None
        if lc is not None and lc["time"].size > 100:
            out.append(lc)
    return out


def _job(a):
    tic, ra, dec, s, cam, ccd = a
    try:
        lc = build(tic, ra, dec, s, cam, ccd)
        return tic, s, (lc["meta"]["scatter_1h_ppm"] if lc else None), None
    except Exception as ex:
        return tic, s, None, repr(ex)[:160]


def main():
    d = pd.read_csv(sys.argv[1])
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 4
    cov = coverage(d.tic.astype(int).tolist(), d.ra.tolist(), d.dec.tolist())
    jobs = []
    for _, r in d.iterrows():
        have = spoc_sectors(int(r.tic))
        for (s, cam, ccd) in cov[int(r.tic)]:
            if s not in have and not (CACHE / f"{int(r.tic)}_s{s}.npz").exists():
                jobs.append((int(r.tic), float(r.ra), float(r.dec), s, cam, ccd))
    print(len(d), "stars;", len(jobs), "light curves to make", flush=True)
    with cf.ProcessPoolExecutor(n) as ex:
        for i, (tic, s, sc, err) in enumerate(ex.map(_job, jobs)):
            print(i, tic, s, "1h scatter ppm" if sc else "", round(sc) if sc else err, flush=True)


if __name__ == "__main__":
    main()
