"""Difference-image centroids: where on the sky the light that disappears in a transit or eclipse comes from.

A Gaussian point-spread function is fitted to each difference image (out of eclipse minus in eclipse); the offset is
measured from the target's own position on the same image, found by fitting the out-of-eclipse image with Gaussians
at the Gaia DR3 positions moved to the observing date. Per-eclipse offsets are combined within and between sectors,
and compared with every Gaia source bright enough to make the observed depth.
"""
import io
import math

import env  # noqa: F401
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

import pixels
import s3io


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
    r = _diff_image(t, F, tc, D)
    if r is None:
        return None
    return fit_diff(*r, w, tc, ra0, dec0, min_snr=min_snr, gf=gf)


def fit_diff(diff, err, base, w, tc, ra0, dec0, min_snr=5.0, gf=None):
    """The PSF fit of _centroid, for any difference image (one eclipse, or several stacked)."""
    from astropy.wcs.utils import proj_plane_pixel_scales
    from scipy.optimize import least_squares
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


def neighbour_offsets(gf, depth, max_sep=150.0, epoch=2016.0):
    """Gaia stars able to make a dip of this depth (an eclipse of at most 100%), as offsets from the target's
    position, both moved with their proper motions to `epoch` (decimal year; Gaia DR3 positions are for 2016.0).
    The centroids are measured from the target's own image on each frame, so the relevant offsets are those at
    the observing epoch; a fast-moving target changes them (TIC 402710280's 2.3" neighbour of 2000 is 3.1" away in
    2016 and 3.4-3.6" in 2021-2026)."""
    tg = gf[gf.is_target]
    if tg.empty:
        return []
    t = tg.iloc[0]
    dt = epoch - 2016.0

    def pm(r):
        return (0.0 if r.pmra != r.pmra else float(r.pmra)) * dt / 1000, (0.0 if r.pmdec != r.pmdec else float(r.pmdec)) * dt / 1000
    tpa, tpd = pm(t)
    lim = float(t.tmag_est) - 2.5 * math.log10(max(depth, 1e-5))
    out = []
    for _, r in gf.iterrows():
        if r.is_target or r.tmag_est >= lim:
            continue
        rpa, rpd = pm(r)
        dra = (float(r.ra) - float(t.ra)) * math.cos(math.radians(float(t.dec))) * 3600 + rpa - tpa
        ddec = (float(r.dec) - float(t.dec)) * 3600 + rpd - tpd
        sep = math.hypot(dra, ddec)
        if 1.0 < sep < max_sep:
            out.append({"name": f"G{r.phot_g_mean_mag:.1f} at {sep:.0f}\"", "dra": dra, "ddec": ddec})
    return out


def _wmean(X, clip=3.5, min_keep=2):
    """Weighted 2-D mean of rows (dra, ddec, sigma), dropping the worst point while it is > clip sigma away;
    the error is scaled up by sqrt(reduced chi-square) when the points scatter more than their errors say."""
    keep = np.ones(len(X), bool)
    while True:
        w = 1 / X[keep, 2] ** 2
        mx, my = np.sum(w * X[keep, 0]) / w.sum(), np.sum(w * X[keep, 1]) / w.sum()
        d = np.hypot(X[:, 0] - mx, X[:, 1] - my) / X[:, 2]
        worst = int(np.argmax(np.where(keep, d, -1)))
        if keep.sum() <= min_keep or d[worst] < clip:
            break
        keep[worst] = False
    n = int(keep.sum())
    err = 1 / math.sqrt(w.sum())
    red = float(np.sum(d[keep] ** 2)) / (2 * (n - 1)) if n >= 2 else None
    if red:
        err *= math.sqrt(max(red, 1.0))
    return mx, my, err, keep, red


def combine_centroids(cents, floor=2.0, clip=3.5):
    """Combine per-eclipse centroid offsets in two steps. Within a sector, a robust weighted mean of its eclipses
    (errors from the fits, scaled up if they scatter more). Between sectors, a robust weighted mean of the sector
    means, each with a 2" systematic floor added (about 0.1 pixel: errors of the PSF model and of the image
    calibration are shared by all eclipses of a sector, so they do not average down within it). Eclipses whose
    difference image a Gaussian fits very badly are dropped first, with an allowance that grows with
    signal-to-noise (a deep eclipse of a bright star never fits a Gaussian exactly)."""
    c = [x for x in cents if x and "dra" in x and not x["edge"] and x["redchi2"] < max(50.0, 0.002 * x["snr"] ** 2)]
    if not c:
        return None
    secs = sorted({x.get("sector") for x in c}, key=lambda v: (v is None, v))
    S, n_ecl = [], 0
    for sec in secs:
        X = np.array([[x["dra"], x["ddec"], x["err"]] for x in c if x.get("sector") == sec])
        mx, my, err, keep, _ = _wmean(X, clip, min_keep=1 if len(X) < 3 else 2)
        S.append([mx, my, math.hypot(err, floor)])
        n_ecl += int(keep.sum())
    S = np.array(S)
    mx, my, err, keep, red = _wmean(S, clip, min_keep=2 if len(S) >= 3 else len(S))
    return {"dra": float(mx), "ddec": float(my), "err": float(err), "n_eclipses": n_ecl, "n_sectors": len(S),
            "n_sectors_used": int(keep.sum()), "rejected_sectors": [secs[i] for i in np.where(~keep)[0]],
            "reduced_chi2_sectors": red,
            "per_sector": [{"sector": secs[i], "dra": float(S[i, 0]), "ddec": float(S[i, 1]), "err": float(S[i, 2])} for i in range(len(S))]}


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
