"""Pixel-level vetting: which star is actually dimming?

gaia_cone      Gaia DR3 sources around a position, from the HATS copy of Gaia DR3 on the MAST S3 bucket
cutout         TESS full-frame-image cutout around a position for one sector (astrocut on the S3 cubes)
diff_image     mean out-of-event minus mean in-event image; flux-weighted centroids of the difference image
               and of the direct image, with a bootstrap over in-event cadences for the uncertainty
"""
import csv
import io
import os
import urllib.request

import env  # noqa: F401
import numpy as np
import pandas as pd
import healpy as hp
from astropy.io import fits
from astropy.wcs import WCS

import s3parquet
from env import CACHE

GAIA = "gaia/gaia_dr3/public/hats/gaia"
_parts = None


def _gaia_partitions():
    global _parts
    if _parts is None:
        rows = list(csv.reader(io.StringIO(urllib.request.urlopen(
            f"https://stpubdata.s3.amazonaws.com/{GAIA}/partition_info.csv").read().decode())))[1:]
        _parts = {(int(o), int(n)) for o, n in rows}
    return _parts


def gaia_cone(ra, dec, radius_arcsec=180.0,
              cols=("source_id", "ra", "dec", "phot_g_mean_mag", "phot_bp_mean_mag", "phot_rp_mean_mag",
                    "parallax", "ruwe", "non_single_star", "phot_variable_flag")):
    parts = _gaia_partitions()
    vec = hp.ang2vec(ra, dec, lonlat=True)
    pix6 = hp.query_disc(2 ** 6, vec, np.radians(radius_arcsec / 3600) + np.radians(0.05), nest=True, inclusive=True)
    keys = set()
    for p in pix6:
        for o in range(6, 1, -1):
            q = int(p) >> (2 * (6 - o))
            if (o, q) in parts:
                keys.add((o, q))
                break
    out = []
    for o, n in sorted(keys):
        key = f"{GAIA}/dataset/Norder={o}/Dir={(n // 10000) * 10000}/Npix={n}.parquet"
        t = s3parquet.read_columns(key, list(cols), threads=12).to_pandas()
        out.append(t)
    g = pd.concat(out, ignore_index=True) if out else pd.DataFrame(columns=cols)
    ra1, de1, ra2, de2 = map(np.radians, (ra, dec, g.ra.values, g.dec.values))
    c = np.sin(de1) * np.sin(de2) + np.cos(de1) * np.cos(de2) * np.cos(ra1 - ra2)
    g["sep_arcsec"] = np.degrees(np.arccos(np.clip(c, -1, 1))) * 3600
    g = g[g.sep_arcsec < radius_arcsec].sort_values("phot_g_mean_mag").reset_index(drop=True)
    # approximate TESS magnitude from Gaia photometry (Stassun et al. 2019 relation)
    bprp = (g.phot_bp_mean_mag - g.phot_rp_mean_mag).fillna(1.0)
    g["tmag_est"] = g.phot_g_mean_mag - 0.00522555 * bprp ** 3 + 0.0891337 * bprp ** 2 - 0.633923 * bprp + 0.0324473
    return g


def cutout(ra, dec, sector, camera, ccd, size=15):
    from astrocut import CutoutFactory
    from astropy.coordinates import SkyCoord
    outdir = CACHE / "cutouts"
    outdir.mkdir(parents=True, exist_ok=True)
    cube = f"s3://stpubdata/tess/public/mast/tess-s{sector:04d}-{camera}-{ccd}-cube.fits"
    path = CutoutFactory().cube_cut(cube, SkyCoord(ra, dec, unit="deg"), cutout_size=size,
                                    output_path=str(outdir), verbose=False)
    return path


def diff_image(path, t0, dur_h, ra=None, dec=None, n_boot=200, rng=None, box=2):
    """Difference image for one event. Centroids are flux-weighted within (2*box+1)^2 pixels around the
    difference-image peak; offsets are measured from the target's catalogue position (via the WCS)."""
    rng = rng or np.random.default_rng(0)
    with fits.open(path) as h:
        d = h[1].data
        t = np.asarray(d["TIME"], float)
        flux = np.asarray(d["FLUX"], float)
        q = np.asarray(d["QUALITY"], int) if "QUALITY" in d.columns.names else np.zeros(t.size, int)
        wcs = WCS(h[2].header)
    good = np.isfinite(t) & (q == 0) & np.isfinite(flux).all(axis=(1, 2))
    t, flux = t[good], flux[good]
    D = dur_h / 24
    inn = np.abs(t - t0) < 0.4 * D
    out = (np.abs(t - t0) > 0.75 * D) & (np.abs(t - t0) < 0.75 * D + max(D, 0.25))
    res = {"n_in": int(inn.sum()), "n_out": int(out.sum())}
    if inn.sum() < 2 or out.sum() < 4:
        return res, None
    img_out = np.median(flux[out], axis=0)
    diff = img_out - np.median(flux[inn], axis=0)
    ny, nx = diff.shape
    yy, xx = np.mgrid[0:ny, 0:nx]
    if ra is not None:
        tx, ty = wcs.all_world2pix([[ra, dec]], 0)[0]
    else:
        tx, ty = (nx - 1) / 2, (ny - 1) / 2
    py, px = np.unravel_index(np.nanargmax(diff), diff.shape)

    def centroid(img, cx, cy):
        m = (np.abs(xx - cx) <= box) & (np.abs(yy - cy) <= box)
        w = np.clip(img, 0, None) * m
        s = w.sum()
        return (np.sum(w * xx) / s, np.sum(w * yy) / s) if s > 0 else (np.nan, np.nan)

    cx_d, cy_d = centroid(diff, px, py)
    cx_o, cy_o = centroid(img_out - np.median(img_out), px, py)  # same box: cancels pixel-phase bias
    boots = []
    idx_in, idx_out = np.where(inn)[0], np.where(out)[0]
    for _ in range(n_boot):
        bi, bo = rng.choice(idx_in, idx_in.size), rng.choice(idx_out, idx_out.size)
        boots.append(centroid(np.median(flux[bo], axis=0) - np.median(flux[bi], axis=0), px, py))
    boots = np.array(boots)
    err = float(np.sqrt(np.nanvar(boots[:, 0]) + np.nanvar(boots[:, 1])))
    off = float(np.hypot(cx_d - tx, cy_d - ty))
    noise = 1.4826 * np.median(np.abs(diff - np.median(diff))) + 1e-12
    res.update({"target_x": float(tx), "target_y": float(ty), "cx_diff": float(cx_d), "cy_diff": float(cy_d),
                "peak_x": int(px), "peak_y": int(py), "offset_pix": off, "offset_arcsec": off * 21.0,
                "offset_err_pix": err, "offset_sigma": off / err if err > 0 else np.nan,
                "offset_direct_pix": float(np.hypot(cx_d - cx_o, cy_d - cy_o)),
                "diff_snr": float(diff[py, px] / noise),
                "diff_at_target_snr": float(diff[int(round(ty)), int(round(tx))] / noise)})
    return res, {"diff": diff, "direct": img_out, "wcs": wcs, "target": (tx, ty)}
