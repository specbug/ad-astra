"""Statistical validation of TIC 402710280.01 with TRICERATOPS (Giacalone & Dressing 2020; Giacalone et al. 2021).

TRICERATOPS weighs the transiting-planet hypothesis against every false-positive scenario it knows: an eclipsing
binary on the target (TEB), an unresolved bound or chance-aligned companion that carries the planet or the
eclipsing binary (PTP/PEB/STP/SEB, DTP/DEB/BTP/BEB), and nearby resolved stars as the source (NTP/NEB). Each
scenario's likelihood comes from fitting the phase-folded light curve; its prior from the TIC/Gaia field around the
star, the photometric apertures, and population statistics. FPP is the total probability of the false-positive
scenarios, NFPP that of the nearby-star ones. Giacalone et al. (2021): FPP < 0.015 and NFPP < 0.001 validates a
planet; FPP < 0.5 and NFPP < 0.001 makes it a likely planet.

Inputs used here (v2, 15 October 2026, after the independent audit):
  - sectors 39, 93, 100 and 102. Sector 66 is left out: its predicted transit falls on a stray-light ramp with data
    starting 0.4 h before ingress (audit_402710280.py), so it is not evidence either way. v1 used 39, 66, 93, 100;
  - apertures: TESS-SPOC optimal apertures (sectors 39, 66, from the light-curve files), and circular apertures of
    the radius used for the photometry in sector 93 (QLP; radius from the file header if given, else 2 px) and
    sector 100 (our own, ffi_new.py);
  - the transits' 1.6-day windows (quadratic baseline divided out, as in the fit), each centred on its own measured
    mid-time (they scatter by tens of minutes about a straight-line ephemeris), stacked and binned to 10 minutes;
    depth 2100 ppm (the four transits' mean, 2097 +- 46); no contrast curve (no high-resolution imaging exists).
    v1 folded on P = 36.33972 d with depth 1449 ppm.
The Monte Carlo is repeated and the mean and scatter of FPP and NFPP reported.

  python fpp_triceratops.py [n_runs]
Output: results/candidates/fpp_tic402710280.json, results/report/cases/TIC402710280_triceratops_field.png
"""
import json
import sys
import warnings

import env  # noqa: F401
import numpy as np

from env import RESULTS

warnings.filterwarnings("ignore")

TIC = 402710280
SECTORS = [39, 93, 100, 102]
TCS = [2384.6094, 3838.1847, 4092.5851, 4128.9449]  # BTJD, measured mid-times (audit_402710280.py)
P_ORB = 36.3399
DEPTH = 2100e-6
import os
BOX5 = os.environ.get("FPP_APERTURES") == "box5"   # sensitivity check: TRICERATOPS's default 5x5 pixel apertures
ERRSCALE = float(os.environ.get("FPP_ERRSCALE", "1"))   # sensitivity check: errors inflated for correlated noise
OUT = RESULTS / "candidates" / ("fpp_tic402710280_v2_box5.json" if BOX5 else
                                f"fpp_tic402710280_v2_err{ERRSCALE:g}.json" if ERRSCALE != 1 else "fpp_tic402710280_v2.json")
FIG = RESULTS / "report" / "cases" / "TIC402710280_triceratops_field.png"


def tess_spoc_aperture(sector):
    """Optimal aperture of the TESS-SPOC full-frame light curve, as [[col, row], ...] (file read from the MAST
    bucket on S3, s3io)."""
    import io
    import s3io
    from astropy.io import fits
    raw = s3io.fetch(s3io.key_for("tspoc", TIC, sector))
    if raw is None:
        return None, "no TESS-SPOC file"
    with fits.open(io.BytesIO(raw)) as h:
        ap = h["APERTURE"].data.astype(int)
        c0, r0 = h["APERTURE"].header["CRVAL1P"], h["APERTURE"].header["CRVAL2P"]
    rows, cols = np.nonzero(np.bitwise_and(ap, 2))
    return np.column_stack([cols + c0, rows + r0]), f"TESS-SPOC optimal aperture ({len(rows)} px)"


def qlp_radius(sector):
    """QLP's aperture radius for this star, if its light-curve header records it."""
    import lightkurve as lk
    sr = lk.search_lightcurve(f"TIC {TIC}", author="QLP", sector=sector)
    if len(sr) == 0:
        return None
    lc = sr[0].download()
    for key in ("BESTAP", "APRAD", "AP_RAD", "APERTURE"):
        v = lc.meta.get(key)
        if v is not None:
            try:
                return float(str(v).split(":")[-1])
            except ValueError:
                continue
    return None


def circle(center, r):
    c, rw = center
    cols = np.arange(np.floor(c - r) - 1, np.ceil(c + r) + 2)
    rows = np.arange(np.floor(rw - r) - 1, np.ceil(rw + r) + 2)
    cc, rr = np.meshgrid(cols, rows)
    m = np.hypot(cc - c, rr - rw) <= r
    return np.column_stack([cc[m], rr[m]]).astype(int)


def light_curve():
    import transit_fit as tf
    from fit_new import lcs_for
    lcs = lcs_for(TIC)
    data = tf.windows_from_lcs(lcs, TCS, 0.8, 5.3)
    t, f = [], []
    for (tt, ff, sig, exp_s, tc) in data:
        t.append(tt - tc)
        f.append(ff)
    t, f = np.concatenate(t), np.concatenate(f)
    o = np.argsort(t)
    t, f = t[o], f[o]
    # 10-minute bins
    edges = np.arange(t.min(), t.max() + 1 / 144, 1 / 144)
    idx = np.digitize(t, edges)
    tb = np.array([t[idx == k].mean() for k in np.unique(idx) if (idx == k).sum() >= 2])
    fb = np.array([f[idx == k].mean() for k in np.unique(idx) if (idx == k).sum() >= 2])
    oot = np.abs(tb) > 0.15
    err = float(np.std(fb[oot]))
    return tb, fb, err


def patch_tic_lookup(tr):
    """TRICERATOPS looks the target up by name (MAST's name resolver, not reachable from here); query the same
    TIC cone around the star's catalogued position instead, nearest star first."""
    import astropy.units as u
    import ticdb
    from astropy.coordinates import SkyCoord
    st = ticdb.query([TIC], ["ra", "dec"]).iloc[0]
    query_region = tr.Catalogs.query_region

    def query_object(name, radius=None, catalog="TIC", **kw):
        res = query_region(SkyCoord(float(st.ra), float(st.dec), unit="deg"), radius=radius, catalog=catalog)
        res.sort("dstArcSec")
        assert int(res[0]["ID"]) == TIC, f"nearest TIC source is {res[0]['ID']}"
        return res
    tr.Catalogs.query_object = query_object

    # TESScut cutouts by coordinates (lightkurve's search goes through the same name resolver)
    class _Cut:
        def __init__(self, hdu):
            self.hdu = hdu

    class _Search:
        def __init__(self, coord, sector):
            self.coord, self.sector = coord, sector

        def download_all(self, cutout_size=(22, 22), **kw):
            from astroquery.mast import Tesscut
            return [_Cut(h) for h in Tesscut.get_cutouts(coordinates=self.coord, size=cutout_size[0], sector=int(self.sector))]

    tr.lightkurve.search_tesscut = lambda target=None, sector=None, **kw: _Search(target, sector)


def main(n_runs=5):
    import matplotlib
    matplotlib.use("Agg")
    from triceratops import triceratops as tr
    patch_tic_lookup(tr)
    target = tr.target(ID=TIC, sectors=np.array(SECTORS))
    aps, notes = [], {}
    for i, s in enumerate(target.sectors):
        ap, note = (None, None)
        if s == 39:
            try:
                ap, note = tess_spoc_aperture(s)
            except Exception as ex:
                ap, note = None, f"TESS-SPOC aperture failed: {ex!r}"[:120]
        if ap is None:
            if s == 93:
                r = None
                try:
                    r = qlp_radius(s)
                except Exception:
                    pass
                r = r or 2.0
                note = (note + "; " if note else "") + f"circle r = {r} px (QLP)"
            else:
                meta = json.loads(str(np.load(f"/home/user/work/cache/ffi_lc/{TIC}_s{s}.npz", allow_pickle=True)["meta"]))
                r = float(meta["aperture_radius_px"])
                note = (note + "; " if note else "") + f"circle r = {r} px (our photometry)"
            ap = circle(target.pix_coords[i][0], r)
        aps.append(ap)
        notes[int(s)] = note
        print("sector", s, note, len(ap), "pixels", flush=True)
    try:
        target.plot_field(sector=SECTORS[-1], ap_pixels=aps[-1], save=True, fname=str(FIG.with_suffix("")))
    except Exception as ex:
        print("field plot failed:", repr(ex)[:120])
    if BOX5:
        aps = None
        notes = {int(x): "5x5 pixels around the target (TRICERATOPS default)" for x in target.sectors}
    target.calc_depths(tdepth=DEPTH, all_ap_pixels=aps)
    near = target.stars[["ID", "Tmag", "sep (arcsec)", "fluxratio", "tdepth"]].head(12)
    print(near.to_string(), flush=True)
    t, f, err = light_curve()
    err *= ERRSCALE
    print("light curve:", len(t), "10-min bins, error", round(err * 1e6), "ppm", flush=True)
    runs = []
    for k in range(n_runs):
        target.calc_probs(time=t, flux_0=f, flux_err_0=err, P_orb=P_ORB, contrast_curve_file=None,
                          exptime=600 / 86400, N=1000000, parallel=True, verbose=0)
        p = target.probs
        runs.append({"FPP": float(target.FPP), "NFPP": float(target.NFPP),
                     "scenarios": {f"{r.scenario}@{r.ID}": float(r.prob) for r in p.itertuples() if r.prob > 1e-4}})
        print(f"run {k + 1}: FPP = {target.FPP:.4f}, NFPP = {target.NFPP:.2e}", flush=True)
        json.dump({"tic": TIC, "sectors": SECTORS, "apertures": notes, "depth": DEPTH, "P": P_ORB,
                   "n_bins": len(t), "flux_err": err, "contrast_curve": None,
                   "nearby": near.to_dict("records"), "runs": runs,
                   "FPP_mean": float(np.mean([r["FPP"] for r in runs])),
                   "FPP_std": float(np.std([r["FPP"] for r in runs])),
                   "NFPP_mean": float(np.mean([r["NFPP"] for r in runs])),
                   "NFPP_std": float(np.std([r["NFPP"] for r in runs]))},
                  open(OUT, "w"), indent=1, default=float)


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 5)
