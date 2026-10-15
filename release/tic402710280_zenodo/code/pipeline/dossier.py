"""Candidate dossiers: one multi-panel figure plus a JSON of every vetting number, per candidate.

make_single_transit_dossier(tic, events)   Track B / D: events = list of (product, sector, t0, dur_h)
"""
import json

import env  # noqa: F401
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import s3io
import ticdb
import known
import pixels
from lcprep import bin_lc, detrend_fast
from vet_events import vet_one, _trap
from env import RESULTS, CATALOGS

DOSSIERS = RESULTS / "dossiers"


def all_sectors(tic):
    import pandas as pd
    out = {}
    for p in ("spoc", "tspoc", "qlp"):
        m = pd.read_parquet(CATALOGS / f"{p}_manifest.parquet", columns=["tic", "sector"],
                            filters=[("tic", "==", int(tic))])
        out[p] = sorted(m.sector.astype(int).tolist())
    return out


def best_lcs(tic, secs):
    """One light curve per sector: SPOC 2-min > TESS-SPOC > QLP."""
    got, lcs = set(), []
    for p in ("spoc", "tspoc", "qlp"):
        for s in secs.get(p, []):
            if s in got:
                continue
            try:
                lc = s3io.load(p, tic, s)
            except Exception:
                lc = None
            if lc is not None and lc["time"].size > 100:
                lcs.append(lc)
                got.add(s)
    return sorted(lcs, key=lambda l: l["meta"]["sector"])


def make_single_transit_dossier(tic, events, label="", pixel=True):
    DOSSIERS.mkdir(parents=True, exist_ok=True)
    st = ticdb.query([tic], ["ra", "dec", "Tmag", "Teff", "logg", "rad", "mass", "d", "contratio", "lumclass",
                             "GAIA", "HIP"]).iloc[0].to_dict()
    secs = all_sectors(tic)
    lcs = best_lcs(tic, secs)
    info = {"tic": int(tic), "star": {k: (None if (isinstance(v, float) and np.isnan(v)) else v) for k, v in st.items()},
            "sectors": secs, "known": known.lookup(tic, st["ra"], st["dec"]), "events": []}
    n_ev = len(events)
    fig = plt.figure(figsize=(14, 4 + 2.2 * (n_ev + 1) + 1.6 * int(np.ceil(len(lcs) / 3))))
    gs = fig.add_gridspec(2 + n_ev + int(np.ceil(len(lcs) / 3)), 3)
    # per-event zooms with trapezoid fits + difference images
    for i, (prod, sector, t0, dur_h) in enumerate(events):
        lc = s3io.load(prod, tic, sector)
        v = vet_one(lc, t0, dur_h)
        ev = {"product": prod, "sector": int(sector), "t0": float(t0), "dur_h": float(dur_h), **v}
        t, f = lc["time"], lc["flux"]
        o = np.argsort(t)
        t, f = t[o], f[o]
        D = dur_h / 24
        w = max(3 * D, 0.6)
        m = np.abs(t - t0) < w
        ax = fig.add_subplot(gs[i, 0:2])
        ax.plot(t[m] - t0, f[m], ".", ms=2, color="0.7")
        tb, fb, _, _ = bin_lc(t[m], f[m], 20.0)
        ax.plot(tb - t0, fb, "o", ms=3, color="k")
        if "fit_depth" in v:
            x = np.linspace(-w, w, 800)
            ax.plot(x, np.median(f[m & (np.abs(t - t0) > D)]) * (1 - _trap(x, v["fit_t0"] - t0, v["fit_depth"], v["fit_t14_h"] / 24,
                                                                       v["fit_ingress_frac"])), "r-", lw=1)
        ax.set_title(f"S{sector} {prod}  t0={t0:.4f}  depth={v.get('fit_depth', np.nan) * 1e3:.2f} ppt  "
                     f"T14={v.get('fit_t14_h', np.nan):.2f} h  ingress={v.get('fit_ingress_frac', np.nan):.2f}  "
                     f"dip/bump={v.get('snr_here', np.nan):.0f}/{v.get('max_bump_snr', np.nan):.1f}", fontsize=8)
        ax.set_xlabel("days from event", fontsize=8)
        axp = fig.add_subplot(gs[i, 2])
        if pixel:
            try:
                path = pixels.cutout(st["ra"], st["dec"], int(sector), lc["meta"]["camera"], lc["meta"]["ccd"], 15)
                pr, imgs = pixels.diff_image(path, v.get("fit_t0", t0), v.get("fit_t14_h", dur_h), st["ra"], st["dec"])
                ev.update({f"pix_{k}": v_ for k, v_ in pr.items()})
                if imgs is not None:
                    axp.imshow(imgs["diff"], origin="lower", cmap="viridis")
                    axp.plot(*imgs["target"], "r+", ms=12, mew=2)
                    g = pixels.gaia_cone(st["ra"], st["dec"], 160)
                    g = g[g.phot_g_mean_mag < st["Tmag"] + 7]
                    xy = imgs["wcs"].all_world2pix(np.c_[g.ra.values, g.dec.values], 0)
                    axp.scatter(xy[:, 0], xy[:, 1], s=np.clip(60 - 6 * (g.phot_g_mean_mag - g.phot_g_mean_mag.min()), 4, 60),
                                facecolors="none", edgecolors="w", lw=0.6)
                    axp.plot(pr["cx_diff"], pr["cy_diff"], "x", color="orange", ms=8)
                    axp.set_xlim(-0.5, imgs["diff"].shape[1] - 0.5)
                    axp.set_ylim(-0.5, imgs["diff"].shape[0] - 0.5)
                    axp.set_title(f"difference image; offset {pr['offset_direct_pix']:.2f} px vs direct", fontsize=8)
                    ev["neighbors_within_60as"] = g[g.sep_arcsec < 60][["source_id", "phot_g_mean_mag", "sep_arcsec"]].to_dict("records")
            except Exception as e:
                axp.text(0.1, 0.5, f"pixel step failed:\n{str(e)[:80]}", fontsize=7)
        info["events"].append(ev)
    # every sector, flattened, events marked
    for j, lc in enumerate(lcs):
        ax = fig.add_subplot(gs[n_ev + 1 + j // 3, j % 3])
        t, f = lc["time"], lc["flux"]
        o = np.argsort(t)
        tb, fb, _, _ = bin_lc(t[o], f[o], 30.0)
        y, _ = detrend_fast(tb, fb, 2.5, knot_step=2.5 / 36)
        ax.plot(tb, y, ".", ms=1.5, color="k")
        for (_, s, t0, _) in events:
            if s == lc["meta"]["sector"]:
                ax.axvline(t0, color="r", alpha=0.4)
        ax.set_title(f"S{lc['meta']['sector']} {lc['meta']['product']}", fontsize=8)
        ax.tick_params(labelsize=6)
    ax = fig.add_subplot(gs[n_ev, :])
    ax.axis("off")
    s = info["star"]
    txt = (f"TIC {tic} {label}   Tmag={s['Tmag']:.2f}  Teff={s['Teff']}  R*={s['rad']}  M*={s['mass']}  d={s['d']} pc  "
           f"contam={s['contratio']}   sectors: " + ", ".join(f"{k}:{v}" for k, v in secs.items() if v) +
           f"\nknown: {[k['source'] + ' ' + k['id'] for k in info['known']] or 'none in TOI lists'}")
    ax.text(0, 0.5, txt, fontsize=8, va="center", family="monospace")
    fig.tight_layout()
    png = DOSSIERS / f"TIC{tic}{('_' + label) if label else ''}.png"
    fig.savefig(png, dpi=80)
    plt.close(fig)
    (DOSSIERS / f"TIC{tic}{('_' + label) if label else ''}.json").write_text(json.dumps(info, indent=1, default=float))
    return png, info


def make_eb_dossier(tic, product, sector, label=""):
    """Track C: EB fold, residuals, second-signal fold, and difference images for each eclipse set."""
    from ebsearch import analyse_eb
    DOSSIERS.mkdir(parents=True, exist_ok=True)
    lc = s3io.load(product, tic, sector)
    res, arr = analyse_eb(lc, return_arrays=True)
    st = ticdb.query([tic], ["ra", "dec", "Tmag", "Teff", "rad", "d", "contratio"])
    star = st.iloc[0].to_dict() if len(st) else {"ra": np.nan, "dec": np.nan, "Tmag": np.nan}
    t, y, tm, rr = arr["t"], arr["y"], arr["template"], arr["resid"]
    P1, t01 = res["P1"], res["t0_1"]
    fig, ax = plt.subplots(3, 2, figsize=(14, 10))
    ax[0, 0].plot(t, y, ".", ms=1.5, color="k"); ax[0, 0].plot(t, tm, "-", color="tab:orange", lw=0.6)
    ax[0, 0].set_title(f"TIC {tic} S{sector} {product}: flattened flux and EB template (P1={P1:.5f} d)", fontsize=9)
    ph1 = ((t - t01) / P1 + 0.25) % 1 - 0.25
    ax[0, 1].plot(ph1, y, ".", ms=1.5, color="k"); ax[0, 1].set_title(f"folded at P1, depth1={res['depth1']:.3f}", fontsize=9)
    ax[1, 0].plot(t, rr, ".", ms=1.5, color="k")
    if "P2" in res:
        P2, t02, D2 = res["P2"], res["t0_2"], res["dur2_h"] / 24
        k = np.arange(np.floor((t.min() - t02) / P2), np.ceil((t.max() - t02) / P2) + 1)
        for tc in t02 + k * P2:
            ax[1, 0].axvline(tc, color="r", alpha=0.3)
        ph2 = ((t - t02) / P2 + 0.25) % 1 - 0.25
        ax[1, 1].plot(ph2 * P2 * 24, rr, ".", ms=1.5, color="k")
        ax[1, 1].set_xlim(-6 * D2 * 24, 6 * D2 * 24)
        ax[1, 1].set_xlabel("hours from P2 eclipse")
        ax[1, 1].set_title(f"residuals folded at P2={P2:.5f} d, depth2={res['depth2']:.4f}, snr2={res['snr2']:.0f}, "
                           f"n_tr2={res['n_tr2']}", fontsize=9)
    ax[1, 0].set_title("residuals (P2 eclipse times in red)", fontsize=9)
    # difference images for both eclipse sets
    pix = {}
    try:
        path = pixels.cutout(star["ra"], star["dec"], int(sector), lc["meta"]["camera"], lc["meta"]["ccd"], 15)
        for j, (Pj, t0j, Dj) in enumerate([(P1, t01, res["dur1_h"] / 24)] + ([(res["P2"], res["t0_2"], res["dur2_h"] / 24)] if "P2" in res else [])):
            # stack every eclipse of this set: use the first in-data epoch as t0 and many epochs via a phase mask
            from astropy.io import fits as _f
            with _f.open(path) as h:
                d = h[1].data
                tt = np.asarray(d["TIME"], float); fl = np.asarray(d["FLUX"], float)
                q = np.asarray(d["QUALITY"], int)
            g = np.isfinite(tt) & (q == 0) & np.isfinite(fl).all(axis=(1, 2))
            tt, fl = tt[g], fl[g]
            ph = np.abs(((tt - t0j + 0.5 * Pj) % Pj) - 0.5 * Pj)
            inn = ph < 0.4 * Dj
            out = (ph > 0.75 * Dj) & (ph < 0.75 * Dj + max(Dj, 0.25))
            if inn.sum() > 2 and out.sum() > 4:
                diff = np.median(fl[out], axis=0) - np.median(fl[inn], axis=0)
                a = ax[2, j]
                a.imshow(diff, origin="lower", cmap="viridis")
                from astropy.wcs import WCS
                with _f.open(path) as h:
                    w = WCS(h[2].header)
                tx, ty = w.all_world2pix([[star["ra"], star["dec"]]], 0)[0]
                a.plot(tx, ty, "r+", ms=12, mew=2)
                py, px = np.unravel_index(np.nanargmax(diff), diff.shape)
                a.plot(px, py, "x", color="orange", ms=10)
                pix[f"P{j + 1}_peak_offset_pix"] = float(np.hypot(px - tx, py - ty))
                a.set_title(f"difference image, eclipse set {j + 1}: peak {pix[f'P{j + 1}_peak_offset_pix']:.1f} px from target", fontsize=9)
    except Exception as e:
        ax[2, 0].text(0.1, 0.5, f"pixel step failed: {str(e)[:80]}", fontsize=8)
    res.update(pix)
    fig.suptitle(f"TIC {tic}  T={star.get('Tmag', np.nan):.2f}  Teff={star.get('Teff')}  R*={star.get('rad')}  "
                 f"known: {[k['id'] for k in known.lookup(tic, star.get('ra'), star.get('dec'))] or 'none'}", fontsize=10)
    fig.tight_layout()
    png = DOSSIERS / f"EB_TIC{tic}_S{sector}{('_' + label) if label else ''}.png"
    fig.savefig(png, dpi=80)
    plt.close(fig)
    res.pop("extra_events", None)
    (DOSSIERS / (png.stem + ".json")).write_text(json.dumps(res, indent=1, default=float))
    return png, res
