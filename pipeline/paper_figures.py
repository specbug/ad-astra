"""Figures for the TIC 402710280.01 preprint (paper/), from the stored analysis products.

  Figure 1  the four transits used (S39 TESS-SPOC, S93 QLP, S100 and S102 our FFI photometry), each centred on its
            measured mid-time, 20-minute bins, with the shape model; and the sector-66 window at its predicted time.
  Figure 2  sector 66: raw (SAP) and corrected (PDCSAP) TESS-SPOC flux and background around the predicted transit.
  Figure 3  transit times against the linear ephemeris (O-C), with the 1-sigma prediction band.
  Figure 4  the diluted-eclipse test for the Gaia neighbour: delta-chi2 of eclipse models against the planet model.
  Figure 5  the sector-39 difference image, and the centroid offsets of the four transits.
  Figure 6  the 3.78-day signal: folded light curve with its transit model, and the depth in each sector.

  python paper_figures.py
Output: paper/figures/fig{1..6}.pdf (and .png)
"""
import json
import math

import env  # noqa: F401
import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

import transit_fit as tf  # noqa: E402
from env import RESULTS  # noqa: E402

ROOT = RESULTS.parent
OUT = ROOT / "paper" / "figures"
plt.rcParams.update({"font.family": "serif", "font.serif": ["STIXGeneral"], "mathtext.fontset": "stix",
                     "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8, "xtick.labelsize": 7,
                     "ytick.labelsize": 7, "legend.fontsize": 7, "axes.linewidth": 0.6, "xtick.direction": "in",
                     "ytick.direction": "in", "xtick.top": True, "ytick.right": True, "xtick.major.width": 0.6,
                     "ytick.major.width": 0.6, "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
                     "pdf.fonttype": 42, "ps.fonttype": 42})
INK, DATA, MODEL, GREY = "#1a1a1a", "#3b5b92", "#b0413e", "#9a9a9a"
COL1, COL2 = 3.4, 7.1


def save(fig, name):
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / f"{name}.pdf")
    fig.savefig(OUT / f"{name}.png", dpi=200)
    plt.close(fig)


def binned(x, y, w):
    edges = np.arange(x.min(), x.max() + w, w)
    k = np.digitize(x, edges)
    xb, yb, eb = [], [], []
    for i in np.unique(k):
        s = k == i
        if s.sum() >= 2:
            xb.append(x[s].mean())
            yb.append(y[s].mean())
            eb.append(y[s].std(ddof=1) / math.sqrt(s.sum()))
    return np.array(xb), np.array(yb), np.array(eb)


def fig1(fit, lcs):
    f = fit["TIC402710280"]
    times = {t["sector"]: t["tc"] for t in f["times"]}
    P = f["period_d"][0]
    k, b, aR, u = f["rp_rs"][0], f["b"][0], f["a_rs"][0], f["u"]
    inc = math.degrees(math.acos(b / aR))
    s66 = f["tref_btjd"][0] + (20 - f["ref_epoch"]) * P
    panels = [(39, times[39], "Sector 39\nTESS-SPOC, 600 s", True), (93, times[93], "Sector 93\nQLP, 200 s", True),
              (100, times[100], "Sector 100\nthis work, 200 s", True), (102, times[102], "Sector 102\nthis work, 200 s", True),
              (66, s66, "Sector 66 (not used)\nTESS-SPOC, 200 s", False)]
    by = {lc["meta"]["sector"]: lc for lc in lcs}
    fig, axes = plt.subplots(2, 5, figsize=(COL2, 2.5), sharex=True, sharey="row",
                             gridspec_kw={"height_ratios": [3, 1], "hspace": 0.05, "wspace": 0.06})
    for j, (s, tc, lab, used) in enumerate(panels):
        ax, rx_ = axes[0, j], axes[1, j]
        w = tf.windows_from_lcs([by[s]], [tc], 0.8, 5.3)[0]
        t, fl = w[0], w[1]
        x = (t - tc) * 24
        sel = np.abs(x) < 9
        m = tf._model(t, tc, P, k, aR, inc, u, exp_s=w[3])
        ax.plot(x[sel], (fl[sel] - 1) * 1e6, ".", ms=0.8, color=GREY if used else "#c9c9c9", rasterized=True, zorder=1)
        xb, yb, eb = binned(x[sel], (fl[sel] - 1) * 1e6, 20 / 60)
        ax.errorbar(xb, yb, eb, fmt="o", ms=2.2, lw=0.6, color=DATA if used else GREY, zorder=3)
        xm = np.linspace(-9, 9, 800)
        mm = tf._model(tc + xm / 24, tc, P, k, aR, inc, u, exp_s=w[3])
        ax.plot(xm, (mm - 1) * 1e6, lw=0.9, color=MODEL if used else GREY, zorder=4, ls="-" if used else "--")
        xr, yr, er = binned(x[sel], (fl[sel] - m[sel]) * 1e6, 20 / 60)
        rx_.errorbar(xr, yr, er, fmt="o", ms=1.8, lw=0.5, color=DATA if used else GREY)
        rx_.axhline(0, color=MODEL if used else GREY, lw=0.7)
        ax.set_title(lab, fontsize=6.8)
        rx_.set_xlim(-9, 9)
        rx_.set_xticks([-6, 0, 6])
        if not used:
            for a in (ax, rx_):
                a.axvline(-3.18, color=INK, lw=0.5, ls=":")
            ax.text(-3.5, -2300, "good data\nbegin", fontsize=5.5, ha="right", va="top")
    axes[0, 0].set_ylim(-3600, 2200)
    axes[1, 0].set_ylim(-900, 900)
    axes[1, 0].set_yticks([-500, 0, 500])
    axes[0, 0].set_ylabel("Relative flux (ppm)")
    axes[1, 0].set_ylabel("Resid.", fontsize=7)
    fig.supxlabel("Hours from mid-transit", fontsize=8, y=-0.03)
    save(fig, "fig1")


def fig2():
    import io
    import s3io
    from astropy.io import fits
    d = fits.open(io.BytesIO(s3io.fetch(s3io.key_for("tspoc", 402710280, 66))))[1].data
    t = np.asarray(d["TIME"], float)
    tc = 3111.40
    w = np.isfinite(t) & (np.abs(t - tc) < 0.7)
    h = (t[w] - tc) * 24
    q = np.asarray(d["QUALITY"])[w]
    sap, pdc, bkg = (np.asarray(d[c], float)[w] for c in ("SAP_FLUX", "PDCSAP_FLUX", "SAP_BKG"))
    q0 = q == 0
    sl = (q & 4096) > 0
    ref_s, ref_p = np.nanmedian(sap[q0 & (h > 10)]), np.nanmedian(pdc[q0 & (h > 10)])
    fig, ax = plt.subplots(3, 1, figsize=(COL1, 3.6), sharex=True)
    ax[0].plot(h[sl], (sap[sl] / ref_s - 1) * 1e6, ".", ms=1, color="#c9c9c9", label="scattered light (4096)")
    ax[0].plot(h[~q0 & ~sl], (sap[~q0 & ~sl] / ref_s - 1) * 1e6, ".", ms=1, color="#e3a3a1", label="manual exclude (128)")
    ax[0].plot(h[q0], (sap[q0] / ref_s - 1) * 1e6, ".", ms=1, color=INK, label="quality 0")
    ax[0].set_ylabel("SAP (ppm)")
    ax[0].legend(loc="lower left", bbox_to_anchor=(0, 1.0), ncol=3, frameon=False, markerscale=4, fontsize=6,
                 handletextpad=0.1, columnspacing=0.8)
    ax[1].plot(h[q0], (pdc[q0] / ref_p - 1) * 1e6, ".", ms=1, color=INK)
    ax[1].set_ylabel("PDCSAP (ppm)")
    ax[2].plot(h, bkg / 1e3, "-", lw=0.8, color=MODEL)
    ax[2].set_ylabel(r"Bkg. ($10^3$ e$^-$ s$^{-1}$)")
    ax[2].set_xlabel("Hours from predicted mid-transit (BTJD 3111.40)")
    for a in ax:
        a.axvspan(-2.64, 2.64, color=DATA, alpha=0.10, lw=0)
        for v in (-2.64, 2.64):
            a.axvline(v, color=DATA, lw=0.6, ls="--")
    save(fig, "fig2")


def fig3(fit, ex):
    f = fit["TIC402710280"]
    ref = f["ref_epoch"]
    eP, eT = f["period_d"][1], f["tref_btjd"][1]
    rms = math.sqrt(np.mean([t["o_c_min"] ** 2 for t in f["times"]]))
    fig, ax = plt.subplots(figsize=(COL1, 2.2))
    n = np.linspace(-5, 64, 300)
    band = np.sqrt(eT ** 2 + ((n - ref) * eP) ** 2) * 1440
    ax.fill_between(n, -band, band, color=DATA, alpha=0.18, lw=0, label="ephemeris, 1σ")
    tot = np.sqrt(band ** 2 + rms ** 2)
    ax.plot(n, tot, ":", color=DATA, lw=0.7, label="incl. observed scatter")
    ax.plot(n, -tot, ":", color=DATA, lw=0.7)
    for t in f["times"]:
        mk = {"tspoc": "s", "qlp": "D", "ffi": "o"}[t["product"]]
        ax.errorbar(t["epoch"], t["o_c_min"], t["err_min"], fmt=mk, ms=3.5, color=INK, lw=0.7, mfc="white" if t["product"] == "ffi" else INK)
        ax.annotate(f"S{t['sector']}", (t["epoch"], t["o_c_min"]), textcoords="offset points", xytext=(4, 3), fontsize=6)
    for e in (57, 59):
        ax.axvline(e, color=MODEL, lw=0.6, ls="--")
    ax.text(57 - 0.8, -56, "S114", fontsize=5.5, color=MODEL, ha="right", va="bottom", rotation=90)
    ax.text(59 + 0.9, -56, "La Silla", fontsize=5.5, color=MODEL, ha="left", va="bottom", rotation=90)
    ax.axhline(0, color=GREY, lw=0.5)
    ax.set_xlim(-5, 64)
    ax.set_ylim(-60, 60)
    ax.set_xlabel("Epoch (cycles from sector 39)")
    ax.set_ylabel("O − C (min)")
    ax.legend(loc="lower left", frameon=False, fontsize=6)
    save(fig, "fig3")


def fig4():
    bl = json.load(open(RESULTS / "candidates" / "blend_shape_402710280.json"))
    ks = sorted({r["k"] for r in bl["eb"]})
    lo = [min(r["dchi2"] for r in bl["eb"] if r["k"] == k) for k in ks]
    hi = [max(r["dchi2"] for r in bl["eb"] if r["k"] == k) for k in ks]
    fig, ax = plt.subplots(figsize=(COL1, 1.6))
    ax.axvspan(0.25, 1.2, color=GREY, alpha=0.12, lw=0)
    ax.text(0.55, 12, "main-sequence pairs", fontsize=5.5, ha="center", color=INK)
    ax.fill_between(ks, lo, hi, color=DATA, alpha=0.25, lw=0, label="range over limb darkening")
    ax.plot(ks, lo, "-o", ms=2.5, lw=0.8, color=INK, mfc="white", label="best case for the blend")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(0.25, 12)
    ax.set_ylim(8, 3000)
    ax.set_xticks([0.3, 0.5, 1, 2, 5, 10])
    ax.set_xticklabels(["0.3", "0.5", "1", "2", "5", "10"])
    ax.set_xlabel(r"Radius ratio of the eclipsing companion, $k$")
    ax.set_ylabel(r"$\Delta\chi^2$ (EB $-$ planet)")
    ax.legend(frameon=False, fontsize=6, loc="upper right")
    save(fig, "fig4")


def fig5(fit):
    """Difference image of the sector-39 transit (the highest signal-to-noise one) and the centroid offsets."""
    import pixels
    from matplotlib.patches import Circle
    ra, dec = 248.5242034, -63.5496187                  # Gaia DR3, epoch 2016.0
    pm_t, pm_n = (59.909, 13.685), (4.152, 2.095)
    nb = (248.5227018, -63.5501627)
    tc = {t["sector"]: t["tc"] for t in fit["TIC402710280"]["times"]}[39]
    ep = 2016.0 + (tc + 2457000 - 2457389.0) / 365.25
    dt = ep - 2016.0

    def at(p, pm):
        return (p[0] + pm[0] * dt / 3.6e6 / math.cos(math.radians(p[1])), p[1] + pm[1] * dt / 3.6e6)
    tpos, npos = at((ra, dec), pm_t), at(nb, pm_n)
    path = pixels.cutout(tpos[0], tpos[1], 39, 2, 4, size=11)
    try:
        res, im = pixels.diff_image(path, tc, 5.27, ra=tpos[0], dec=tpos[1])
    finally:
        import os
        os.remove(path)
    wcs = im["wcs"]
    tx, ty = wcs.all_world2pix([tpos], 0)[0]
    nx_, ny_ = wcs.all_world2pix([npos], 0)[0]
    cen = json.load(open(RESULTS.parent / "release" / "tic402710280_zenodo" / "data" / "results" / "centroids.json"))
    fig, ax = plt.subplots(1, 3, figsize=(COL2, 2.25), gridspec_kw={"width_ratios": [1, 1, 1.05], "wspace": 0.55})
    for a, img, lab in ((ax[0], im["direct"], "Out of transit (e$^-$ s$^{-1}$)"),
                        (ax[1], im["diff"], "Out minus in transit (e$^-$ s$^{-1}$)")):
        ext = (-0.5, img.shape[1] - 0.5, -0.5, img.shape[0] - 0.5)
        h = a.imshow(img, origin="lower", cmap="Greys" if a is ax[0] else "RdBu_r", extent=ext,
                     vmin=None if a is ax[0] else -np.nanmax(np.abs(img)), vmax=None if a is ax[0] else np.nanmax(np.abs(img)))
        a.plot(tx, ty, "+", color=MODEL, ms=7, mew=1.0)
        a.plot(nx_, ny_, "x", color="#2b8a3e", ms=4, mew=0.8)
        a.set_xlabel("Column (pixels)")
        a.set_title(lab, fontsize=6.8)
        cb = fig.colorbar(h, ax=a, fraction=0.046, pad=0.02)
        cb.ax.tick_params(labelsize=5.5)
    ax[0].set_ylabel("Row (pixels)")
    c = cen["combined"]["centroid"]
    a = ax[2]
    a.plot(0, 0, "+", color=MODEL, ms=9, mew=1.2, label="TIC 402710280")
    off = [(nb[0] - ra) * math.cos(math.radians(dec)) * 3600 + (pm_n[0] - pm_t[0]) * (2024.775 - 2016) / 1000,
           (nb[1] - dec) * 3600 + (pm_n[1] - pm_t[1]) * (2024.775 - 2016) / 1000]
    a.plot(off[0], off[1], "x", color="#2b8a3e", ms=6, mew=1.0, label="Gaia DR3 neighbor ($G$ = 16.15)")
    for ps in c["per_sector"]:
        a.errorbar(ps["dra"], ps["ddec"], xerr=ps["err"], yerr=ps["err"], fmt="o", ms=2.5,
                   color=GREY, lw=0.6, mfc="white")
        off_ = {39: (3, -8), 93: (2, 4), 100: (3, -8), 102: (-4, 4)}.get(ps["sector"], (3, 3))
        a.annotate(f"S{ps['sector']}", (ps["dra"], ps["ddec"]), textcoords="offset points", xytext=off_, fontsize=5,
                   color="#6f6f6f", ha="left" if off_[0] > 0 else "right")
    a.plot(c["dra"], c["ddec"], "o", ms=3.5, color=DATA, label="combined centroid")
    for r_, ls in ((1, "-"), (3, ":")):
        a.add_patch(Circle((c["dra"], c["ddec"]), r_ * c["err"], fill=False, ec=DATA, lw=0.7, ls=ls))
    a.set_xlim(7, -7)
    a.set_ylim(-7, 7)
    a.set_aspect("equal")
    a.set_xlabel(r"$\Delta\alpha\cos\delta$ (arcsec, east to the left)")
    a.set_ylabel(r"$\Delta\delta$ (arcsec)")
    a.set_title("Difference-image centroids", fontsize=6.8)
    a.legend(loc="upper center", bbox_to_anchor=(0.5, -0.24), frameon=False, fontsize=5.3, handletextpad=0.3,
             ncol=1)
    save(fig, "fig5")
    return res


def fig6():
    """The 3.78-day signal: folded light curve with the transit model, and the depth in each sector."""
    from second_signal_402710280 import load, TEFF
    ss = json.load(open(RESULTS / "candidates" / "second_signal_402710280.json"))
    P, t0 = ss["bls_joint"]["P"], ss["bls_joint"]["t0_btjd"]
    _, rows = load(P, t0)
    x = np.concatenate([((r[1] - t0 + P / 2) % P) - P / 2 for r in rows]) * 24
    y = (np.concatenate([r[2] for r in rows]) - 1) * 1e6
    m = np.abs(x) < 6
    xb, yb, eb = binned(x[m], y[m], 0.25)
    sf = ss["fit_free_rho"]
    k, b, rho, dt = sf["rp_rs"][0], sf["b"][0], sf["rho_solar"][0], sf["dt_min"][0]
    aR = (rho * tf.RHO_SUN * tf.G * (P * 86400) ** 2 / (3 * np.pi)) ** (1 / 3)
    xm = np.linspace(-6, 6, 600)
    mm = tf._model(xm / 24, dt / 1440, P, k, aR, math.degrees(math.acos(b / aR)), tf.tess_ld(TEFF), exp_s=200.0)
    fig, ax = plt.subplots(2, 1, figsize=(COL1, 3.3), gridspec_kw={"height_ratios": [1.5, 1], "hspace": 0.45})
    ax[0].errorbar(xb, yb, eb, fmt="o", ms=2.2, lw=0.6, color=DATA)
    ax[0].plot(xm, (mm - 1) * 1e6, color=MODEL, lw=0.9)
    ax[0].axhline(0, color=GREY, lw=0.5)
    ax[0].set_xlim(-6, 6)
    ax[0].set_xlabel("Hours from mid-transit")
    ax[0].set_ylabel("Relative flux (ppm)")
    sec = ss["sectors_depth"]
    xs = np.arange(len(sec))
    ax[1].errorbar(xs, [r["depth_ppm"] for r in sec], [r["err_ppm"] for r in sec], fmt="o", ms=3, color=INK, lw=0.7)
    ax[1].axhline(ss["depth_all"][0], color=MODEL, lw=0.8)
    ax[1].axhline(0, color=GREY, lw=0.5)
    ax[1].set_xticks(xs)
    ax[1].set_xticklabels([str(r["sector"]) for r in sec])
    ax[1].set_xlabel("Sector")
    ax[1].set_ylabel("Depth (ppm)")
    ax[1].set_ylim(-100, 800)
    save(fig, "fig6")


def main():
    from lcaccess import lcs_for
    fit = json.load(open(RESULTS / "report" / "transit_fits.json"))
    ex = json.load(open(RESULTS / "report" / "extras.json"))["TIC402710280"]
    lcs = [lc for lc in lcs_for(402710280) if lc["time"].size > 100]
    fig1(fit, lcs)
    fig2()
    fig3(fit, ex)
    fig4()
    fig5(fit)
    fig6()
    print("figures written to", OUT)


if __name__ == "__main__":
    main()
