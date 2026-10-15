"""Figures for the TIC 402710280.01 preprint (paper/), from the stored analysis products.

  Figure 1  the four transits used (S39 TESS-SPOC, S93 QLP, S100 and S102 our FFI photometry), each centred on its
            measured mid-time, 20-minute bins, with the shape model; and the sector-66 window at its predicted time.
  Figure 2  sector 66: raw (SAP) and corrected (PDCSAP) TESS-SPOC flux and background around the predicted transit.
  Figure 3  transit times against the linear ephemeris (O-C), with the 1-sigma prediction band.
  Figure 4  the diluted-eclipse test for the Gaia neighbour: delta-chi2 of eclipse models against the planet model.

  python paper_figures.py
Output: paper/figures/fig{1..4}.pdf (and .png)
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
plt.rcParams.update({"font.family": "serif", "font.serif": ["DejaVu Serif"], "mathtext.fontset": "dejavuserif",
                     "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8, "xtick.labelsize": 7,
                     "ytick.labelsize": 7, "legend.fontsize": 7, "axes.linewidth": 0.6, "xtick.direction": "in",
                     "ytick.direction": "in", "xtick.top": True, "ytick.right": True, "xtick.major.width": 0.6,
                     "ytick.major.width": 0.6, "savefig.bbox": "tight", "savefig.pad_inches": 0.02})
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
    fig, axes = plt.subplots(1, 5, figsize=(COL2, 1.9), sharey=True)
    for ax, (s, tc, lab, used) in zip(axes, panels):
        w = tf.windows_from_lcs([by[s]], [tc], 0.8, 5.3)[0]
        t, fl = w[0], w[1]
        x = (t - tc) * 24
        sel = np.abs(x) < 9
        ax.plot(x[sel], (fl[sel] - 1) * 1e6, ".", ms=0.8, color=GREY if used else "#c9c9c9", rasterized=True, zorder=1)
        xb, yb, eb = binned(x[sel], (fl[sel] - 1) * 1e6, 20 / 60)
        ax.errorbar(xb, yb, eb, fmt="o", ms=2.2, lw=0.6, color=DATA if used else GREY, zorder=3)
        xm = np.linspace(-9, 9, 800)
        m = tf._model(tc + xm / 24, tc, P, k, aR, inc, u, exp_s=w[3])
        ax.plot(xm, (m - 1) * 1e6, "-", lw=0.9, color=MODEL if used else GREY, zorder=4, ls="-" if used else "--")
        ax.set_title(lab, fontsize=6.8)
        ax.set_xlim(-9, 9)
        ax.set_xticks([-6, 0, 6])
        if not used:
            ax.axvline(-3.18, color=INK, lw=0.5, ls=":")
            ax.text(-3.0, 1700, "good data\nbegin", fontsize=5.5, va="top")
    axes[0].set_ylim(-3600, 2200)
    axes[0].set_ylabel("Relative flux (ppm)")
    fig.supxlabel("Hours from mid-transit", fontsize=8, y=-0.04)
    fig.subplots_adjust(wspace=0.06)
    save(fig, "fig1")


def fig2():
    au = json.load(open(RESULTS / "candidates" / "audit_402710280.json"))
    import io
    import s3io
    from astropy.io import fits
    d = fits.open(io.BytesIO(s3io.fetch(s3io.key_for("tspoc", 402710280, 66))))[1].data
    t = np.asarray(d["TIME"], float)
    w = np.isfinite(t) & (np.abs(t - 3111.40) < 0.7)
    tc = 3111.40
    h = (t[w] - tc) * 24
    q = np.asarray(d["QUALITY"])[w]
    sap, pdc, bkg = (np.asarray(d[c], float)[w] for c in ("SAP_FLUX", "PDCSAP_FLUX", "SAP_BKG"))
    q0 = q == 0
    ref_s, ref_p = np.nanmedian(sap[q0 & (h > 10)]), np.nanmedian(pdc[q0 & (h > 10)])
    fig, ax = plt.subplots(3, 1, figsize=(COL1, 3.6), sharex=True)
    ax[0].plot(h[~q0], (sap[~q0] / ref_s - 1) * 1e6, ".", ms=1, color="#c9c9c9", label="flagged")
    ax[0].plot(h[q0], (sap[q0] / ref_s - 1) * 1e6, ".", ms=1, color=INK, label="quality = 0")
    ax[0].set_ylabel("SAP (ppm)")
    ax[0].legend(loc="lower left", frameon=False, markerscale=4)
    ax[1].plot(h[q0], (pdc[q0] / ref_p - 1) * 1e6, ".", ms=1, color=INK)
    ax[1].set_ylabel("PDCSAP (ppm)")
    ax[2].plot(h, bkg / 1e3, "-", lw=0.8, color=MODEL)
    ax[2].set_ylabel(r"Background ($10^3$ e$^-$ s$^{-1}$)")
    ax[2].set_xlabel("Hours from predicted mid-transit (BTJD 3111.40)")
    for a in ax:
        a.axvspan(-2.64, 2.64, color=DATA, alpha=0.10, lw=0)
    save(fig, "fig2")
    return au


def fig3(fit, ex):
    f = fit["TIC402710280"]
    P, T, ref = f["period_d"][0], f["tref_btjd"][0], f["ref_epoch"]
    eP, eT = f["period_d"][1], f["tref_btjd"][1]
    rms = math.sqrt(np.mean([t["o_c_min"] ** 2 for t in f["times"]]))
    fig, ax = plt.subplots(figsize=(COL1, 2.2))
    n = np.linspace(-5, 62, 300)
    band = np.sqrt(eT ** 2 + ((n - ref) * eP) ** 2) * 1440
    ax.fill_between(n, -band, band, color=DATA, alpha=0.18, lw=0, label="ephemeris, 1σ")
    tot = np.sqrt(band ** 2 + rms ** 2)
    ax.plot(n, tot, ":", color=DATA, lw=0.7, label="incl. observed scatter")
    ax.plot(n, -tot, ":", color=DATA, lw=0.7)
    for t in f["times"]:
        mk = {"tspoc": "s", "qlp": "D", "ffi": "o"}[t["product"]]
        ax.errorbar(t["epoch"], t["o_c_min"], t["err_min"], fmt=mk, ms=3.5, color=INK, lw=0.7, mfc="white" if t["product"] == "ffi" else INK)
        ax.annotate(f"S{t['sector']}", (t["epoch"], t["o_c_min"]), textcoords="offset points", xytext=(4, 3), fontsize=6)
    for e in ex["all_next"]:
        if e["epoch"] in (57, 59):
            ax.axvline(e["epoch"], color=MODEL, lw=0.6, ls="--")
    ax.axhline(0, color=GREY, lw=0.5)
    ax.set_xlim(-5, 62)
    ax.set_ylim(-60, 60)
    ax.set_xlabel("Epoch (cycles from sector 39)")
    ax.set_ylabel("O − C (min)")
    ax.legend(loc="lower left", frameon=False, fontsize=6)
    save(fig, "fig3")


def fig4():
    bl = json.load(open(RESULTS / "candidates" / "blend_shape_402710280.json"))
    fig, ax = plt.subplots(figsize=(COL1, 2.0))
    lds = sorted({tuple(r["u"]) for r in bl["eb"]})
    for u, mk in zip(lds, ("o", "s", "^")):
        rr = sorted([r for r in bl["eb"] if tuple(r["u"]) == u], key=lambda r: r["k"])
        ax.plot([r["k"] for r in rr], [r["dchi2"] for r in rr], "-" + mk, ms=3, lw=0.7, color=INK,
                mfc="white", label=f"$u$ = ({u[0]:.2f}, {u[1]:.2f})")
    ax.axhline(0, color=GREY, lw=0.5)
    ax.set_yscale("log")
    ax.set_ylim(80, 3000)
    ax.set_xlabel(r"Radius ratio of the eclipsing companion, $k$")
    ax.set_ylabel(r"$\chi^2_{\rm EB} - \chi^2_{\rm planet}$")
    ax.legend(frameon=False, fontsize=6, title="background-star limb darkening", title_fontsize=6)
    save(fig, "fig4")


def main():
    from fit_new import lcs_for
    fit = json.load(open(RESULTS / "report" / "transit_fits.json"))
    ex = json.load(open(RESULTS / "report" / "extras.json"))["TIC402710280"]
    lcs = [lc for lc in lcs_for(402710280) if lc["time"].size > 100]
    fig1(fit, lcs)
    fig2()
    fig3(fit, ex)
    fig4()
    print("figures written to", OUT)


if __name__ == "__main__":
    main()
