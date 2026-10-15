"""Diagnostic figures for individual events and stars."""
import numpy as np

import env  # noqa: F401
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import s3io
from lcprep import bin_lc, detrend_fast, robust_std


def event_figure(product, tic, sector, t0, dur_h, out, title=None, lc=None):
    """Full sector light curve, zoom on the event, background and centroids around it."""
    lc = lc or s3io.load(product, tic, sector)
    if lc is None:
        return None
    t, f = lc["time"], lc["flux"]
    o = np.argsort(t)
    t, f = t[o], f[o]
    tb, fb, _, ex = bin_lc(t, f, 10.0, {k: lc[k][o] for k in ("bkg", "cx", "cy", "flux_sml", "flux_lag") if k in lc})
    yl, trend = detrend_fast(tb, fb, 2.5, knot_step=2.5 / 36)
    D = dur_h / 24
    fig, ax = plt.subplots(5, 1, figsize=(10, 12), gridspec_kw={"height_ratios": [1.2, 1.6, 0.8, 0.8, 0.8]})
    ax[0].plot(t, f, ",", color="0.7")
    ax[0].plot(tb, fb, ".", ms=1.5, color="k")
    ax[0].plot(tb, trend, "-", color="tab:orange", lw=0.8)
    ax[0].axvline(t0, color="r", alpha=0.4)
    ax[0].set_ylabel("flux (raw, binned 10m)")
    w = max(4 * D, 1.0)
    m = np.abs(t - t0) < w
    mb = np.abs(tb - t0) < w
    ax[1].plot(t[m] - t0, f[m] / np.interp(t[m], tb, trend), ".", ms=2, color="0.6")
    ax[1].plot(tb[mb] - t0, yl[mb], "o", ms=3, color="k")
    for k, c in (("flux_sml", "tab:blue"), ("flux_lag", "tab:green")):
        if k in ex:
            v = ex[k][mb] / np.nanmedian(ex[k][mb])
            ax[1].plot(tb[mb] - t0, v - 0.0, "-", lw=0.8, color=c, alpha=0.7, label=k)
    ax[1].axvspan(-D / 2, D / 2, color="r", alpha=0.08)
    ax[1].set_ylabel("flattened flux")
    ax[1].legend(fontsize=7)
    for i, (k, lab) in enumerate((("bkg", "background"), ("cx", "centroid x"), ("cy", "centroid y"))):
        a = ax[2 + i]
        if k in ex:
            a.plot(tb[mb] - t0, ex[k][mb], ".", ms=2, color="k")
            a.axvspan(-D / 2, D / 2, color="r", alpha=0.08)
        a.set_ylabel(lab)
    ax[-1].set_xlabel("time from event (d)")
    meta = lc["meta"]
    fig.suptitle(title or f"TIC {tic}  S{sector} {product}  cam{meta['camera']} ccd{meta['ccd']}  "
                 f"T={meta['tmag']:.2f}  t0={t0:.3f}  D={dur_h:g}h", fontsize=10)
    fig.tight_layout()
    fig.savefig(out, dpi=90)
    plt.close(fig)
    return out
