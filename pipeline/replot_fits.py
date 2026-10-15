"""Regenerate the plotted transit panels of stored fits (transit_fit.plot_data) without refitting.

plot_data used to drop every 20-minute bin with fewer than 3 points, which empties the panels of 10- and 30-minute
full-frame sectors (sector 39 of TIC 402710280, for one). This rebuilds the panels from the same data windows and
the stored median parameters.

  python replot_fits.py
Updates "plot" in results/report/transit_fits.json (TIC402710280) and results/candidates/corrections_fit.json.
"""
import json

import env  # noqa: F401
import numpy as np

import transit_fit as tf
from env import RESULTS
from lcaccess import lcs_for

FITS = RESULTS / "report" / "transit_fits.json"
CORR = RESULTS / "candidates" / "corrections_fit.json"


def panels(lcs, times, win, t14_h, f, epochs):
    data = tf.windows_from_lcs(lcs, times, win, t14_h)
    P, tref, rp, aR = f["period_d"][0], f["tref_btjd"][0], f["rp_rs"][0], f["a_rs"][0]
    inc = float(np.degrees(np.arccos(min(f["b"][0] / aR, 1.0))))
    tcs = [tref + (ep - f["ref_epoch"]) * P for ep in epochs]
    return tf.plot_data([(w[0], w[1], w[2]) for w in data], tcs, P, rp, aR, inc, f["u"])


def main():
    fits = json.load(open(FITS))
    f = fits["TIC402710280"]
    f["plot"] = panels(lcs_for(402710280), [2384.607, 3111.411, 3838.214, 4092.595, 4128.936], 0.8, 5.3, f,
                       [0, 20, 40, 47, 48])
    print("TIC402710280", [len(p["x"]) for p in f["plot"]])
    json.dump(fits, open(FITS, "w"), indent=1, default=float)

    from corrections_fit import FITS as CF
    corr = json.load(open(CORR))
    for key, tic, times, dur_h, P0, gap in CF:
        c = corr.get(key)
        if not c or "period_d" not in c:
            continue
        c["plot"] = panels([lc for lc in lcs_for(tic) if lc["time"].size > 100], c["times"],
                           max(0.6, 2.5 * dur_h / 24), dur_h, c, c["epochs"])
        print(key, [len(p["x"]) for p in c["plot"]])
    json.dump(corr, open(CORR, "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
