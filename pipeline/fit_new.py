"""Transit fits of the two candidates found by the period hunt, for the report (results/report/transit_fits.json,
keys TIC402710280 and TOI2493), and their derived quantities and upcoming transits (results/report/extras.json).

TIC 402710280.01: five transits (S39 TESS-SPOC, S66 TESS-SPOC, S93 QLP, S100 and S102 from our own full-frame
photometry, ffi_new.py). The S102 transit lies 1.4 d from the start of an orbit, where the light curve has a
ramp that the window's quadratic baseline does not follow; it alone pulls the fitted stellar density to 0.19
solar (all five: 0.19; without it: 0.78, and the K dwarf's 1.61 costs dchi2 = 0.2). The shape is therefore fitted
on the other four; the S102 window is shown with that model.
TOI-2493.01: three transits in SPOC 2-min data (S33, S98, S107).

  python fit_new.py
"""
import json

import env  # noqa: F401
import numpy as np

import report_extras as rx
import ticdb
import transit_fit as tf
from dossier import all_sectors, best_lcs
from env import RESULTS
from ffi_new import new_lcs

FITS = RESULTS / "report" / "transit_fits.json"
EXTRAS = RESULTS / "report" / "extras.json"


def lcs_for(tic):
    st = ticdb.query([tic], ["ra", "dec"]).iloc[0]
    lcs = best_lcs(tic, all_sectors(tic))
    have = {lc["meta"]["sector"] for lc in lcs}
    return lcs + [lc for lc in new_lcs(tic, float(st.ra), float(st.dec), build_missing=False) if lc["meta"]["sector"] not in have]


def main():
    fits = json.load(open(FITS))
    ex = json.load(open(EXTRAS))
    # TIC 402710280.01
    tic = 402710280
    lcs = lcs_for(tic)
    tcs = [2384.607, 3111.411, 3838.214, 4092.595, 4128.936]
    epochs = [0, 20, 40, 47, 48]
    data = tf.windows_from_lcs(lcs, tcs, 0.8, 5.3)
    out, ch = tf.fit_periodic_windows(data[:4], epochs[:4], 36.34001, 5182.0, nstep=6000)
    med = np.median(ch, axis=0)
    P, tref = med[1], med[0]
    aR = ((10 ** med[4]) * tf.RHO_SUN * tf.G * (P * 86400) ** 2 / (3 * np.pi)) ** (1 / 3)
    inc = np.degrees(np.arccos(med[3] / aR))
    out["plot"] = tf.plot_data([(w[0], w[1], w[2]) for w in data], [tref + (ep - out["ref_epoch"]) * P for ep in epochs],
                               P, med[2], aR, inc, out["u"])
    out["sectors"] = [39, 66, 93, 100, 102]
    out["fit_windows"] = [39, 66, 93, 100]
    fits["TIC402710280"] = out
    print("TIC 402710280", {k: v for k, v in out.items() if k != "plot"}, flush=True)
    # TOI-2493.01
    tic = 123664207
    lcs = lcs_for(tic)
    data = tf.windows_from_lcs(lcs, [2225.82, 4043.534, 4280.627], 0.8, 5.0)
    out, ch = tf.fit_periodic_windows(data, [0, 23, 26], 79.03019, 4242.0, nstep=6000)
    out["sectors"] = [33, 98, 107]
    fits["TOI2493"] = out
    print("TOI-2493", {k: v for k, v in out.items() if k != "plot"}, flush=True)
    json.dump(fits, open(FITS, "w"), indent=1, default=float)
    ex["TIC402710280"] = rx.periodic(fits, "TIC402710280", 0.88, 0.817, 5182.0, 248.524, -63.55, now=4325.5)
    ex["TOI2493"] = rx.periodic(fits, "TOI2493", 0.66, 0.652, 4242.0, 89.74609, -21.57558, now=4325.5)
    json.dump(ex, open(EXTRAS, "w"), indent=1, default=float)
    for k in ("TIC402710280", "TOI2493"):
        print(k, {a: b for a, b in ex[k].items() if a not in ("next", "all_next")})
        for x in ex[k]["next"]:
            print("   ", x["utc"], x["note"])


if __name__ == "__main__":
    main()
