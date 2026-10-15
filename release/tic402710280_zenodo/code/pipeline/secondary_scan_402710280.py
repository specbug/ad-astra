"""Secondary-eclipse search for TIC 402710280.01 over orbital phases 0.3-0.7 (eccentric orbits with |e cos w| up to
about 0.3 put the secondary there), in all eight sectors.

Each light curve has the 36.34-d transits and the 3.78-d signal's transits masked and is divided by a 1.5-d running
median; a box of the transit duration (5.27 h) is stepped across phase in 0.001 steps, and at each phase the mean
depth of all covered cycles is measured, with an error from the out-of-box scatter scaled by the red-noise factor of
the binned light curve.

  python secondary_scan_402710280.py
Output: results/candidates/secondary_scan_402710280.json
"""
import json

import env  # noqa: F401
import numpy as np

from env import RESULTS
from fit_new import lcs_for
from second_signal_402710280 import phase_dist, running_median

P, T0 = 36.33989590709609, 2460838.20190249 - 2457000
P2, T2 = 3.7837059787718332, 4095.308738819944
DUR = 5.27 / 24
OUT = RESULTS / "candidates" / "secondary_scan_402710280.json"


def main():
    secs = []
    for lc in lcs_for(402710280):
        t, f = np.asarray(lc["time"], float), np.asarray(lc["flux"], float)
        ok = np.isfinite(t) & np.isfinite(f)
        t, f = t[ok], f[ok]
        keep = (phase_dist(t, P, T0) > 0.3) & (phase_dist(t, P2, T2) > 0.05)
        t, f = t[keep], f[keep]
        tr = running_median(t, f, np.zeros(t.size, bool), w=1.5)
        ok = np.isfinite(tr)
        t, r = t[ok], f[ok] / tr[ok] - 1
        sig = 1.4826 * np.median(np.abs(np.diff(r))) / np.sqrt(2)          # white noise per cadence
        idx = np.floor((t - t.min()) / DUR)
        bins = [(r[idx == k].mean(), (idx == k).sum()) for k in np.unique(idx) if (idx == k).sum() > 10]
        bm, bn = np.array(bins).T
        beta = float(max(1.0, 1.4826 * np.median(np.abs(bm - np.median(bm))) / (sig / np.sqrt(np.median(bn)))))
        secs.append((int(lc["meta"]["sector"]), t, r, sig, beta))
    rows = []
    for p0 in np.arange(0.3, 0.7001, 0.001):
        ds, ws = [], []
        for s_, t, r, sig, beta in secs:
            ph = ((t - T0) / P) % 1
            inn = np.abs(((ph - p0 + 0.5) % 1) - 0.5) * P < DUR / 2
            if inn.sum() < 10:
                continue
            ds.append(-r[inn].mean() * 1e6)
            ws.append(1 / (sig * beta / np.sqrt(inn.sum()) * 1e6) ** 2)
        if ds:
            ds, ws = np.array(ds), np.array(ws)
            rows.append((float(p0), float(np.sum(ds * ws) / ws.sum()), float(1 / np.sqrt(ws.sum())), len(ds)))
    rows = np.array(rows)
    z = rows[:, 1] / rows[:, 2]
    i = int(np.argmax(z))
    out = {"beta": {s_: b for s_, _, _, _, b in secs}, "n_phases": len(rows),
           "fraction_of_0.3_0.7_covered": float(len(rows) / 401),
           "deepest": {"phase": rows[i, 0], "depth_ppm": rows[i, 1], "err_ppm": rows[i, 2], "sigma": float(z[i])},
           "median_3sigma_limit_ppm": float(np.median(rows[:, 1] + 3 * rows[:, 2])),
           "max_3sigma_limit_ppm": float(np.max(rows[:, 1] + 3 * rows[:, 2])),
           "rows": rows.tolist()}
    json.dump(out, open(OUT, "w"), indent=1)
    print({k: v for k, v in out.items() if k != "rows"})


if __name__ == "__main__":
    main()
