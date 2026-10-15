"""TRICERATOPS scenario probabilities for the 3.78-day signal of TIC 402710280 (second_signal_402710280.py).

Same machinery as fpp_triceratops.py (TIC field, the photometric apertures of sectors 39, 93, 100 and 102, no
contrast curve), with the light curve replaced by all eight sectors folded on the 3.78-d ephemeris (36.34-d transits
masked, 0.75-d running median divided out), binned to 10 minutes within 0.3 d of mid-transit.

  python fpp_second_signal.py [n_runs]
Output: results/candidates/fpp_tic402710280_3.78d.json
"""
import json
import sys

import env  # noqa: F401
import numpy as np

import fpp_triceratops as F
from env import RESULTS

SS = RESULTS / "candidates" / "second_signal_402710280.json"


def folded():
    from second_signal_402710280 import load, phase_dist  # noqa: F401
    ss = json.load(open(SS))
    P, T0 = ss["bls_joint"]["P"], ss["bls_joint"]["t0_btjd"]
    _, rows = load(P, T0)
    t = np.concatenate([r[1] for r in rows])
    f = np.concatenate([r[2] for r in rows])
    x = ((t - T0 + P / 2) % P) - P / 2
    m = np.abs(x) < 0.3
    x, f = x[m], f[m]
    edges = np.arange(-0.3, 0.3 + 1e-9, 1 / 144)
    k = np.digitize(x, edges)
    tb = np.array([x[k == i].mean() for i in np.unique(k) if (k == i).sum() >= 5])
    fb = np.array([f[k == i].mean() for i in np.unique(k) if (k == i).sum() >= 5])
    err = float(np.std(fb[np.abs(tb) > 0.06]))
    return P, float(ss["depth_all"][0]) * 1e-6, tb, fb, err


def main(n_runs=5):
    P, depth, t, f, err = folded()
    F.P_ORB, F.DEPTH = P, depth
    F.OUT = RESULTS / "candidates" / "fpp_tic402710280_3.78d.json"
    F.FIG = RESULTS / "report" / "cases" / "TIC402710280_3.78d_triceratops_field.png"
    F.light_curve = lambda: (t, f, err)
    print("P", P, "depth", depth, "bins", len(t), "err ppm", round(err * 1e6), flush=True)
    F.main(n_runs)


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 5)
