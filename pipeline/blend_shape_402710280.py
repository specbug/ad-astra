"""Could the 3.1-3.6" background star (Gaia DR3 5828441289812582656, G = 16.15) make TIC 402710280.01's transits?

Diluted by the target (Delta G = 6.37, a flux fraction of about 0.28%), a ~2100 ppm dip needs an eclipse of about
74% of that star's light: a companion nearly as large as the star, and a transit whose ingress and egress last a
large fraction of its duration. The four clean transits (sectors 39, 93, 100, 102; each centred on its measured
mid-time, quadratic baseline divided out, errors scaled by the red-noise factor) are fitted with
  - a planet: radius ratio k, impact parameter b and a/R* free, undiluted;
  - an eclipse of the background star: k fixed on a grid (0.3-10), b and a/R* free, diluted by the fraction F
    (free, but at most the star's flux share, generously 0.4% to allow for the TESS-band colour difference),
and the chi-square difference is reported for each k. A smaller k cannot reach the depth within the flux budget,
which the fit enforces. For k > 1 the eclipsed star would have to supply most of the neighbour's light while being
smaller than its companion, which no main-sequence pair does; those k are fitted to show how the shape behaves.

  python blend_shape_402710280.py
Output: results/candidates/blend_shape_402710280.json
"""
import json
import math

import env  # noqa: F401
import numpy as np
from scipy.optimize import minimize

import transit_fit as tf
from audit_402710280 import TRANSITS, beta_factor
from env import RESULTS

OUT = RESULTS / "candidates" / "blend_shape_402710280.json"
P = 36.3399
F_MAX = 0.004
U_BG = (0.40, 0.25)        # limb darkening for a background star of unknown type; varied below


def data():
    from fit_new import lcs_for
    au = json.load(open(RESULTS / "candidates" / "audit_402710280.json"))
    tcs = {x["sector"]: x["tc"] for x in au["times_noS102_shape"]}
    lcs = [lc for lc in lcs_for(402710280) if lc["time"].size > 100]
    out = []
    for (t, f, sig, exp_s, tc) in tf.windows_from_lcs(lcs, [tcs[s] for s in sorted(TRANSITS)], 0.8, 5.2):
        oot = np.abs(t - tc) > 0.75 * 5.2 / 24
        out.append((t - tc, f, sig * beta_factor(t[oot], f[oot] - 1, sig), exp_s))
    return out


def chi2(windows, k, b, aR, F, u):
    if not (0 <= b < 1 + k and aR > 1.5 and 0 < F <= 1):
        return 1e12
    inc = math.degrees(math.acos(min(b / aR, 1)))
    c = 0.0
    for (t, f, s, exp_s) in windows:
        m = tf._model(t, 0.0, P, k, aR, inc, u, exp_s=exp_s)
        c += float(np.sum(((f - (1 - F * (1 - m))) / s) ** 2))
    return c


def best(fun, starts):
    r = None
    for x0 in starts:
        q = minimize(fun, x0, method="Nelder-Mead", options={"maxiter": 6000, "xatol": 1e-6, "fatol": 1e-3})
        if r is None or q.fun < r.fun:
            r = q
    return r


def main():
    W = data()
    u_t = tf.tess_ld(5180.0)
    pl = best(lambda p: chi2(W, p[0], p[1], p[2], 1.0, u_t),
              [[0.042, 0.3, 54], [0.045, 0.6, 45], [0.05, 0.85, 35]])
    out = {"planet": {"chi2": float(pl.fun), "k": float(pl.x[0]), "b": float(pl.x[1]), "a_rs": float(pl.x[2])},
           "n_points": int(sum(len(w[0]) for w in W)), "eb": []}
    print("planet chi2 %.1f  k %.4f b %.2f a/R %.1f" % (pl.fun, *pl.x), flush=True)
    for u_bg in (U_BG, (0.25, 0.30), (0.55, 0.15)):
        for k in (0.3, 0.5, 0.7, 0.85, 1.0, 1.2, 1.5, 2.0, 3.0, 5.0, 10.0):
            def fun(p):
                F = F_MAX / (1 + math.exp(-p[2]))
                return chi2(W, k, abs(p[0]), p[1], F, u_bg)
            r = best(fun, [[0.2, 30, 2], [0.6, 20, 2], [0.9 + k * 0.3, 15, 4], [0.3, 60, 0]])
            F = F_MAX / (1 + math.exp(-r.x[2]))
            row = {"k": k, "u": u_bg, "chi2": float(r.fun), "dchi2": float(r.fun - pl.fun), "b": float(abs(r.x[0])),
                   "a_rs": float(r.x[1]), "F": F, "eclipse_depth_needed": None}
            out["eb"].append(row)
            print("EB k=%.2f u=%s chi2 %.1f dchi2 %+.1f b %.2f a/R %.1f F %.4f" % (k, u_bg, r.fun, r.fun - pl.fun,
                                                                                   abs(r.x[0]), r.x[1], F), flush=True)
    json.dump(out, open(OUT, "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
