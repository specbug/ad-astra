"""Transit fits of the headline candidates that are not new discoveries (results/report/transit_fits.json keys
HIP71639, HIP56966, HD11506, TOI2493). Until 15 October 2026 these were run from scratch scripts outside the repository
(found while answering the independent audits); this file reproduces them, now with emcee convergence diagnostics
(transit_fit.convergence) and 30,000-step chains (the b-rho-k degeneracy gives autocorrelation times
near 450 steps, so 8,000 steps was only ~18 of them). TIC 402710280 is fitted in fit_new.py.

  - HIP 71639 (TIC 128703021): SPOC transit times from its DV epoch at P = 740.4158/31 d in sectors 11, 38, 65,
    and our sector-102 transit; periodic fit.
  - HIP 56966 (TIC 266376971): sectors 10, 99, 100 (2-minute data); the sector-64 transit is full-frame only and is
    used for timing and depth elsewhere; periodic fit.
  - HD 11506 (TIC 266680498): two transits (sectors 30 and 107), period unknown: shape-only fit.
  - TOI-2493.01 (TIC 123664207): three transits in SPOC 2-min data (sectors 33, 98, 107); periodic fit, with its
    derived quantities and upcoming transits written to results/report/extras.json.

  python fit_headline.py [nstep] [HIP71639 HIP56966 HD11506 TOI2493]
"""
import json
import sys

import env  # noqa: F401

import transit_fit as tf
from env import RESULTS

FITS = RESULTS / "report" / "transit_fits.json"
EXTRAS = RESULTS / "report" / "extras.json"


def save(key, res):
    """Merge one fit into transit_fits.json, re-reading it first (several fits may run in parallel)."""
    fits = json.load(open(FITS))
    fits[key] = {**fits.get(key, {}), **res}
    json.dump(fits, open(FITS, "w"), indent=1, default=float)


def main(nstep=30000, which=("HIP71639", "HIP56966", "HD11506", "TOI2493")):
    if "HIP71639" in which:
        P = 740.4157984815635 / 31
        t0 = 1601.8441280984866
        ev = [(11, t0), (38, t0 + 31 * P), (65, t0 + 62 * P), (102, 4133.5817)]
        res, _ = tf.fit_periodic(128703021, ev, [0, 31, 62, 106], P, 7.3, 6281, nstep=nstep)
        save("HIP71639", res)
        print("HIP71639", res["chain_summary"], flush=True)
    if "HIP56966" in which:
        res, _ = tf.fit_periodic(266376971, [(10, 1586.767), (99, 4065.044), (100, 4084.870)], [0, 125, 126], 19.82621,
                                 3.0, 5967, nstep=nstep, win=0.6)
        save("HIP56966", res)
        print("HIP56966", res["chain_summary"], flush=True)
    if "HD11506" in which:
        res, _ = tf.fit_shape(266680498, [(30, 2121.6914), (107, 4271.0803)], 3.3, 6030, nstep=nstep)
        save("HD11506", res)
        print("HD11506", res["chain_summary"], flush=True)
    if "TOI2493" in which:
        import predictions
        from lcaccess import lcs_for
        data = tf.windows_from_lcs(lcs_for(123664207), [2225.82, 4043.534, 4280.627], 0.8, 5.0)
        res, _ = tf.fit_periodic_windows(data, [0, 23, 26], 79.03019, 4242.0, nstep=nstep)
        res["sectors"] = [33, 98, 107]
        save("TOI2493", res)
        ex = json.load(open(EXTRAS))
        ex["TOI2493"] = predictions.periodic(json.load(open(FITS)), "TOI2493", 0.66, 0.652, 4242.0, 89.74609, -21.57558,
                                             now=4325.5)
        json.dump(ex, open(EXTRAS, "w"), indent=1, default=float)
        print("TOI2493", res["chain_summary"], flush=True)


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 30000, tuple(sys.argv[2:]) or ("HIP71639", "HIP56966", "HD11506", "TOI2493"))
