"""Transit fits for the corrected ephemerides (posted_cases.py) and for the separated community TOIs.

  - CTOI 438225366.01: all transits at the corrected period (20.86 d, half the posted 41.72 d);
  - TOI-2633.01: the transits on the corrected ephemeris (through the community TOI's transit, 18.215 d);
  - CTOI 156716001.01, 369532319.01, 67711670.01: the community TOI's own transits, which are not transits of
    the TOI they were merged into. Their periods are only known to be the gap between transits divided by a
    whole number. For a fixed transit shape the stellar density a fit implies grows in proportion to the assumed
    period (rho = 3 pi (a/R)^3 / (G P^2), with a/R in proportion to P), so one fit at a reference alias gives the
    implied density at every alias, and the alias whose implied density matches the star's is the likeliest
    period for a circular orbit.

  python corrections_fit.py
Output: results/candidates/corrections_fit.json
"""
import json

import env  # noqa: F401
import numpy as np

import ticdb
import transit_fit as tf
from env import RESULTS
from fit_new import lcs_for

C = RESULTS / "candidates"
OUT = C / "corrections_fit.json"

# (key, tic, event times [BTJD], duration [h], reference period, gap for the alias list or None)
FITS = [
    ("CTOI 438225366.01", 438225366, [1470.513, 2221.520, 2492.718, 2534.440, 3243.725, 3264.586, 3285.448, 3681.813], 5.0,
     20.86132, None),
    ("TOI-2633.01", 386131764, [1506.771, 2235.367, 2253.582, 3674.344, 3692.559, 3710.774], 4.0, 18.2149, None),
    ("CTOI 156716001.01", 156716001, [1485.406, 2224.520, 3702.770], 8.0, 739.12032 / 8, 739.12032),
    ("CTOI 369532319.01", 369532319, [1755.810, 2831.560], 5.6, 1075.74954 / 20, 1075.74954),
    ("CTOI 67711670.01", 67711670, [1798.985, 3622.800], 6.0, 1823.8042 / 30, 1823.8042),
]


# light-curve strips around a community TOI's transit and the nearby transits of the TOI it was merged into:
# (key, tic, [(t_lo, t_hi), ...], [(time, label, kind), ...]); kind "ctoi", "toi" or "listed" (no transit there)
STRIPS = [
    ("CTOI 156716001.01", 156716001, [(1484.2, 1487.6), (3699.6, 3711.6)],
     [(1485.406, "CTOI", "ctoi"), (1486.403, "TOI-7602.01", "toi"), (3702.770, "CTOI", "ctoi"),
      (3700.428, "TOI-7602.01", "toi"), (3710.726, "TOI-7602.01", "toi")]),
    ("CTOI 369532319.01", 369532319, [(1752.2, 1756.8), (2826.3, 2843.4)],
     [(1753.114, "TOI-5967.01", "toi"), (1755.810, "CTOI", "ctoi"), (2827.165, "TOI-5967.01", "toi"),
      (2831.560, "CTOI", "ctoi"), (2842.508, "TOI-5967.01", "toi")]),
    ("CTOI 67711670.01", 67711670, [(1797.6, 1805.0), (3619.0, 3627.0)],
     [(1798.985, "CTOI", "ctoi"), (1803.572, "TOI-6033.01 (listed)", "listed"), (3622.800, "CTOI", "ctoi"),
      (3625.632, "TOI-6033.01 (listed)", "listed")]),
    ("TOI-2633.01", 386131764, [(2228.6, 2254.4), (3664.0, 3688.6), (3689.0, 3712.5)],
     [(2229.362, "listed", "listed"), (2247.558, "listed", "listed"), (2235.367, "", "ctoi"), (2253.582, "", "ctoi"),
      (3666.849, "listed", "listed"), (3685.046, "listed", "listed"), (3703.242, "listed", "listed"),
      (3674.344, "", "ctoi"), (3692.559, "", "ctoi"), (3710.774, "", "ctoi")]),
]


def strips(lcs, windows):
    """Binned flux (ppm from the window median) in each window, from the light curve that covers most of it."""
    out = []
    for lo, hi in windows:
        best = None
        for lc in lcs:
            t = lc["time"]
            n = int(((t >= lo) & (t <= hi)).sum())
            if n and (best is None or n > best[0]):
                best = (n, lc)
        if best is None:
            continue
        lc = best[1]
        t, f = lc["time"], lc["flux"]
        w = (t >= lo) & (t <= hi) & np.isfinite(f)
        t, f = t[w], f[w] / np.nanmedian(f[w])
        cad = float(np.median(np.diff(t))) if t.size > 1 else 1 / 48
        bin_d = max(1 / 48, 1.01 * cad, (hi - lo) / 400)
        scatter = 1.4826 * float(np.median(np.abs(np.diff(f)))) / np.sqrt(2)   # per point
        edges = np.arange(lo, hi + bin_d, bin_d)
        idx = np.digitize(t, edges)
        xs, ys, es = [], [], []
        for k in np.unique(idx):
            sel = idx == k
            xs.append(float(t[sel].mean()))
            ys.append(float((f[sel].mean() - 1) * 1e6))
            es.append(float(scatter / np.sqrt(sel.sum()) * 1e6))
        out.append({"lo": lo, "hi": hi, "sector": lc["meta"].get("sector"), "product": lc["meta"].get("product"),
                    "x": xs, "y": ys, "e": es})
    return out


def main():
    res = json.load(open(OUT)) if OUT.exists() else {}
    for key, tic, windows, marks in STRIPS:
        if key in res and "strips" in res[key] and res[key].get("strips_v") == 2:
            continue
        lcs = [lc for lc in lcs_for(tic) if lc["time"].size > 100]
        res.setdefault(key, {})["strips"] = strips(lcs, windows)
        res[key]["marks"] = marks
        res[key]["strips_v"] = 2
        json.dump(res, open(OUT, "w"), indent=1, default=float)
    for key, tic, times, dur_h, P0, gap in FITS:
        if key in res and "period_d" in res[key]:
            continue
        st = ticdb.query([tic], ["Teff", "rho", "rad", "mass"]).iloc[0]
        teff = float(st.Teff) if st.Teff == st.Teff else 5500.0
        rho = float(st.rho) if st.rho == st.rho else None
        lcs = [lc for lc in lcs_for(tic) if lc["time"].size > 100]
        data = tf.windows_from_lcs(lcs, times, max(0.6, 2.5 * dur_h / 24), dur_h)
        got = [w[4] for w in data]
        ep = [int(round((g - times[0]) / P0)) for g in got]
        out, ch = tf.fit_periodic_windows(data, ep, P0, teff, nstep=4000)
        out.update(tic=tic, times=got, epochs=ep, rho_tic=rho, rad_tic=float(st.rad) if st.rad == st.rad else None)
        r = out["rho_solar"]
        if gap and rho:
            # implied density at every alias gap/n (proportional to the period), and its distance from the star's
            als = []
            for n in range(1, 200):
                P = gap / n
                f = P / out["period_d"][0]
                rr = [r[0] * f, r[1] * f, r[2] * f]
                z = (rr[0] - rho) / (rr[2] if rr[0] > rho else rr[1])
                als.append({"n": n, "P": P, "rho": rr[0], "z": float(z)})
                if P < 2:
                    break
            out["aliases"] = als
            out["P_match"] = out["period_d"][0] * rho / r[0]
        elif rho:
            out["rho_z"] = float((r[0] - rho) / (r[2] if r[0] > rho else r[1]))
        res[key] = {**res.get(key, {}), **out}
        print(key, "P=%.5f" % out["period_d"][0], "T0=%.4f" % out["tref_btjd"][0], "depth=%.0f" % out["depth_ppm"][0],
              "b=%.2f" % out["b"][0], "rho=%.3f vs TIC %s" % (r[0], rho), "P_match=%s" % out.get("P_match"), flush=True)
        json.dump(res, open(OUT, "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
