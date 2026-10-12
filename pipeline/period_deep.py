"""Transit fits for the new periods of the period review (the stars period_vet.py tested at the pixel level).

For each candidate and each period still allowed (at most three), the transits are fitted jointly
(transit_fit.fit_periodic_windows: one shape, free stellar density, the period held near the alias) and the
density the shape implies is compared with the TIC density of the star. A period at which the implied density
is far from the star's needs an eccentric orbit or a different star; with several aliases, the one that matches
the star is the likelier period.

  python period_deep.py
Output: results/candidates/period_deep.json, keyed "tic|id"
"""
import json
import math

import env  # noqa: F401
import numpy as np

import ticdb
import transit_fit as tf
from env import RESULTS
from fit_new import lcs_for
from period_vet import new_periods

C = RESULTS / "candidates"
OUT = C / "period_deep.json"


def fit_one(o, P, lcs, teff):
    times = [o["ref_epoch"]] + [e["t"] for e in o["events"] if e.get("ok")]
    times = sorted(times)
    ep = [int(round((t - times[0]) / P)) for t in times]
    data = tf.windows_from_lcs(lcs, times, max(0.6, 2.5 * o["dur_h"] / 24), o["dur_h"])
    if len(data) < 2:
        return None
    # keep only windows we got, with their epochs
    got = [w[4] for w in data]
    ep = [e for t, e in zip(times, ep) if any(abs(t - g) < 1e-6 for g in got)]
    out, ch = tf.fit_periodic_windows(data, ep, P, teff, nstep=2500)
    out.pop("plot", None)
    return out


def main():
    rv = json.load(open(C / "period_review.json"))
    pv = json.load(open(C / "period_vet.json")) if (C / "period_vet.json").exists() else {}
    res = json.load(open(OUT)) if OUT.exists() else {}
    for o in new_periods(rv):
        key = f"{o['tic']}|{o['id']}"
        if key in res:
            continue
        st = ticdb.query([int(o["tic"])], ["Teff", "rho", "rad", "mass", "Tmag"]).iloc[0]
        teff = float(st.Teff) if st.Teff == st.Teff else 5500.0
        rho = float(st.rho) if st.rho == st.rho else None
        lcs = lcs_for(int(o["tic"]))
        fits = []
        for P in (o.get("aliases") or [])[:3]:
            try:
                f = fit_one(o, P, lcs, teff)
            except Exception as ex:
                f = {"error": repr(ex)[:150]}
            if f and "rho_solar" in f and rho:
                r = f["rho_solar"]
                z = (r[0] - rho) / (r[2] if r[0] > rho else r[1])
                f["rho_tic"], f["rho_z"] = rho, float(z)
            fits.append({"P": P, **(f or {})})
        res[key] = {"tic": int(o["tic"]), "id": o["id"], "tmag": float(st.Tmag), "rad": float(st.rad) if st.rad == st.rad else None,
                    "fits": fits, "pixels": pv.get(key, {}).get("verdict")}
        print(key, [(round(f["P"], 3), round(f.get("rho_solar", [np.nan])[0], 2), round(f.get("rho_z", np.nan), 1)) for f in fits], flush=True)
        json.dump(res, open(OUT, "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
