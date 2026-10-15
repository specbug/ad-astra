"""Derived quantities for the report's headline candidates (results/report/extras.json).

HIP 71639: orbit, insolation, expected RV amplitude, upcoming transits with night-time coverage from La Silla.
HD 11506: period aliases of the duo (duo.alias_table on SPOC sectors 3, 30, 107) and, for each allowed period,
the impact parameter a circular orbit would need, and the minimum eccentricity at the fitted impact parameter.
"""
import json
import math

import env  # noqa: F401

import s3io
from duo import alias_table, flatten_sectors
from env import RESULTS
from predictions import a_over_rs, periodic

def hip71639(fit):
    def dv(f, T, P, ref):
        # S102 transit vs the SPOC DV ephemeris (s0001-s0065: epoch 1601.84413, period 740.41580 = 31 P)
        p_dv = 740.4157984815635 / 31
        return {"s102_offset_min": (T + (106 - ref) * P - (1601.8441280984866 + 106 * p_dv)) * 1440,
                "alias_min_sigma": 38.8}
    return periodic(fit, "HIP71639", 1.23, 1.596, 6281.0, 219.795192626349, -46.5844751296518, dv)


def hip56966(fit):
    return periodic(fit, "HIP56966", 1.09, 1.069, 5967.0, 175.169514349685, -58.1043454687562,
                    lambda f, T, P, ref: {"alias_min_sigma": 16.4})


def hd11506(fit):
    f = fit["HD11506"]
    t1, t2 = f["tc_0"][0], f["tc_1"][0]
    dur = f["t14_h"][0]
    depth = f["depth_ppm"][0] * 1e-6
    tic = 266680498
    lcs = [lc for lc in (s3io.load("spoc", tic, s) for s in (3, 30, 107)) if lc is not None]
    flat = flatten_sectors(lcs)
    rows = alias_table([t1, t2], depth, dur, flat, p_min=5.0)
    rho = 1.12 / 1.35 ** 3
    k = f["rp_rs"][0]
    b50 = f["b"][0]
    for a in rows:
        P = a["P"]
        aR = a_over_rs(P, rho)
        x = math.pi * (dur / 24) * aR / P
        a["b_circ"] = math.sqrt(max((1 + k) ** 2 - x ** 2, 0.0))
        tcirc = P / (math.pi * aR) * math.sqrt(max((1 + k) ** 2 - b50 ** 2, 1e-9)) * 24
        rr = min(dur / tcirc, 1.0)
        a["e_min_b50"] = (1 - rr ** 2) / (1 + rr ** 2)
    per = f["per_transit"]
    d = abs(per[0]["rp_rs"] - per[1]["rp_rs"])
    s = math.hypot((per[0]["rp_rs_hi"] - per[0]["rp_rs_lo"]) / 2, (per[1]["rp_rs_hi"] - per[1]["rp_rs_lo"]) / 2)
    return {"dt": t2 - t1, "aliases": rows, "depth_diff_sigma": d / s}


def main():
    fit = json.load(open(RESULTS / "report" / "transit_fits.json"))
    p = RESULTS / "report" / "extras.json"
    ex = json.load(open(p)) if p.exists() else {}          # keep the entries fit_new.py writes
    ex.update({"HIP71639": hip71639(fit), "HD11506": hd11506(fit)})
    if "HIP56966" in fit:
        ex["HIP56966"] = hip56966(fit)
    json.dump(ex, open(RESULTS / "report" / "extras.json", "w"), indent=1, default=float)
    for key in ("HIP71639", "HIP56966"):
        h = ex.get(key)
        if not h:
            continue
        print(key, {k: v for k, v in h.items() if k not in ("next", "all_next")})
        for x in h["next"]:
            print("   ", x["epoch"], x["utc"], f"{x['night_frac']:.2f}", f"{x['sun_sep']:.0f}")
    al = ex["HD11506"]["aliases"]
    ok = [a for a in al if a["status"] == "allowed"]
    print("HD 11506 aliases:", len(al), "allowed:", len(ok), "shortest:", ok[-1]["P"] if ok else None)
    for a in ok[::-1][:8]:
        print(f"  n={a['n']} P={a['P']:.3f} b_circ={a['b_circ']:.2f} e_min={a['e_min_b50']:.2f}")


if __name__ == "__main__":
    main()
