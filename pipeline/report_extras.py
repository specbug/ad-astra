"""Derived quantities for the report's headline candidates (results/report/extras.json).

HIP 71639: orbit, insolation, expected RV amplitude, upcoming transits with night-time coverage from La Silla.
HD 11506: period aliases of the duo (duo.alias_table on SPOC sectors 3, 30, 107) and, for each allowed period,
the impact parameter a circular orbit would need, and the minimum eccentricity at the fitted impact parameter.
"""
import json
import math

import env  # noqa: F401
import numpy as np

import s3io
from duo import alias_table, flatten_sectors
from env import RESULTS

G = 6.674e-11
RHO_SUN = 1408.0
AU_RSUN = 215.032


def a_over_rs(P_d, rho_solar):
    return (G * rho_solar * RHO_SUN * (P_d * 86400) ** 2 / (3 * math.pi)) ** (1 / 3)


def observability(tc_btjd, t14_h, ra, dec, site):
    """Fraction of the transit (+-T14/2) when the Sun is below -12 deg and the star above 30 deg."""
    from astropy.coordinates import AltAz, EarthLocation, SkyCoord, get_sun
    from astropy.time import Time
    import astropy.units as u
    loc = EarthLocation(lat=site[0] * u.deg, lon=site[1] * u.deg, height=site[2] * u.m)
    ts = Time(tc_btjd + 2457000.0 + np.linspace(-t14_h / 48, t14_h / 48, 25), format="jd", scale="tdb")
    fr = AltAz(obstime=ts, location=loc)
    star = SkyCoord(ra * u.deg, dec * u.deg).transform_to(fr)
    sun = get_sun(ts).transform_to(fr)
    ok = (sun.alt.deg < -12) & (star.alt.deg > 30)
    sun_sep = float(get_sun(ts[12]).separation(SkyCoord(ra * u.deg, dec * u.deg)).deg)
    return float(ok.mean()), sun_sep


def hip71639(fit):
    f = fit["HIP71639"]
    P, Pe = f["period_d"][0], max(f["period_d"][1:])
    T, Te = f["tref_btjd"][0], max(f["tref_btjd"][1:])
    ref = f["ref_epoch"]
    M, R, Teff = 1.23, 1.596, 6281.0
    a_au = (M * (P / 365.25) ** 2) ** (1 / 3)
    L = R ** 2 * (Teff / 5772.0) ** 4
    teq = Teff * math.sqrt(R / AU_RSUN / (2 * a_au))
    rp = f["rp_rs"][0] * R * 109.08
    mass = (rp / 0.808) ** (1 / 0.589)  # Chen & Kipping (2017) Neptunian branch
    k = 28.43 * (mass / 317.83) * M ** (-2 / 3) * (P / 365.25) ** (-1 / 3)
    # S102 transit vs the SPOC DV ephemeris (s0001-s0065: epoch 1601.84413, period 740.41580 = 31 P)
    p_dv = 740.4157984815635 / 31
    s102_pred = 1601.8441280984866 + 106 * p_dv
    s102_obs = T + (106 - ref) * P
    nxt = []
    now = 4323.0  # 2026-10-10
    n0 = int(math.ceil((now - T) / P)) + ref
    for n in range(n0, n0 + 30):
        tc = T + (n - ref) * P
        sig = math.hypot(Te, (n - ref) * Pe) * 1440
        frac, sun_sep = observability(tc, f["t14_h"][0], 219.795192626349, -46.5844751296518, (-29.2563, -70.7380, 2400))
        from report import btjd_to_utc
        note = ("full transit at night" if frac > 0.95 else f"{frac:.0%} of the transit at night" if frac > 0.05
                else "daytime or too low" if sun_sep > 50 else "too close to the Sun")
        nxt.append({"epoch": n, "btjd": tc, "utc": btjd_to_utc(tc).strftime("%Y-%m-%d %H:%M"), "sigma_min": sig,
                    "night_frac": frac, "sun_sep": sun_sep, "note": note})
    keep = [x for x in nxt if x["night_frac"] > 0.05][:8]
    rho_tic = M / R ** 3
    r = f["rho_solar"]
    z = (r[0] - rho_tic) / (r[2] if r[0] > rho_tic else r[1])
    return {"a_au": a_au, "s_earth": L / a_au ** 2, "teq": teq, "mass_assumed": mass, "k_ms": k,
            "s102_offset_min": (s102_obs - s102_pred) * 1440, "alias_min_sigma": 38.8, "rho_tic": rho_tic,
            "rho_verdict": "consistent" if abs(z) < 2 else "in tension", "rho_z": z, "next": keep, "all_next": nxt}


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
    ex = {"HIP71639": hip71639(fit), "HD11506": hd11506(fit)}
    json.dump(ex, open(RESULTS / "report" / "extras.json", "w"), indent=1, default=float)
    h = ex["HIP71639"]
    print({k: v for k, v in h.items() if k not in ("next", "all_next")})
    for x in h["all_next"][:30]:
        print(x["epoch"], x["utc"], f"{x['night_frac']:.2f}", f"{x['sun_sep']:.0f}")
    al = ex["HD11506"]["aliases"]
    ok = [a for a in al if a["status"] == "allowed"]
    print("HD 11506 aliases:", len(al), "allowed:", len(ok), "shortest:", ok[-1]["P"] if ok else None)
    for a in ok[::-1][:8]:
        print(f"  n={a['n']} P={a['P']:.3f} b_circ={a['b_circ']:.2f} e_min={a['e_min_b50']:.2f}")


if __name__ == "__main__":
    main()
