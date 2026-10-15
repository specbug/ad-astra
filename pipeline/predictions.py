"""Orbit, insolation, radial-velocity amplitude and upcoming transits of a fitted candidate, with each transit's
time as seen from La Silla (barycentric light-travel time and TDB-UTC removed) and its night-time coverage there.
"""
import math

import env  # noqa: F401
import numpy as np

G = 6.674e-11
RHO_SUN = 1408.0
AU_RSUN = 215.032


def a_over_rs(P_d, rho_solar):
    return (G * rho_solar * RHO_SUN * (P_d * 86400) ** 2 / (3 * math.pi)) ** (1 / 3)


LA_SILLA = (-29.2563, -70.7380, 2400)


def site_time(tc_btjd, ra, dec, site=LA_SILLA):
    """The time a transit at barycentric BJD_TDB = tc_btjd + 2457000 is seen at the site, as an astropy Time: the
    barycentric light-travel time to the site, towards the star, is removed by iteration. (Before 15 October 2026
    the report printed BJD_TDB directly as UTC, about 6 minutes off; found by an independent audit.)"""
    from astropy.coordinates import EarthLocation, SkyCoord
    from astropy.time import Time
    import astropy.units as u
    loc = EarthLocation(lat=site[0] * u.deg, lon=site[1] * u.deg, height=site[2] * u.m)
    star = SkyCoord(ra * u.deg, dec * u.deg)
    bjd = Time(np.atleast_1d(tc_btjd) + 2457000.0, format="jd", scale="tdb", location=loc)
    t = bjd
    for _ in range(3):
        t = Time(bjd.jd - t.light_travel_time(star, "barycentric").to_value("day"), format="jd", scale="tdb",
                 location=loc)
    return t


def observability(tc_btjd, t14_h, ra, dec, site):
    """Fraction of the transit (+-T14/2, as seen at the site) when the Sun is below -12 deg and the star above 30 deg."""
    from astropy.coordinates import AltAz, EarthLocation, SkyCoord, get_sun
    from astropy.time import Time
    import astropy.units as u
    loc = EarthLocation(lat=site[0] * u.deg, lon=site[1] * u.deg, height=site[2] * u.m)
    t0 = site_time(tc_btjd, ra, dec, site)[0]
    ts = Time(t0.jd + np.linspace(-t14_h / 48, t14_h / 48, 25), format="jd", scale="tdb")
    fr = AltAz(obstime=ts, location=loc)
    star = SkyCoord(ra * u.deg, dec * u.deg).transform_to(fr)
    sun = get_sun(ts).transform_to(fr)
    ok = (sun.alt.deg < -12) & (star.alt.deg > 30)
    sun_sep = float(get_sun(ts[12]).separation(SkyCoord(ra * u.deg, dec * u.deg)).deg)
    return float(ok.mean()), sun_sep


def periodic(fit, key, M, R, Teff, ra, dec, dv_check=None, now=4323.0, n_show=8):
    """Orbit, insolation, RV amplitude and upcoming transits seen at night from La Silla."""
    f = fit[key]
    P, Pe = f["period_d"][0], max(f["period_d"][1:])
    T, Te = f["tref_btjd"][0], max(f["tref_btjd"][1:])
    ref = f["ref_epoch"]
    a_au = (M * (P / 365.25) ** 2) ** (1 / 3)
    L = R ** 2 * (Teff / 5772.0) ** 4
    teq = Teff * math.sqrt(R / AU_RSUN / (2 * a_au))
    rp = f["rp_rs"][0] * R * 109.08
    mass = (rp / 0.808) ** (1 / 0.589)  # Chen & Kipping (2017) Neptunian branch
    k = 28.43 * (mass / 317.83) * M ** (-2 / 3) * (P / 365.25) ** (-1 / 3)
    nxt = []
    n0 = int(math.ceil((now - T) / P)) + ref
    for n in range(n0, n0 + 40):
        tc = T + (n - ref) * P
        sig = math.hypot(Te, (n - ref) * Pe) * 1440
        frac, sun_sep = observability(tc, f["t14_h"][0], ra, dec, LA_SILLA)
        site_utc = site_time(tc, ra, dec)[0].utc.datetime
        note = ("full transit at night" if frac > 0.95 else f"{frac:.0%} of the transit at night" if frac > 0.05
                else "daytime or too low" if sun_sep > 50 else "too close to the Sun")
        nxt.append({"epoch": n, "btjd": tc, "utc": site_utc.strftime("%Y-%m-%d %H:%M"),
                    "utc_note": "UTC at La Silla (light-travel time and TDB-UTC removed)", "sigma_min": sig,
                    "night_frac": frac, "sun_sep": sun_sep, "note": note})
    rho_tic = M / R ** 3
    r = f["rho_solar"]
    z = (r[0] - rho_tic) / (r[2] if r[0] > rho_tic else r[1])
    aR = a_over_rs(P, rho_tic)
    kk = f["rp_rs"][0]
    x = math.pi * (f["t14_h"][0] / 24) * aR / P
    out = {"a_au": a_au, "s_earth": L / a_au ** 2, "teq": teq, "mass_assumed": mass, "k_ms": k, "rho_tic": rho_tic,
           "rho_verdict": "consistent" if abs(z) < 2 else "in tension", "rho_z": z,
           "b_circ": math.sqrt(max((1 + kk) ** 2 - x ** 2, 0.0)),
           "next": [x for x in nxt if x["night_frac"] > 0.05][:n_show], "all_next": nxt}
    if dv_check:
        out.update(dv_check(f, T, P, ref))
    return out
