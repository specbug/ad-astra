"""Turn reduced survey events into per-track candidate lists.

Track B  single / duo transits: isolated, clean, transit-shaped dips with planet-sized implied radii.
Track C  eclipsing binaries to analyse for extra eclipses (light curves with several strong periodic dips).
Track D  anomalies: deep non-periodic dips, very long or strongly asymmetric dips around dwarf stars.
"""
import numpy as np
import pandas as pd

R_EARTH_PER_RSUN = 109.1
R_JUP_PER_RSUN = 9.73


def period_from_duration(dur_h, rho_sun):
    """Circular, central-transit period implied by a duration (days); rough, b=0 gives the minimum."""
    return 365.25 * (dur_h / 13.0) ** 3 * rho_sun


def attach_star(ev, tic):
    cols = ["tic", "Tmag", "Teff", "logg", "rad", "mass", "rho", "lumclass", "d", "contratio", "ra", "dec"]
    s = tic[[c for c in cols if c in tic]].copy()
    if "ra" in ev:
        s = s.drop(columns=["ra", "dec"])
    return ev.merge(s, on="tic", how="left")


def lc_context(ev):
    """Per light curve: number of clean dips and the strongest clean dip, used for uniqueness."""
    c = ev[ev.clean]
    g = c.groupby(["invert", "tic", "sector"]).agg(n_clean=("snr", "size"), snr_max=("snr", "max"))
    return ev.join(g, on=["invert", "tic", "sector"])


def track_b(ev, min_snr=7.5):
    e = ev[ev.clean & (ev.snr >= min_snr) & (ev.local_snr >= 6) & (ev.dur_h >= 1.0) & (ev.depth < 0.25)].copy()
    # uniqueness within its light curve: no other clean dip of comparable strength unless it looks the same
    e = e[(e.n_clean <= 2)]
    rstar = e.rad.fillna(1.0)
    dil = 1 + e.contratio.fillna(0).clip(0, 10)
    e["rp_rearth"] = np.sqrt(e.depth * dil) * rstar * R_EARTH_PER_RSUN
    e["rp_rjup"] = e.rp_rearth / 11.2
    rho = e.rho.where(e.rho > 0, e.mass / e.rad ** 3)
    e["p_dur_d"] = period_from_duration(e.dur_h, rho.fillna(1.0))
    e["shape_vu"] = (e.depth_core - e.depth_wing) / e.depth.where(e.depth > 0)
    e["dwarf"] = (e.rad.fillna(1) < 2.0) & (e.logg.fillna(4.4) > 3.6)
    e["planet_sized"] = e.rp_rjup < 2.2
    return e


def group_track_b(e):
    """Group candidate single events per star across sectors."""
    rows = []
    for tic, g in e[~e.invert].groupby("tic"):
        g = g.sort_values("t0")
        best = g.loc[g.snr.idxmax()]
        # consistent repeats: depth within 35% and duration within one grid step of the best event
        same = g[(np.abs(g.depth / best.depth - 1) < 0.35) & (np.abs(np.log(g.dur_h / best.dur_h)) < 0.45)]
        rows.append({"tic": tic, "n_events": len(g), "n_consistent": len(same),
                     "t0s": ";".join(f"{x:.4f}" for x in same.t0), "sectors": ";".join(str(int(s)) for s in g.sector),
                     "best_t0": best.t0, "best_snr": best.snr, "best_local_snr": best.local_snr,
                     "depth": best.depth, "dur_h": best.dur_h, "rp_rearth": best.rp_rearth,
                     "p_dur_d": best.p_dur_d, "shape_vu": best.shape_vu, "Tmag": best.Tmag, "Teff": best.Teff,
                     "rad": best.rad, "dwarf": best.dwarf, "planet_sized": best.planet_sized,
                     "contratio": best.contratio})
    return pd.DataFrame(rows)


def track_c_targets(ev, min_events=3, min_snr=10):
    c = ev[ev.clean & (ev.snr >= min_snr) & ~ev.invert]
    g = c.groupby(["tic", "sector"]).agg(n=("snr", "size"), snr=("snr", "max"), depth=("depth", "max"),
                                         product=("product", "first")).reset_index()
    return g[g.n >= min_events]


def track_d(ev, min_snr=10):
    """Deep, non-repeating or oddly shaped dips around dwarfs."""
    c = ev[ev.clean & (ev.snr >= min_snr) & ~ev.invert].copy()
    c["asym"] = (c.depth_first - c.depth_second) / c.depth.where(c.depth > 0)
    c["dwarf"] = (c.rad.fillna(1) < 2.0) & (c.logg.fillna(4.4) > 3.6)
    long_deep = (c.dur_h >= 13) & (c.depth >= 0.01)
    asym = c.asym.abs() > 0.6
    multi = c.n_clean.between(2, 8)
    c["d_reason"] = np.select([long_deep, asym, multi], ["long_deep", "asymmetric", "multi"], "")
    return c[(c.d_reason != "") & c.dwarf]
