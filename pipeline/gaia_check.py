"""ESA Gaia DR3 checks for every star the report names, from the Gaia archive (gea.esac.esa.int, TAP).

Each TIC star is matched to the Gaia DR3 source within 15" whose position, moved back to J2000 with its proper
motion, is closest to the TIC position (TIC positions are J2000) and whose G magnitude fits the TESS magnitude.
Then, per star:
  - gaia_source: RUWE (above 1.4: the astrometry does not fit a single star, often an unresolved companion),
    ipd_frac_multi_peak (a second peak in Gaia's image windows: a close resolved double), non_single_star,
    radial-velocity scatter (rv_chisq_pvalue < 0.01 with rv_renormalised_gof > 4: RV variable), phot_variable_flag;
  - nss_two_body_orbit: Gaia's orbital solutions (astrometric, spectroscopic, eclipsing), with periods;
  - vari_eclipsing_binary: Gaia's eclipsing-binary model, with its period;
  - vari_classifier_result: Gaia's best variability class.
For the Track C pixel test, the neighbour that an eclipse set was traced to is also looked up, to tell a chance
blend from a bound wide pair (same parallax, tangential velocity difference of a few km/s).

  python gaia_check.py
Output: results/candidates/gaia_checks.csv, gaia_nss.csv, gaia_neighbours.csv
"""
import io
import json
import math
import re
import time
import urllib.parse
import urllib.request

import env  # noqa: F401
import numpy as np
import pandas as pd

import ticdb
from env import RESULTS

C = RESULTS / "candidates"
TAP = "https://gea.esac.esa.int/tap-server/tap/sync"
COLS = ("source_id, ra, dec, pmra, pmdec, pmra_error, pmdec_error, parallax, parallax_error, phot_g_mean_mag, "
        "bp_rp, ruwe, ipd_frac_multi_peak, ipd_gof_harmonic_amplitude, non_single_star, radial_velocity, "
        "radial_velocity_error, rv_chisq_pvalue, rv_renormalised_gof, rv_nb_transits, phot_variable_flag")


def tap(q, tries=4):
    data = urllib.parse.urlencode({"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv", "QUERY": q}).encode()
    for k in range(tries):
        try:
            with urllib.request.urlopen(urllib.request.Request(TAP, data=data), timeout=300) as r:
                return pd.read_csv(io.BytesIO(r.read()))
        except Exception:
            if k == tries - 1:
                raise
            time.sleep(5 * (k + 1))


def match(ra, dec, tmag, radius=15.0):
    g = tap(f"SELECT {COLS} FROM gaiadr3.gaia_source WHERE 1=CONTAINS(POINT('ICRS', ra, dec), "
            f"CIRCLE('ICRS', {ra}, {dec}, {radius / 3600}))")
    if g.empty:
        return None
    pmra, pmdec = g.pmra.fillna(0), g.pmdec.fillna(0)
    # back to J2000 (Gaia DR3 epoch is 2016.0)
    ra0 = g.ra - 16 * pmra / 3.6e6 / np.cos(np.radians(g.dec))
    de0 = g.dec - 16 * pmdec / 3.6e6
    sep = np.hypot((ra0 - ra) * math.cos(math.radians(dec)), de0 - dec) * 3600
    g = g.assign(sep_j2000=sep)
    # T is about G - 0.4 for sun-like stars; allow a generous range
    ok = g[(g.sep_j2000 < 3) & ((g.phot_g_mean_mag - tmag).between(-0.5, 2.0))]
    if ok.empty:
        return None
    return ok.sort_values("sep_j2000").iloc[0]


def ids_in(ids):
    return "(" + ",".join(str(int(i)) for i in ids) + ")"


def trackc_neighbours():
    """The neighbour each Track C eclipse set was traced to by eb_pixels.py, by G magnitude and separation."""
    p = C / "eb_pixels.json"
    if not p.exists():
        return []
    out = []
    for tic, o in json.load(open(p)).items():
        for j, s in o.get("sets", {}).items():
            alt = s.get("best_alternative")
            if not alt or not str(s.get("verdict", "")).startswith("favours"):
                continue
            m = re.match(r'G([\d.]+) at (\d+)"', alt)
            if m:
                out.append({"tic": int(tic), "set": int(j), "name": alt, "g": float(m.group(1)), "sep": float(m.group(2))})
    return out


def main():
    from exofop_check import named_stars
    tics = set(named_stars())
    p = C / "eb_pixels.json"
    if p.exists():
        tics |= {int(t) for t in json.load(open(p))}
    p = C / "kelt_targets.json"
    if p.exists():
        tics |= {int(x["tic"]) for x in json.load(open(p))}
    tics = sorted(tics)
    st = ticdb.query(tics, ["ra", "dec", "Tmag"]).set_index("tic")
    rows = []
    for i, tic in enumerate(tics):
        if tic not in st.index:
            rows.append({"tic": tic, "gaia_match": "no TIC row"})
            continue
        s = st.loc[tic]
        try:
            g = match(float(s.ra), float(s.dec), float(s.Tmag))
        except Exception as ex:
            rows.append({"tic": tic, "gaia_match": "error " + repr(ex)[:80]})
            continue
        if g is None:
            rows.append({"tic": tic, "gaia_match": "none"})
            continue
        r = {"tic": tic, "gaia_match": "ok", "Tmag": float(s.Tmag)}
        r.update({k: (g[k] if not isinstance(g[k], float) or math.isfinite(g[k]) else None) for k in g.index})
        rows.append(r)
        if i % 20 == 0:
            print(i, tic, r.get("source_id"), r.get("ruwe"), flush=True)
    d = pd.DataFrame(rows)
    sid = d.source_id.dropna().astype("int64").tolist()
    nss, veb, cls = [], [], []
    for k in range(0, len(sid), 100):
        chunk = ids_in(sid[k:k + 100])
        nss.append(tap("SELECT source_id, nss_solution_type, period, period_error, eccentricity, significance, "
                       f"goodness_of_fit FROM gaiadr3.nss_two_body_orbit WHERE source_id IN {chunk}"))
        veb.append(tap("SELECT source_id, frequency, model_type, global_ranking, derived_primary_ecl_depth, "
                       f"derived_secondary_ecl_depth FROM gaiadr3.vari_eclipsing_binary WHERE source_id IN {chunk}"))
        cls.append(tap("SELECT source_id, best_class_name, best_class_score FROM gaiadr3.vari_classifier_result "
                       f"WHERE source_id IN {chunk}"))
    nss, veb, cls = (pd.concat(x, ignore_index=True) for x in (nss, veb, cls))
    nss = nss.merge(d[["tic", "source_id"]], on="source_id")
    nss.to_csv(C / "gaia_nss.csv", index=False)
    veb = veb.assign(veb_period=1 / veb.frequency)
    d = d.merge(veb[["source_id", "veb_period", "model_type", "derived_primary_ecl_depth"]], on="source_id", how="left")
    d = d.merge(cls, on="source_id", how="left")
    agg = nss.groupby("source_id").apply(
        lambda g: "; ".join(f"{a} P={p:.4g} d" if p == p else str(a) for a, p in zip(g.nss_solution_type, g.period)),
        include_groups=False).rename("nss_orbits")
    d = d.merge(agg, left_on="source_id", right_index=True, how="left")
    d["rv_variable"] = (d.rv_chisq_pvalue < 0.01) & (d.rv_renormalised_gof > 4) & (d.rv_nb_transits >= 10)
    d["flags"] = d.apply(flags, axis=1)
    d.to_csv(C / "gaia_checks.csv", index=False)
    print(len(d), "stars;", int((d.gaia_match == "ok").sum()), "matched;", int((d.ruwe > 1.4).sum()), "RUWE>1.4;",
          len(nss), "orbits;", int(d.veb_period.notna().sum()), "Gaia EBs")
    # neighbours blamed for Track C eclipse sets: co-moving or not?
    nb = []
    for x in trackc_neighbours():
        if x["tic"] not in st.index:
            continue
        s = st.loc[x["tic"]]
        g = tap(f"SELECT {COLS} FROM gaiadr3.gaia_source WHERE 1=CONTAINS(POINT('ICRS', ra, dec), "
                f"CIRCLE('ICRS', {float(s.ra)}, {float(s.dec)}, {(x['sep'] + 3) / 3600})) "
                f"AND phot_g_mean_mag BETWEEN {x['g'] - 0.06} AND {x['g'] + 0.06}")
        if g.empty:
            continue
        n = g.iloc[0]
        t = d[d.tic == x["tic"]]
        t = t.iloc[0] if len(t) else None
        r = dict(x, neighbour_source_id=int(n.source_id), n_parallax=n.parallax, n_parallax_error=n.parallax_error,
                 n_pmra=n.pmra, n_pmdec=n.pmdec, n_ruwe=n.ruwe)
        if t is not None and t.get("gaia_match") == "ok":
            dplx = (n.parallax - t.parallax) / math.hypot(n.parallax_error, t.parallax_error)
            dmu = math.hypot(n.pmra - t.pmra, n.pmdec - t.pmdec)
            plx = max(t.parallax, 0.05)
            r.update(t_parallax=t.parallax, t_parallax_error=t.parallax_error, t_pmra=t.pmra, t_pmdec=t.pmdec,
                     parallax_diff_sigma=dplx, dpm_masyr=dmu, dvt_kms=4.74 * dmu / plx,
                     comoving=bool(abs(dplx) < 3 and 4.74 * dmu / plx < 5))
        nb.append(r)
    pd.DataFrame(nb).to_csv(C / "gaia_neighbours.csv", index=False)
    print(len(nb), "neighbours looked up;", sum(1 for r in nb if r.get("comoving")), "co-moving")


def flags(r):
    f = []
    if r.get("gaia_match") != "ok":
        return "no Gaia match"
    if r.ruwe == r.ruwe and r.ruwe > 1.4:
        f.append(f"RUWE {r.ruwe:.1f}")
    if r.ipd_frac_multi_peak == r.ipd_frac_multi_peak and r.ipd_frac_multi_peak > 10:
        f.append(f"double peak in {r.ipd_frac_multi_peak:.0f}% of Gaia scans")
    if r.rv_variable:
        f.append("RV variable")
    if isinstance(r.nss_orbits, str):
        f.append("Gaia orbit: " + r.nss_orbits)
    if r.veb_period == r.veb_period:
        f.append(f"Gaia EB P={r.veb_period:.4f} d")
    if isinstance(r.best_class_name, str):
        f.append(f"Gaia class {r.best_class_name}")
    return "; ".join(f)


if __name__ == "__main__":
    main()
