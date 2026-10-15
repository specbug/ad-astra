"""What is already on record for a star, from every check run on 11-12 October 2026.

Sources (all in results/candidates/):
  exofop_status.csv            ExoFOP target pages: TOIs, community TOIs (CTOIs), user-uploaded periods
  confirmed_matches.csv        NASA confirmed-planet table, position match (via match_confirmed.py)
  neighbour_planet_transits.csv  dips that fall on transits of known planets within 2'
  simbad_named.csv             SIMBAD object types and names
  vizier_named.csv             TESS and Gaia eclipsing-binary catalogues, VSX, published quadruple lists
  ads_mentions.json            ADS full-text search by every identifier (when present)
  gaia_checks.csv              ESA Gaia DR3: orbits (NSS), eclipsing-binary model, RUWE, resolved doubles, RV scatter

prior(tic) returns {"kind": ..., "notes": [...]}, kind being the strongest reason the star is not new:
"blend" (a known planet on a neighbour), "confirmed", "multiple" (a published triple/quadruple), "binary"
(catalogued eclipsing or spectroscopic binary), "ctoi" (or TOI, or a period posted by someone else), or None.
"""
import json

import pandas as pd

from env import RESULTS

C = RESULTS / "candidates"
_d = {}


def _load():
    if _d:
        return _d

    def csv(name):
        p = C / name
        return pd.read_csv(p) if p.exists() else pd.DataFrame()
    _d["exofop"] = csv("exofop_status.csv").set_index("tic") if (C / "exofop_status.csv").exists() else pd.DataFrame()
    _d["confirmed"] = csv("confirmed_matches.csv")
    _d["neighbour"] = csv("neighbour_planet_transits.csv")
    _d["simbad"] = csv("simbad_named.csv").set_index("tic") if (C / "simbad_named.csv").exists() else pd.DataFrame()
    _d["vizier"] = csv("vizier_named.csv")
    _d["ctoi_all"] = csv("exofop_ctoi_matches.csv")  # ExoFOP CTOI table of 12 Oct 2026, every candidate star
    p = C / "ads_mentions.json"
    _d["ads"] = json.load(open(p)) if p.exists() else {}
    _d["gaia"] = csv("gaia_checks.csv").set_index("tic") if (C / "gaia_checks.csv").exists() else pd.DataFrame()
    return _d


MULTIPLE_CATS = ("Kostov+2022 quadruple candidates", "Zasche+2022 multiply eclipsing", "Zasche+2022 doubly eclipsing",
                 "Zasche+2025 doubly eclipsing quadruples", "Borkovits+2025 triply eclipsing")
EB_CATS = ("TESS EBs (Prsa+2022)", "Gaia DR3 eclipsing binaries")
# confirmed planets on the star that cannot be the survey's signal
NOT_THE_SIGNAL = {266680498: "radial-velocity planets; none fits the transits"}


def prior(tic):
    d = _load()
    tic = int(tic)
    notes, kinds = [], []
    nb = d["neighbour"]
    if len(nb):
        own = nb[(nb.tic == tic)]
        far = own[own.sep_arcsec > 10]
        if len(far):
            r = far.iloc[0]
            notes.append(f"dip falls on a transit of {r.planet}, {r.sep_arcsec:.0f}″ away")
            kinds.append("blend")
    cm = d["confirmed"]
    if len(cm) and tic in set(cm.tic):
        names = str(cm[cm.tic == tic].known_planet.iloc[0])
        if tic in NOT_THE_SIGNAL:
            notes.append(f"NASA lists {names} ({NOT_THE_SIGNAL[tic]})")
        else:
            notes.append("NASA lists " + names)
            kinds.append("confirmed")
    ex = d["exofop"]
    if len(ex) and tic in ex.index:
        e = ex.loc[tic]
        if isinstance(e.get("tois"), str) and e.tois:
            notes.append(e.tois)
            kinds.append("ctoi")
        if isinstance(e.get("ctois"), str) and e.ctois:
            who, _, when = str(e.ctoi_users).partition(" ")
            notes.append(f"CTOI since {when[:7]} ({who})")
            kinds.append("ctoi")
        if pd.notna(e.get("period")) and e.get("period_user") not in ("Exoplanet Archive", "exoplanetarchive", "TESS project"):
            pn = str(e.get("period_note") or "")
            pn = f" ('{pn.strip()}')" if 0 < len(pn) <= 70 and e.period_user not in str(e.ctoi_users) else ""
            notes.append(f"period {e.period:.4f} d posted by {e.period_user}, {str(e.period_date)[:10]}{pn}")
            kinds.append("ctoi")
    elif len(d["ctoi_all"]) and tic in set(d["ctoi_all"].tic):
        r = d["ctoi_all"][d["ctoi_all"].tic == tic].iloc[0]
        notes.append(f"CTOI since {str(r['CTOI lastmod'])[:7]} ({r.User})")
        kinds.append("ctoi")
    vz = d["vizier"]
    if len(vz):
        v = vz[vz.tic == tic]
        for cat in sorted(set(v.catalog)):
            if cat in MULTIPLE_CATS:
                notes.append(cat.replace("+", " et al. ").replace(" candidates", " candidate"))
                kinds.append("multiple")
            elif cat in EB_CATS:
                close = v.loc[v.catalog == cat, "info"].str.extract(r'"_r": ([0-9.]+)', expand=False).astype(float).min()
                name = cat.split(" (")[0].replace("TESS EBs", "TESS eclipsing-binary catalogue")
                if close <= 3:
                    notes.append(name)
                    kinds.append("binary")
                else:
                    notes.append(f"{name}: an eclipsing binary {close:.0f}″ away")
    sb = d["simbad"]
    if len(sb) and tic in sb.index and float(sb.loc[tic].sep) > 3:
        ot = str(sb.loc[tic].otypes)
        if "EB*" in ot or "SB*" in ot:
            notes.append(f"SIMBAD: a binary ({sb.loc[tic].main_id.strip()}) {float(sb.loc[tic].sep):.0f}″ away")
    elif len(sb) and tic in sb.index:
        ot = str(sb.loc[tic].otypes)
        if "EB*" in ot:
            notes.append("eclipsing binary in SIMBAD")
            kinds.append("binary")
        elif "SB*" in ot:
            notes.append("spectroscopic binary in SIMBAD")
            kinds.append("binary")
        elif "**" in ot:
            notes.append("double star in SIMBAD")
    gd = d["gaia"]
    if len(gd) and tic in gd.index and gd.loc[tic].get("gaia_match") == "ok":
        g = gd.loc[tic]
        orbit = g.get("nss_orbits")
        if isinstance(orbit, str) and orbit:
            notes.append("Gaia DR3 orbit: " + orbit.replace("Orbital", "astrometric").replace("SB2C", "SB2"))
            if any(k in orbit for k in ("SB1", "SB2", "Eclipsing")):
                kinds.append("binary")
        if pd.notna(g.get("veb_period")):
            if any(n.startswith("Gaia DR3 eclipsing binaries") for n in notes):
                notes = [n + f" (P = {g.veb_period:.4f} d)" if n == "Gaia DR3 eclipsing binaries" else n for n in notes]
            else:
                notes.append(f"Gaia DR3 eclipsing binary, P = {g.veb_period:.4f} d")
                kinds.append("binary")
        if pd.notna(g.get("ruwe")) and g.ruwe > 1.4:
            notes.append(f"Gaia RUWE {g.ruwe:.1f}")
        if pd.notna(g.get("ipd_frac_multi_peak")) and g.ipd_frac_multi_peak > 10:
            notes.append(f"a second star within ~1″ in {g.ipd_frac_multi_peak:.0f}% of Gaia scans")
        if bool(g.get("rv_variable")) and not (isinstance(orbit, str) and orbit):
            notes.append("Gaia RV variable")
    a = d["ads"].get(str(tic))
    if a and a.get("n"):
        if a["n"] <= 2:  # show what the papers are, since a stray mention in an old catalogue means little
            notes.append("ADS: " + "; ".join(f"{p['year']}, \u201c{p['title'][:60]}\u2026\u201d" for p in a["papers"][:2]))
        else:
            notes.append(f"named in {a['n']} papers (ADS)")
    order = ("blend", "confirmed", "multiple", "binary", "ctoi")  # a catalogued binary says more than a CTOI does
    kind = next((k for k in order if k in kinds), None)
    return {"kind": kind, "notes": notes}


def summary(tic):
    p = prior(tic)
    return "; ".join(p["notes"]) if p["notes"] else "nothing on record"
