"""Literature check of the stars behind the new periods of the period review (period_review.py), in NASA ADS.

A period found here counts as new only if no paper has published it. For each star:
  1. identifiers from SIMBAD (all names of the closest SIMBAD object within 3" of the TIC position), plus the
     TIC name and the TOI or CTOI designation;
  2. ADS full-text search on all of them at once (ads_check.search), which also covers tables;
  3. for stars with papers, a second full-text query that adds the period itself (the alias to two decimals,
     e.g. "29.76", and to one decimal), so papers that print the period are flagged.

Needs the environment variable ADS_API_KEY (never printed).

  python ads_periods.py
Output: results/candidates/ads_periods.json, keyed "tic|id"
"""
import json
import os
import sys
import time
import urllib.parse
import urllib.request

import env  # noqa: F401

import ticdb
from ads_check import search, variants
from env import RESULTS
from period_vet import new_periods

C = RESULTS / "candidates"
OUT = C / "ads_periods.json"
TAP = "https://simbad.cds.unistra.fr/simbad/sim-tap/sync"
KEEP = ("HD ", "HIP ", "TYC ", "2MASS ", "Gaia DR3 ", "BD", "CD", "CPD", "GJ ", "LHS ", "UCAC4 ", "TOI", "KELT", "WASP",
        "HAT", "K2", "EPIC ", "KIC ", "Kepler", "NGTS", "LP ", "Wolf ", "Ross ", "G ", "PM ", "LSPM ", "NLTT ",
        "WISE", "CoRoT", "ASASSN", "ATO ", "ZTF", "Gaia DR2 ")


def simbad_ids(ra, dec, radius=3.0):
    q = f"""SELECT b.main_id, i.id, DISTANCE(POINT('ICRS', b.ra, b.dec), POINT('ICRS', {ra}, {dec})) AS d
            FROM basic AS b JOIN ident AS i ON i.oidref = b.oid
            WHERE CONTAINS(POINT('ICRS', b.ra, b.dec), CIRCLE('ICRS', {ra}, {dec}, {radius / 3600})) = 1"""
    data = urllib.parse.urlencode({"request": "doQuery", "lang": "adql", "format": "json", "query": q}).encode()
    for attempt in range(4):
        try:
            with urllib.request.urlopen(urllib.request.Request(TAP, data=data), timeout=60) as r:
                rows = json.loads(r.read())["data"]
            break
        except Exception:
            time.sleep(5 * (attempt + 1))
    else:
        return None, []
    if not rows:
        return None, []
    best = min(rows, key=lambda r: r[2])[0]
    ids = sorted({" ".join(r[1].split()) for r in rows if r[0] == best})
    return best, [i for i in ids if i.startswith(KEEP)]


def names_for(o):
    """TIC, TOI/CTOI names in the forms papers use."""
    tic = int(o["tic"])
    out = [f"TIC {tic}"]
    i = str(o["id"])
    if i.startswith("TOI-"):
        base = i.split(".")[0]
        out += [base, base.replace("-", " "), i]
    elif i.startswith("CTOI"):
        out += [f"TIC {i.split()[1]}"]   # CTOIs are written as TIC NNN.01
    return out


def period_terms(aliases):
    out = set()
    for P in aliases[:3]:
        out.add(f"{P:.2f}")
        out.add(f"{P:.1f}")
    return sorted(out)


def main():
    token = os.environ.get("ADS_API_KEY")
    if not token:
        sys.exit("ADS_API_KEY is not set")
    rv = json.load(open(C / "period_review.json"))
    res = json.load(open(OUT)) if OUT.exists() else {}
    todo = new_periods(rv)
    pos = ticdb.query([int(o["tic"]) for o in todo], ["ra", "dec"]).set_index("tic")
    for o in todo:
        key = f"{o['tic']}|{o['id']}"
        if key in res:
            continue
        ra, dec = float(pos.loc[int(o["tic"]), "ra"]), float(pos.loc[int(o["tic"]), "dec"])
        main_id, ids = simbad_ids(ra, dec)
        names = sorted(set(names_for(o) + ids))
        terms = sorted({v for n in names for v in variants(n)})
        q = "full:(" + " OR ".join(f'"{t}"' for t in terms) + ")"
        r = search(q, token)
        papers = [{"bibcode": d["bibcode"], "year": d.get("year"), "title": (d.get("title") or [""])[0],
                   "pub": d.get("pub"), "doctype": d.get("doctype")} for d in r.get("docs", [])]
        al = o.get("aliases")
        al = json.loads(al) if isinstance(al, str) else (al or [])
        with_period = []
        if papers and al:
            pt = period_terms(al)
            r2 = search(q + " AND full:(" + " OR ".join(f'"{t}"' for t in pt) + ")", token)
            with_period = [d["bibcode"] for d in r2.get("docs", [])]
        res[key] = {"tic": int(o["tic"]), "id": o["id"], "simbad": main_id, "names": names, "n": r.get("numFound", 0),
                    "papers": papers, "period_terms": period_terms(al) if al else [], "papers_with_period": with_period}
        print(key, main_id, res[key]["n"], len(with_period), flush=True)
        json.dump(res, open(OUT, "w"), indent=1)
        time.sleep(0.5)
    n = sum(1 for v in res.values() if v["n"])
    print("done:", len(res), "candidates;", n, "with papers;", sum(1 for v in res.values() if v["papers_with_period"]),
          "with a paper that also prints the period")


if __name__ == "__main__":
    main()
