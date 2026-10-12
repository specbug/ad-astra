"""Full-text literature search in NASA ADS for every star the report names.

Each star is searched by all its identifiers at once (TIC, HD, HIP, TYC, Gaia DR3, 2MASS, BD/CD/CPD, TOI; from
SIMBAD), in ADS's full-text field, which also covers tables in the body of papers. Input:
results/candidates/ads_targets.json ({tic: [identifiers]}); output: results/candidates/ads_mentions.json.

Needs the environment variable ADS_API_KEY (never printed). Standard library only.

  python ads_check.py
"""
import json
import os
import pathlib
import sys
import time
import urllib.parse
import urllib.request

C = pathlib.Path(__file__).resolve().parent.parent / "results" / "candidates"
API = "https://api.adsabs.harvard.edu/v1/search/query"


def variants(ident):
    """'HD 101544' -> also 'HD101544' (papers write both)."""
    out = {ident}
    head, _, tail = ident.partition(" ")
    if head in ("TIC", "HD", "HIP", "TOI", "TYC", "BD", "CD", "CPD") and tail:
        out.add(head + tail.replace(" ", ""))
    return sorted(out)


def search(q, token, rows=200):
    url = API + "?" + urllib.parse.urlencode({"q": q, "fl": "bibcode,title,year,pub,doctype", "rows": rows, "sort": "date desc"})
    req = urllib.request.Request(url, headers={"Authorization": "Bearer " + token})
    for attempt in range(4):
        try:
            with urllib.request.urlopen(req, timeout=60) as r:
                return json.loads(r.read())["response"]
        except urllib.error.HTTPError as e:
            if e.code == 429 or e.code >= 500:
                time.sleep(10 * (attempt + 1))
                continue
            raise
        except urllib.error.URLError:  # connection resets
            time.sleep(10 * (attempt + 1))
    raise RuntimeError("ADS kept failing for: " + q[:80])


def main():
    token = os.environ.get("ADS_API_KEY")
    if not token:
        sys.exit("ADS_API_KEY is not set")
    targets = json.load(open(C / "ads_targets.json"))
    out = {}
    for i, (tic, ids) in enumerate(targets.items()):
        terms = sorted({v for ident in ids for v in variants(ident)})
        q = "full:(" + " OR ".join(f'"{t}"' for t in terms) + ")"
        res = search(q, token)
        papers = [{"bibcode": d["bibcode"], "year": d.get("year"), "title": (d.get("title") or [""])[0],
                   "pub": d.get("pub"), "doctype": d.get("doctype")} for d in res.get("docs", [])]
        # which identifiers are actually named, one query each, only for stars with hits
        per_id = {}
        if papers:
            for t in terms:
                r = search(f'full:"{t}"', token, rows=1)
                per_id[t] = r.get("numFound", 0)
                time.sleep(0.3)
        out[tic] = {"ids": ids, "n": res.get("numFound", 0), "per_id": per_id, "papers": papers}
        print(i, tic, out[tic]["n"], flush=True)
        time.sleep(0.3)
    json.dump(out, open(C / "ads_mentions.json", "w"), indent=1)
    print("done:", len(out), "stars;", sum(1 for v in out.values() if v["n"]), "with papers")


if __name__ == "__main__":
    main()
