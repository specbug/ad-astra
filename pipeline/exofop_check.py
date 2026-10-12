"""What ExoFOP already holds for each candidate star: TOIs, community TOIs (CTOIs), every user-uploaded set
of planet parameters (with its period, author and date), and follow-up observations.

Target pages are fetched as JSON (target.php?id=<TIC>&json), cached under catalogs/exofop/targets/, and
summarised into results/candidates/exofop_status.csv.

  python exofop_check.py [--named] [--shortlist]      (default: both)
"""
import html
import json
import re
import sys
import time
import urllib.request

import env  # noqa: F401
import pandas as pd

from env import CATALOGS, RESULTS

CACHE = CATALOGS / "exofop" / "targets"
URL = "https://exofop.ipac.caltech.edu/tess/target.php?id={}&json"
C = RESULTS / "candidates"


def fetch(tic, pause=0.4):
    p = CACHE / f"{tic}.json"
    if p.exists() and p.stat().st_size > 100:
        return json.load(open(p))
    for attempt in range(3):
        try:
            with urllib.request.urlopen(URL.format(tic), timeout=60) as r:
                data = r.read()
            d = json.loads(data)
            p.write_bytes(data)
            time.sleep(pause)
            return d
        except Exception as ex:  # transient network errors: back off and retry
            err = ex
            time.sleep(5 * (attempt + 1))
    return {"error": repr(err)[:200]}


def summarise(tic, d):
    if "error" in d:
        return {"tic": tic, "error": d["error"]}
    uploads = [e for e in d.get("planet_parameters", []) if "name" in e and "puser" in e]
    with_p = [e for e in uploads if e.get("per") not in (None, "", "0", "0.0")]
    with_p.sort(key=lambda e: e.get("pdate") or "")
    last = with_p[-1] if with_p else {}
    confirmed = [e["name"] for e in d.get("planet_parameters", []) if e.get("puser") == "Exoplanet Archive"]
    return {
        "tic": tic,
        "star_names": d.get("basic_info", {}).get("star_names", ""),
        "tois": "; ".join(str(e.get("toi")) for e in d.get("tois", [])),
        "ctois": "; ".join(e.get("ctoi", "") for e in d.get("ctois", [])),
        "ctoi_users": "; ".join(f"{e.get('cuser')} {e.get('ctoi_created')}" for e in d.get("ctois", [])),
        "ctoi_notes": " | ".join(html.unescape(e.get("cnotes") or "") for e in d.get("ctois", [])),
        "n_uploads": len(uploads),
        "period": float(last["per"]) if last else None,
        "period_user": last.get("puser"),
        "period_date": last.get("pdate"),
        "period_note": html.unescape(last.get("pnotes") or "") if last else None,
        "confirmed": "; ".join(confirmed),
        "n_imaging": len(d.get("imaging", [])),
        "n_spectra": len(d.get("spectroscopy", [])),
        "n_timeseries": len(d.get("time_series", [])),
        "n_companions": len(d.get("stellar_companions", [])),
    }


def named_stars():
    s = html.unescape(open(RESULTS / "report" / "index.html").read())
    return sorted({int(t) for t in re.findall(r"TIC\s*(\d{5,10})", s)} | {128703021, 266680498, 266376971})


def main(args):
    import report
    CACHE.mkdir(parents=True, exist_ok=True)
    tics = set()
    if "--named" in args or len(args) == 0:
        tics |= set(named_stars())
    if "--shortlist" in args or len(args) == 0:
        tics |= set(report.trackb_shortlist(pd.read_csv(C / "trackb_final.csv")).tic.astype(int))
    import concurrent.futures as cf
    tics = sorted(tics)
    with cf.ThreadPoolExecutor(3) as ex:  # ExoFOP is slow per page; three at a time is still gentle
        data = list(ex.map(fetch, tics))
    rows = [summarise(t, d) for t, d in zip(tics, data)]
    out = pd.DataFrame(rows)
    out.to_csv(C / "exofop_status.csv", index=False)
    print(f"{len(out)} stars; with TOI {int((out.tois != '').sum())}, with CTOI {int((out.ctois != '').sum())}, "
          f"with a user-uploaded period {int(out.period.notna().sum())}, with errors {int(out.get('error', pd.Series()).notna().sum())}")


if __name__ == "__main__":
    main(sys.argv[1:])
