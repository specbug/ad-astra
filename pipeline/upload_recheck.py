"""Live ExoFOP recheck of every target in the draft upload, immediately before it is used: has anyone posted a
matching period (within 0.1%, or an integer ratio 2 or 3) or created a TOI on the star since our last look?

Reads results/candidates/exofop_upload/params_planet_DRAFT.txt; fetches each host's target page fresh (no cache);
writes results/candidates/exofop_upload/live_recheck.csv with any matches, and the date of the check.

  python upload_recheck.py
"""
import time
from datetime import datetime, timezone

import env  # noqa: F401
import pandas as pd
import requests

from env import RESULTS

D = RESULTS / "candidates" / "exofop_upload"
URL = "https://exofop.ipac.caltech.edu/tess/target.php?id={}&json"


def tic_of(target, tois):
    if target.startswith("TIC"):
        return int(target[3:].split(".")[0])
    if target.startswith("TOI"):
        return tois.get(target[3:])
    return None


def main():
    df = pd.read_csv(D / "params_planet_DRAFT.txt", sep="|", comment="\\", dtype=str)
    from env import CATALOGS
    toi = pd.read_csv(CATALOGS / "exofop" / "toi_20261012.csv", low_memory=False)
    tois = {f"{float(t):.2f}": int(tic) for t, tic in zip(toi["TOI"], toi["TIC ID"])}
    out = []
    now = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    for _, r in df.iterrows():
        tic = tic_of(r["target"], tois)
        P = float(r["period"]) if isinstance(r["period"], str) and r["period"] else None
        rec = {"target": r["target"], "tic": tic, "our_period": P, "checked": now}
        try:
            d = requests.get(URL.format(tic), timeout=60).json()
        except Exception as e:
            rec["error"] = repr(e)[:100]
            out.append(rec)
            continue
        posts = []
        ephem = str(r.get("notes", "")).startswith("Alt. ephemeris")
        T = float(r["epoch"]) if isinstance(r["epoch"], str) and r["epoch"] else None
        for e in d.get("planet_parameters", []) or []:
            try:
                per = float(e.get("per"))
            except (TypeError, ValueError):
                continue
            if not (per > 0 and P):
                continue
            who = f"{e.get('name') or e.get('toi')} P={per} by {e.get('puser')} {str(e.get('pdate'))[:10]}"
            if ephem:
                # an ephemeris correction is already posted if someone's ephemeris puts the 1 Dec 2026 transit
                # within an hour of ours
                try:
                    ep = float(e.get("epoch"))
                except (TypeError, ValueError):
                    continue
                ref = 2461375.5
                ours = T + round((ref - T) / P) * P
                theirs = ep + round((ours - ep) / per) * per
                if abs(theirs - ours) * 24 < 1.0:
                    posts.append(who + f" (1 Dec 2026 within {abs(theirs - ours) * 1440:.0f} min)")
            else:
                for ratio in (1, 2, 3, 0.5, 1 / 3):
                    if abs(per / (P * ratio) - 1) < 1e-3:
                        posts.append(who + ("" if ratio == 1 else f" (x{ratio:g})"))
        rec["matching_posts"] = "; ".join(posts)
        rec["tois_now"] = "; ".join(str(t.get("toi")) for t in (d.get("tois") or []))
        out.append(rec)
        time.sleep(0.4)
    pd.DataFrame(out).to_csv(D / "live_recheck.csv", index=False)
    hits = [o for o in out if o.get("matching_posts")]
    print(len(out), "targets rechecked;", len(hits), "with a matching period already posted")
    for o in hits:
        print("  ", o["target"], "|", o["matching_posts"])


if __name__ == "__main__":
    main()
