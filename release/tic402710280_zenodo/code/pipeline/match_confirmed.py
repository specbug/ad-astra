"""Cross-match every survey candidate with NASA's confirmed-planet table (known.flag_confirmed).

Adds known_planet / known_planet_P / known_planet_sep to trackb_final.csv, trackc_candidates.csv and
trackd_final.csv, and writes every match to results/candidates/confirmed_matches.csv.

  python match_confirmed.py [--dry-run]
"""
import sys

import env  # noqa: F401
import pandas as pd

import known
from env import RESULTS

C = RESULTS / "candidates"
FILES = {"B": "trackb_final.csv", "C": "trackc_candidates.csv", "D": "trackd_final.csv"}


def main(dry=False):
    ps = known.confirmed_planets()
    if ps is None:
        sys.exit("no confirmed-planet pages in catalogs/hf_ps/")
    print(f"confirmed-planet table: {len(ps)} planets, rows {ps.row.min()}..{ps.row.max()}")
    out = []
    for track, name in FILES.items():
        d = pd.read_csv(C / name)
        d = known.flag_confirmed(d.drop(columns=[c for c in ("known_planet", "known_planet_P", "known_planet_sep") if c in d]))
        m = d[d.known_planet.notna()]
        keep = [c for c in ("tic", "Tmag", "verdict", "is_toi", "known_planet", "known_planet_P", "known_planet_sep") if c in m]
        out.append(m[keep].assign(track=track))
        print(f"Track {track}: {len(d)} stars, {len(m)} with a confirmed planet within 10″ or by TIC name "
              f"({int((~m.is_toi.astype(bool)).sum()) if 'is_toi' in m else '?'} not TOIs)")
        if not dry:
            d.to_csv(C / name, index=False, float_format="%.6g")
    res = pd.concat(out, ignore_index=True)
    if not dry:
        res.to_csv(C / "confirmed_matches.csv", index=False, float_format="%.6g")
    show = res[~res.is_toi.astype(bool)] if "is_toi" in res else res
    print(show.to_string(index=False, max_colwidth=60))


if __name__ == "__main__":
    main(dry="--dry-run" in sys.argv)
