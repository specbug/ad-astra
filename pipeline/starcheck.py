"""Star-level check of a single-dip candidate against every sector of data for that star.

Runs the survey dip finder on all available light curves (SPOC 2-min > TESS-SPOC > QLP per sector) and
compares every significant dip with the candidate:
  - deeper dips elsewhere (> 1.6x the candidate depth): an eclipsing binary whose secondary eclipse
    (or a diluted eclipse) masquerades as a transit
  - dips with consistent depth and duration: repeat transits (duo-transit / period constraints)
"""
import numpy as np

import dossier
from survey import process_lc


def check_star(tic, cand_depth, cand_dur_h, cand_t0s=(), snr_min=8.0):
    secs = dossier.all_sectors(tic)
    lcs = dossier.best_lcs(tic, secs)
    deeper, similar, all_ev = [], [], []
    for lc in lcs:
        try:
            _, evs = process_lc(lc, snr_min=snr_min)
        except Exception:
            continue
        for e in evs:
            if e["snr"] < snr_min or e.get("local_snr", 0) < 6:
                continue
            rec = {"sector": e["sector"], "product": e["product"], "t0": e["t0"], "depth": e["depth"],
                   "dur_h": e["dur_h"], "snr": e["snr"], "local_snr": e.get("local_snr")}
            all_ev.append(rec)
            if any(abs(e["t0"] - t) < 0.3 for t in cand_t0s):
                continue
            if e["depth"] > 1.6 * cand_depth:
                deeper.append(rec)
            elif abs(e["depth"] / cand_depth - 1) < 0.4 and abs(np.log(e["dur_h"] / cand_dur_h)) < 0.6:
                similar.append(rec)
    return {"n_sectors": len(lcs), "sectors": [l["meta"]["sector"] for l in lcs],
            "n_deeper": len(deeper), "deeper": deeper[:20], "n_similar": len(similar), "similar": similar[:20],
            "n_events": len(all_ev)}, lcs
