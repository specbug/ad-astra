"""Reduce a survey run (one or more complete sectors).

1. Systematics shared by many stars:
   - per (sector, camera, ccd) 0.1-d bins whose rate of stars-with-dips is far above that CCD's norm,
   - spacecraft-wide 0.05-d bins (all cameras together) with excess rates.
2. Blends: dips at the same time (|dt| < half a duration) in stars within 2.5 arcmin of each other.
   The star with the deepest dip in the group is the most likely source; the others are flagged.
3. Event quality flags and a per-event 'clean' summary flag used by the downstream tracks.

Usage: python survey_reduce.py SURVEY_DIR [SURVEY_DIR ...]
Writes SURVEY_DIR/events_reduced.parquet and SURVEY_DIR/hot_bins.parquet.
"""
import glob
import sys

import numpy as np
import pandas as pd

import env  # noqa: F401
from env import CATALOGS

CCD_BIN_D = 0.1
SC_BIN_D = 0.05
HOT_FACTOR = 5.0
BLEND_ARCMIN = 2.5


def load(survey_dir):
    evs = [pd.read_parquet(f) for f in sorted(glob.glob(f"{survey_dir}/events_[0-9]*.parquet"))]
    ev = pd.concat([e for e in evs if len(e)], ignore_index=True)
    st = pd.concat([pd.read_parquet(f) for f in sorted(glob.glob(f"{survey_dir}/lcstats_*.parquet"))],
                   ignore_index=True)
    if "invert" not in ev:
        ev["invert"] = False
    if "invert" not in st:
        st["invert"] = False
    return ev, st


def _poisson_hot(n, lam):
    z = 2 * (np.sqrt(n + 3 / 8) - np.sqrt(lam + 3 / 8))
    return (n >= HOT_FACTOR * lam) & (z > 4)


def hot_bins(ev, st, min_snr=6.0):
    e = ev[(ev.snr >= min_snr)].copy()
    out = []
    # per CCD
    e["bin"] = np.floor(e.t0 / CCD_BIN_D).astype(np.int64)
    g = e.groupby(["invert", "sector", "camera", "ccd", "bin"]).tic.nunique().rename("n").reset_index()
    med = g.groupby(["invert", "sector", "camera", "ccd"]).n.transform("median")
    # the norm: median over bins, floored at 1 and at the CCD's mean over all time bins of the sector
    ns = st[st.status == "ok"].groupby(["invert", "sector", "camera", "ccd"]).size().rename("nstars")
    g = g.join(ns, on=["invert", "sector", "camera", "ccd"])
    lam = np.maximum(med, 1.0)
    g["hot"] = _poisson_hot(g.n, lam)
    g["kind"] = "ccd"
    out.append(g)
    # spacecraft-wide (all cameras)
    e["sbin"] = np.floor(e.t0 / SC_BIN_D).astype(np.int64)
    s = e.groupby(["invert", "sector", "sbin"]).tic.nunique().rename("n").reset_index()
    lam = np.maximum(s.groupby(["invert", "sector"]).n.transform("median"), 1.0)
    s["hot"] = _poisson_hot(s.n, lam)
    s["kind"] = "spacecraft"
    out.append(s)
    return out[0], out[1]


def flag_systematics(ev, ccd_hot, sc_hot):
    ev = ev.copy()
    ev["bin"] = np.floor(ev.t0 / CCD_BIN_D).astype(np.int64)
    h = ccd_hot[ccd_hot.hot][["invert", "sector", "camera", "ccd", "bin"]]
    h = pd.concat([h, h.assign(bin=h.bin - 1), h.assign(bin=h.bin + 1)]).drop_duplicates()
    h["sys_ccd"] = True
    ev = ev.merge(h, on=["invert", "sector", "camera", "ccd", "bin"], how="left")
    ev["sbin"] = np.floor(ev.t0 / SC_BIN_D).astype(np.int64)
    s = sc_hot[sc_hot.hot][["invert", "sector", "sbin"]]
    s = pd.concat([s, s.assign(sbin=s.sbin - 1), s.assign(sbin=s.sbin + 1)]).drop_duplicates()
    s["sys_sc"] = True
    ev = ev.merge(s, on=["invert", "sector", "sbin"], how="left")
    for c in ("sys_ccd", "sys_sc"):
        ev[c] = ev[c].astype("boolean").fillna(False).astype(bool)
    ev["sys_hot"] = ev.sys_ccd | ev.sys_sc
    return ev.drop(columns=["bin", "sbin"])


def flag_blends(ev, coords):
    """Mark events coincident in time with a dip in a star within BLEND_ARCMIN."""
    ev = ev.merge(coords, on="tic", how="left")
    ev["blend_n"] = 0
    ev["blend_src"] = np.nan
    ev["blend_is_deepest"] = True
    strong = ev[(ev.snr >= 6) & ~ev.sys_hot].copy()
    for (inv, sec), grp in strong.groupby(["invert", "sector"]):
        if len(grp) < 2:
            continue
        grp = grp.sort_values("t0")
        t = grp.t0.values
        ra = np.radians(grp.ra.values)
        dec = np.radians(grp.dec.values)
        dur = grp.dur_h.values / 24
        dep = grp.depth.values
        tic = grp.tic.values
        idx = grp.index.values
        j0 = 0
        for i in range(len(grp)):
            while t[i] - t[j0] > 1.0:
                j0 += 1
            j1 = np.searchsorted(t, t[i] + 1.0)
            cand = np.arange(j0, j1)
            cand = cand[(cand != i) & (tic[cand] != tic[i])]
            if cand.size == 0:
                continue
            dt = np.abs(t[cand] - t[i])
            cand = cand[dt < 0.5 * np.maximum(dur[cand], dur[i])]
            if cand.size == 0:
                continue
            sep = np.degrees(np.arccos(np.clip(np.sin(dec[i]) * np.sin(dec[cand]) +
                                               np.cos(dec[i]) * np.cos(dec[cand]) * np.cos(ra[i] - ra[cand]), -1, 1))) * 60
            cand = cand[sep < BLEND_ARCMIN]
            if cand.size == 0:
                continue
            ev.loc[idx[i], "blend_n"] = cand.size
            deepest = cand[np.argmax(dep[cand])]
            if dep[deepest] > dep[i]:
                ev.loc[idx[i], "blend_is_deepest"] = False
                ev.loc[idx[i], "blend_src"] = tic[deepest]
    return ev


def quality_flags(ev):
    ev = ev.copy()
    D = ev.dur_h / 24
    ev["f_edge"] = ev.edge_dist_d < np.maximum(0.25, D)
    ev["f_cov"] = ev.coverage < 0.8
    ev["f_bkg"] = ev.z_bkg > 4          # background rises during the dip: over-subtraction artefact
    ev["f_cen"] = (ev.z_cx.abs() > 8) | (ev.z_cy.abs() > 8)
    # the dip must be in the raw (SAP) flux too, at >= 40% of the flattened depth
    ev["f_sap"] = ~(ev.depth_sap_flux > 0.4 * ev.depth)
    ev["f_local"] = ~(ev.local_snr >= 5)
    ev["clean"] = ~(ev.sys_hot | ev.f_edge | ev.f_cov | ev.f_bkg | ev.f_sap | ev.f_local) & ev.blend_is_deepest
    return ev


def reduce_dirs(dirs):
    ev, st = zip(*[load(d) for d in dirs])
    ev = pd.concat(ev, ignore_index=True)
    st = pd.concat(st, ignore_index=True)
    ccd_hot, sc_hot = hot_bins(ev, st)
    ev = flag_systematics(ev, ccd_hot, sc_hot)
    coords = pd.read_parquet(CATALOGS / "tic_master.parquet", columns=["tic", "ra", "dec"])
    coords = coords[coords.tic.isin(ev.tic.unique())]
    ev = flag_blends(ev, coords)
    ev = quality_flags(ev)
    return ev, st, ccd_hot, sc_hot


if __name__ == "__main__":
    dirs = sys.argv[1:]
    ev, st, ccd_hot, sc_hot = reduce_dirs(dirs)
    print("light curves:", len(st), st.status.value_counts().to_dict())
    print("events:", len(ev), " systematics:", int(ev.sys_hot.sum()), " clean:", int(ev.clean.sum()))
    for inv in (False, True):
        e = ev[(ev.invert == inv) & ev.clean]
        n_lc = int(((st.invert == inv) & (st.status == "ok")).sum())
        if n_lc:
            print(f"invert={inv}: {n_lc} LCs, clean events per 1000 LCs at snr>=7/10/15/25:",
                  [round(1000 * (e.snr >= x).sum() / n_lc, 2) for x in (7, 10, 15, 25)])
    out = dirs[0]
    ev.to_parquet(f"{out}/events_reduced.parquet", index=False)
    pd.concat([ccd_hot, sc_hot]).to_parquet(f"{out}/hot_bins.parquet", index=False)
