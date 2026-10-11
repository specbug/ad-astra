"""Track A triage: automatic red flags for multi-sector BLS signals (real and flipped).

Hard flags (evidence of an impostor): odd/even mismatch, secondary eclipse, too long for the orbit,
one or two events dominate, edge dominated, inconsistent sectors, broad depression, too deep, period
re-found after masking (variability), common period across stars (instrument).
Soft flags (context): TESS-orbit alias, weak against empirical (red) noise.

Usage: python triage_a.py OUTDIR TARGETS.csv
Reads OUTDIR/stars/*.json, re-derives each light curve, writes OUTDIR/triage.csv.
"""
import glob
import json
import sys
import multiprocessing as mp

import env  # noqa: F401
import numpy as np
import pandas as pd
from astropy.timeseries import BoxLeastSquares
from scipy.stats import chi2 as chi2_dist

from periodic import load_star, prepare_star, duration_days

TESS_ORBIT_D = 13.7
SDE_MIN, SNR_MIN, MIN_TR = 7.0, 6.0, 3
SOFT = {"tess_orbit_alias", "low_snr_red_noise"}


def same_period(p1, p2, tol=0.005):
    return any(abs(p1 / p2 - h) < tol * h for h in (0.5, 1.0, 2.0))


def box_depth(f, e, sel):
    if sel.sum() <= 3:
        return np.nan
    w = 1 / e[sel] ** 2
    return float(np.sum((1 - f[sel]) * w) / w.sum())


def diagnostics(lc, sig, rho, flipped=False):
    t, f, e, sec = lc["time"], lc["flux"].copy(), lc["flux_err"], lc["sector"]
    if flipped:
        f = 2 - f
    P, t0, dur = sig["period"], sig["t0"], sig["duration_hr"] / 24
    phase = ((t - t0 + 0.5 * P) % P) - 0.5 * P
    intr = np.abs(phase) < dur / 2
    epoch = np.round((t - t0) / P).astype(int)
    w = (1 - f) / e ** 2
    per_event = pd.Series(w[intr]).groupby(epoch[intr]).sum()
    total = per_event.clip(lower=0).sum()
    top_frac = float(per_event.max() / total) if total > 0 else 1.0

    def box_snr(sel):
        if sel.sum() <= 3:
            return 0.0
        wk = 1 / e[sel] ** 2
        return float(np.sum((1 - f[sel]) * wk) / wk.sum() * np.sqrt(wk.sum()))

    snr_all = box_snr(intr)
    snr_wo = box_snr(intr & ~np.isin(epoch, per_event.nlargest(2).index))
    offsets = np.concatenate([np.linspace(0.1, 0.4, 20), np.linspace(0.6, 0.9, 20)]) * P
    null = np.array([box_depth(f, e, np.abs(((t - t0 - o + 0.5 * P) % P) - 0.5 * P) < dur / 2) for o in offsets])
    null = null[np.isfinite(null)]
    sig_emp = 1.4826 * np.median(np.abs(null - np.median(null))) if null.size > 10 else np.nan
    depth = box_depth(f, e, intr)
    # per-sector consistency with empirical per-sector noise
    sd, ss = [], []
    for s in np.unique(sec):
        ms = sec == s
        if np.sum(intr & ms) < 3:
            continue
        sd.append(box_depth(f, e, intr & ms))
        nn = np.array([box_depth(f, e, ms & (np.abs(((t - t0 - o + 0.5 * P) % P) - 0.5 * P) < dur / 2))
                       for o in offsets])
        nn = nn[np.isfinite(nn)]
        ss.append(1.4826 * np.median(np.abs(nn - np.median(nn))) if nn.size > 10 else np.nan)
    sd, ss = np.array(sd), np.array(ss)
    ok = np.isfinite(sd) & np.isfinite(ss) & (ss > 0)
    if ok.sum() >= 3:
        ws = 1 / ss[ok] ** 2
        md = np.sum(sd[ok] * ws) / ws.sum()
        sector_p = float(chi2_dist.sf(np.sum(((sd[ok] - md) / ss[ok]) ** 2), ok.sum() - 1))
    else:
        sector_p = np.nan
    wings = (np.abs(phase) > dur) & (np.abs(phase) < 3 * dur)
    gaps = np.where(np.diff(t) > 0.5)[0]
    edges = np.concatenate([[t[0]], t[gaps], t[gaps + 1], [t[-1]]])
    near_edge = np.min(np.abs(t[:, None] - edges[None, :]), axis=1) < 0.5
    bls = BoxLeastSquares(t, f, e)
    st = bls.compute_stats(P, dur, t0)
    d_odd, e_odd = st["depth_odd"]
    d_even, e_even = st["depth_even"]
    d_half, e_half = st["depth_half"]
    return {"depth_ppm": depth * 1e6, "sigma_emp_ppm": sig_emp * 1e6,
            "snr_emp": depth / sig_emp if sig_emp > 0 else np.nan,
            "top_event_frac": top_frac, "snr_all": snr_all, "snr_wo_top2": snr_wo, "sector_p": sector_p,
            "wing_ratio": box_depth(f, e, wings) / depth if depth > 0 else np.nan,
            "edge_frac": float(near_edge[intr].mean()) if intr.any() else 1.0,
            "oddeven_sigma": abs(d_odd - d_even) / np.hypot(e_odd, e_even) if np.isfinite(e_odd + e_even) else np.nan,
            "secondary_sigma": d_half / e_half if e_half > 0 else np.nan,
            "dur_ratio": dur / duration_days(P, rho), "n_tr_data": int(np.sum(st["per_transit_count"] > 0))}


def flags(row, rad):
    fl = []
    if row["oddeven_sigma"] > 3:
        fl.append("odd_even_mismatch")
    if row["secondary_sigma"] > 3:
        fl.append("secondary_eclipse")
    if row["dur_ratio"] > 1.6:
        fl.append("long_duration")
    if row["top_event_frac"] > 0.5:
        fl.append("one_event_dominated")
    if row["snr_wo_top2"] < 0.5 * row["snr_all"]:
        fl.append("few_events_dominate")
    if row["edge_frac"] > 0.5:
        fl.append("edge_dominated")
    if row["sector_p"] < 1e-4:
        fl.append("inconsistent_sectors")
    if row["wing_ratio"] > 0.4:
        fl.append("broad_depression")
    rp = np.sqrt(max(row["depth_ppm"], 0) * 1e-6) * rad * 109.1
    if rp > 25 or row["depth_ppm"] > 3e4:
        fl.append("too_deep")
    if row["n_tr_data"] < MIN_TR:
        fl.append("few_transits")
    for k in (1, 2, 3):
        if abs(row["period"] / (k * TESS_ORBIT_D) - 1) < 0.03:
            fl.append("tess_orbit_alias")
    if row["snr_emp"] < 6:
        fl.append("low_snr_red_noise")
    return fl


def triage_star(args):
    path, rad, mass = args
    rec = json.load(open(path))
    if rec.get("status") != "ok":
        return []
    out = []
    lc = None
    for kind in ("signals", "signals_inverted"):
        sigs = rec.get(kind, [])
        for i, s in enumerate(sigs):
            if not (s["sde"] >= SDE_MIN and s["snr"] >= SNR_MIN and s["n_transits"] >= MIN_TR):
                continue
            if lc is None:
                sec = rec["sectors"]
                lc = prepare_star(load_star(rec["tic"], sec))
            try:
                d = diagnostics(lc, s, rec["rho"], flipped=(kind == "signals_inverted"))
            except ValueError as ex:  # e.g. no in-transit points left after quality cuts
                bad = f"diagnostics_failed:{str(ex)[:40]}"
                out.append({"tic": rec["tic"], "kind": kind, **s, "flags": bad, "hard_flags": bad})
                continue
            row = {"tic": rec["tic"], "kind": kind, **s, **d}
            # variability: same period found again after masking
            row["variability_repeat"] = any(same_period(s["period"], o["period"]) for o in sigs[i + 1:])
            fl = flags(row, rad)
            if row["variability_repeat"]:
                fl.append("variability_repeat")
            row["flags"] = ";".join(fl)
            row["hard_flags"] = ";".join(x for x in fl if x not in SOFT)
            out.append(row)
    return out


def main():
    outdir, targets = sys.argv[1], sys.argv[2]
    tg = pd.read_csv(targets).set_index("tic")
    files = sorted(glob.glob(f"{outdir}/stars/*.json"))
    args = []
    for p in files:
        tic = int(p.split("/")[-1].split(".")[0])
        r = tg.loc[tic] if tic in tg.index else None
        args.append((p, float(r["rad"]) if r is not None else 0.6, float(r["mass"]) if r is not None else 0.6))
    rows = []
    with mp.get_context("fork").Pool() as pool:
        for r in pool.imap_unordered(triage_star, args, chunksize=4):
            rows.extend(r)
    df = pd.DataFrame(rows)
    if len(df):
        # common period: the same period (0.1%) in >= 3 other stars => instrument
        real = df[df.kind == "signals"]
        cp = []
        for _, r in df.iterrows():
            others = real[(real.tic != r.tic) & (np.abs(real.period / r.period - 1) < 1e-3)]
            cp.append(others.tic.nunique() >= 3)
        df["common_period"] = cp
        df.loc[df.common_period, "hard_flags"] = df.loc[df.common_period, "hard_flags"].astype(str) + ";common_period"
    df.to_csv(f"{outdir}/triage.csv", index=False)
    for kind in ("signals", "signals_inverted"):
        d = df[df.kind == kind] if len(df) else df
        if len(d):
            clean = d[d["hard_flags"].fillna("") == ""]
            print(kind, "passing search:", len(d), " no hard flags:", len(clean),
                  " no flags:", int((d["flags"].fillna("") == "").sum()),
                  " max SDE no-hard:", round(clean.sde.max(), 2) if len(clean) else None)


if __name__ == "__main__":
    main()
