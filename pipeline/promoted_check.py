"""Community TOIs merged into a TOI whose ephemeris does not pass through the community TOI's transit.

ExoFOP marks a community TOI as "promoted" when the TESS project makes it a TOI (or finds it is the same object
as one). For each promoted community TOI, the offset of its transit from the TOI's ephemeris is computed from the
catalogue values (ExoFOP tables of 12 Oct 2026), in hours and in units of the ephemeris error at that epoch.
Those more than one transit duration and 5 sigma off are then checked against the TESS light curves:

  - does the transit of the community TOI exist (a dip at its epoch)?
  - does it match the ephemeris of another TOI on the same star (a label swap)?
  - do the transits the TOI's ephemeris predicts show up in the data?
  - if not, which period near the TOI's, through the community TOI's transit, does fit?

  python promoted_check.py
Output: results/candidates/promoted_check.json (and the flagged list, promoted_offsets.csv)
"""
import json
import math

import env  # noqa: F401
import numpy as np
import pandas as pd

from env import CATALOGS, RESULTS
from lcaccess import lcs_for
from period_review import lc_for, measure

C = RESULTS / "candidates"
E = CATALOGS / "exofop"
OUT = C / "promoted_check.json"
BJD0 = 2457000.0


def flagged():
    toi = pd.read_csv(E / "toi_20261012.csv", low_memory=False)
    ctoi = pd.read_csv(E / "ctoi_20261012.csv", low_memory=False)
    toi["key"] = toi["TOI"].map(lambda x: f"{float(x):.2f}")
    c = ctoi[ctoi["Promoted to TOI"].notna()].copy()
    c["key"] = c["Promoted to TOI"].map(lambda x: f"{float(x):.2f}" if str(x).replace(".", "", 1).isdigit() else str(x))
    m = c.merge(toi, on="key", how="inner", suffixes=("_c", "_t"))
    P = pd.to_numeric(m["Period (days)_t"], errors="coerce")
    Pe = pd.to_numeric(m["Period (days) err"], errors="coerce").fillna(0)
    T0 = pd.to_numeric(m["Epoch (BJD)"], errors="coerce")
    T0e = pd.to_numeric(m["Epoch (BJD) err"], errors="coerce").fillna(0)
    tc = pd.to_numeric(m["Transit Epoch (BJD)"], errors="coerce")
    D = pd.to_numeric(m["Duration (hours)"], errors="coerce") / 24
    k = ((tc - T0) / P).round()
    off = tc - (T0 + k * P)
    sig = np.sqrt((k.abs() * Pe) ** 2 + T0e ** 2)
    m = m.assign(P=P, T0=T0 - BJD0, tc=tc - BJD0, D=D, off_h=off * 24, sig_h=sig * 24,
                 off_sig=off.abs() / sig.replace(0, np.nan), off_dur=off.abs() / D, k=k)
    ok = (m.P > 0) & m.tc.notna() & m.T0.notna() & (m.tc > 1000) & m.D.notna()
    m = m[ok]
    bad = m[(m.off_dur > 1.0) & ((m.off_sig > 5) | m.off_sig.isna())].copy()
    bad.to_csv(C / "promoted_offsets.csv", index=False)
    return bad, toi, len(m)


def predicted(lcs, T0, P, D):
    out = []
    for lc in lcs:
        t = lc["time"]
        lo, hi = float(np.nanmin(t)), float(np.nanmax(t))
        if P > 0:
            out += [(n, T0 + n * P) for n in range(math.ceil((lo + D / 2 - T0) / P), math.floor((hi - D / 2 - T0) / P) + 1)]
    return sorted(set(out))


def at(lcs, tc, D):
    lc = lc_for(lcs, tc)
    m = measure(lc, tc, D) if lc is not None else None
    if m is None or m["coverage"] < 0.5:
        return None
    return {"tc": round(tc, 4), "sector": lc["meta"].get("sector"), "depth_ppm": round(m["depth"] * 1e6, 1),
            "err_ppm": round(abs(m["depth"] / m["snr"]) * 1e6, 1) if m["snr"] else None, "snr": round(m["snr"], 2)}


def series(lcs, T0, P, D, dref):
    ev = []
    for n, tc in predicted(lcs, T0, P, D):
        x = at(lcs, tc, D)
        if x and x["err_ppm"]:
            x["n"] = n
            x["dip"] = x["snr"] >= 3 and x["depth_ppm"] > 0.4 * dref
            x["no_dip"] = x["depth_ppm"] < dref / 3 and (dref - x["depth_ppm"]) / x["err_ppm"] > 3
            ev.append(x)
    return ev


def best_near(lcs, T0, P0, D, dref, frac=5e-3):
    """Period within P0 (1 +- frac) through T0 whose predicted transits show the most signal."""
    span = max(abs(float(np.nanmax(lc["time"])) - T0) for lc in lcs) + max(abs(float(np.nanmin(lc["time"])) - T0) for lc in lcs)
    step = max(P0 * (D / 4) / max(span, P0), 1e-6)
    grid = np.arange(P0 * (1 - frac), P0 * (1 + frac) + step, step)
    if grid.size > 3000:
        grid = np.linspace(P0 * (1 - frac), P0 * (1 + frac), 3000)
    best = (-np.inf, P0, [])
    for P in grid:
        ev = series(lcs, T0, P, D, dref)
        z = [e["depth_ppm"] / e["err_ppm"] for e in ev]
        s = float(np.sum(z) / math.sqrt(len(z))) if z else -np.inf
        if s > best[0]:
            best = (s, float(P), ev)
    return best


def summary(ev):
    return {"n_data": len(ev), "n_dip": sum(e["dip"] for e in ev), "n_no_dip": sum(e["no_dip"] for e in ev)}


def main():
    bad, toi, n_all = flagged()
    print(len(bad), "of", n_all, "promoted community TOIs are off their TOI's ephemeris", flush=True)
    res = json.load(open(OUT)) if OUT.exists() else {}
    for _, r in bad.iterrows():
        tic = int(r["TIC ID_c"])
        key = f"{tic}|{r['CTOI']}|{r['key']}"
        if key in res:
            continue
        out = {"tic": tic, "ctoi": str(r["CTOI"]), "toi": r["key"], "disp": str(r["TFOPWG Disposition_t"]),
               "user": str(r["User"]), "toi_source": str(r.get("Detection")), "P": float(r.P), "T0": float(r.T0),
               "dur_h": float(r.D * 24), "toi_depth_ppm": float(r["Depth (ppm)"]) if r["Depth (ppm)"] == r["Depth (ppm)"] else None,
               "ctoi_tc": float(r.tc), "ctoi_depth": r["Depth ppm"] if r["Depth ppm"] == r["Depth ppm"] else None,
               "off_h": float(r.off_h), "off_sig": float(r.off_sig) if r.off_sig == r.off_sig else None, "off_dur": float(r.off_dur)}
        # does the transit match another TOI on the same star?
        others = toi[(toi["TIC ID"] == tic) & (toi["key"] != r["key"])]
        out["matches_other_toi"] = None
        for _, o in others.iterrows():
            P, T0 = float(o["Period (days)"]), float(o["Epoch (BJD)"]) - BJD0
            Dh = float(o["Duration (hours)"]) if o["Duration (hours)"] == o["Duration (hours)"] else 3.0
            if not P > 0:
                if abs(r.tc - T0) < max(Dh / 48, 0.1):
                    out["matches_other_toi"] = o["key"]
                continue
            k = round((r.tc - T0) / P)
            sig = math.hypot(abs(k) * (float(o["Period (days) err"]) if o["Period (days) err"] == o["Period (days) err"] else 0),
                             float(o["Epoch (BJD) err"]) if o["Epoch (BJD) err"] == o["Epoch (BJD) err"] else 0)
            if abs(r.tc - (T0 + k * P)) < max(Dh / 48, 3 * sig):
                out["matches_other_toi"] = o["key"]
        try:
            lcs = [lc for lc in lcs_for(tic) if lc["time"].size > 100]
        except Exception as ex:
            out["error"] = repr(ex)[:200]
            res[key] = out
            continue
        if not lcs:
            out["error"] = "no light curves"
            res[key] = out
            json.dump(res, open(OUT, "w"), indent=1, default=float)
            continue
        D = float(r.D)
        dref = out["toi_depth_ppm"] or 1000.0
        out["ctoi_transit"] = at(lcs, float(r.tc), D)
        ev = series(lcs, float(r.T0), float(r.P), D, dref)
        out["toi_series"] = summary(ev)
        out["toi_events"] = ev
        if out["matches_other_toi"] is None and 1.0 <= float(r.P) < 200:
            s, Pb, evb = best_near(lcs, float(r.tc), float(r.P), D, dref)
            out["through_ctoi"] = {"P": Pb, "snr": s, **summary(evb)}
            out["through_ctoi_events"] = evb
        print(key, out["disp"], f"P={r.P:.4f}", f"off={r.off_h:.1f} h", "other:", out["matches_other_toi"],
              "| ctoi dip:", (out["ctoi_transit"] or {}).get("snr"), "| toi series:", out["toi_series"],
              "| through ctoi:", {k: out.get("through_ctoi", {}).get(k) for k in ("P", "n_data", "n_dip", "n_no_dip")}, flush=True)
        res[key] = out
        json.dump(res, open(OUT, "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
