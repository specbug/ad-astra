"""Review of the period hunt (period_hunt.py): which period-less TOIs and CTOIs gained credible new transits,
what period they give, and whether that period is already on record.

For every target with new events:
  1. each new event, and the catalogued one, is re-measured on the light curve of its own sector (not the
     flattened one used by the hunt): a box of the hunt's duration against a local quadratic baseline fitted
     within 3 durations (+0.3 d) of the event, the event masked. An event is credible when
       - its box SNR is >= 7 with the noise of that window (red noise included through the scatter of
         box means at other positions in the window);
       - >= 70% of the expected in-transit cadences exist and there are data on both sides;
       - it is not within 0.5 d of a sector start or end, or next to a data gap longer than 0.5 d;
       - the background (when the light curve has one) does not jump during the event (< 5 MAD);
       - the SAP flux (when present) shows at least half of the dip;
       - its depth is 0.5-2x that of the catalogued event.
  2. the credible events and the catalogued one go through duo.alias_table: every period that puts them all on
     one ephemeris and is not excluded by the rest of the data. The shortest period a circular orbit allows
     for the measured duration and the TIC stellar density is given too: aliases far below it are unlikely.
  3. on record: the NASA Exoplanet Archive (planets on the star, with periods) and the ExoFOP page of the
     star (periods uploaded by anyone).
Output: results/candidates/period_review.csv and period_review.json

  python period_review.py
"""
import json
import math
import os

import env  # noqa: F401
import numpy as np
import pandas as pd

from env import CATALOGS, RESULTS

C = RESULTS / "candidates"
G = 6.674e-11


def measure(lc, tc, D):
    """Box depth of an event against a local quadratic baseline, with quality diagnostics."""
    t, f = lc["time"], lc["flux"]
    w = np.abs(t - tc) < 3 * D + 0.3
    if w.sum() < 20:
        return None
    tt, ff = t[w], f[w]
    inn = np.abs(tt - tc) < 0.4 * D
    out = np.abs(tt - tc) > 0.6 * D
    if inn.sum() < 3 or out.sum() < 10:
        return None
    c = np.polyfit(tt[out] - tc, ff[out], 2)
    r = ff / np.polyval(c, tt - tc)
    s = 1.4826 * np.median(np.abs(r[out] - np.median(r[out])))
    depth = float(1 - np.mean(r[inn]))
    n_in = int(inn.sum())
    # red noise: scatter of box means (as many points as in transit) at other places in the window
    ro = r[out]
    if ro.size >= 3 * n_in:
        k = ro.size // n_in
        means = ro[: k * n_in].reshape(k, n_in).mean(axis=1)
        s_box = max(1.4826 * np.median(np.abs(means - np.median(means))), s / math.sqrt(n_in))
    else:
        s_box = s / math.sqrt(n_in)
    cad = np.median(np.diff(t))
    out_d = {"tc": float(tc), "depth": depth, "snr": depth / s_box if s_box > 0 else 0.0, "n_in": n_in,
             "coverage": float(n_in / max(0.8 * D / cad, 1)),
             "left": int((tt < tc - 0.6 * D).sum()), "right": int((tt > tc + 0.6 * D).sum()),
             "edge_days": float(min(tc - t.min(), t.max() - tc))}
    near = t[np.abs(t - tc) < D / 2 + 0.5]
    out_d["max_gap_near"] = float(np.max(np.diff(near))) if near.size > 1 else 99.0
    for key in ("bkg", "cx", "cy"):
        v = lc.get(key)
        if v is None:
            continue
        v = np.asarray(v)[w]
        vo, vi = v[out], v[inn]
        if not (np.isfinite(vo).sum() > 5 and np.isfinite(vi).sum() > 1):
            continue
        mad = 1.4826 * np.nanmedian(np.abs(vo - np.nanmedian(vo))) + 1e-12
        if key == "bkg":
            out_d["bkg_jump_mad"] = float((np.nanmedian(vi) - np.nanmedian(vo)) / mad)
        else:
            out_d[f"{key}_shift"] = float(np.nanmean(vi) - np.nanmean(vo))
            out_d[f"{key}_shift_err"] = float(mad / math.sqrt(np.isfinite(vi).sum()))
    sap = lc.get("sap_flux")
    if sap is not None and lc["meta"].get("product") != "ffi":
        sp = np.asarray(sap)[w]
        ok = np.isfinite(sp)
        if ok[out].sum() > 10 and ok[inn].sum() > 2:
            cs = np.polyfit((tt - tc)[out & ok], sp[out & ok], 2)
            rs = sp / np.polyval(cs, tt - tc)
            out_d["sap_depth"] = float(1 - np.nanmean(rs[inn]))
    return out_d


def credible(m, d_ref):
    why = []
    if m is None:
        return False, ["not measurable"]
    if m["snr"] < 7:
        why.append(f"SNR {m['snr']:.1f}")
    if m["coverage"] < 0.7 or m["left"] < 5 or m["right"] < 5:
        why.append("partial coverage")
    if m["edge_days"] < 0.5 or m["max_gap_near"] > 0.5:
        why.append("at a sector edge or data gap")
    if abs(m.get("bkg_jump_mad", 0)) > 5:
        why.append(f"background jump {m['bkg_jump_mad']:.0f} MAD")
    # (centroid moves are recorded but not used here: an on-target transit also moves the flux-weighted
    # centroid when neighbours share the aperture, and differently in every sector; the shortlisted events get
    # a difference-image centroid test instead)
    if "sap_depth" in m and m["sap_depth"] < 0.5 * m["depth"]:
        why.append("weak in SAP flux")
    if d_ref and not (0.5 * d_ref <= m["depth"] <= 2 * d_ref):
        why.append(f"depth {m['depth'] * 1e6:.0f} vs {d_ref * 1e6:.0f} ppm")
    return not why, why


def p_min_circular(dur_h, rho_sun, k=0.03):
    """Shortest period for which a central transit of a circular orbit lasts dur_h (stellar density in solar units)."""
    rho = rho_sun * 1408.0
    T = dur_h * 3600
    return (T * math.pi / (1 + k)) ** 3 * G * rho / (3 * math.pi) / 86400


def lc_for(lcs, tc):
    for lc in lcs:
        if lc["time"].min() - 0.1 <= tc <= lc["time"].max() + 0.1:
            return lc
    return None


def on_record(tic, toi_id, nasa, ra, dec):
    from exofop_check import fetch
    out = {}
    sep = np.hypot((nasa.ra.values - ra) * math.cos(math.radians(dec)), nasa.dec.values - dec) * 3600
    n = nasa[sep < 10]
    if len(n):
        out["nasa"] = "; ".join(f"{a} P={p:.4f}" for a, p in zip(n.pl_name, n.pl_orbper))
    d = fetch(tic)
    ups = []
    for e in d.get("planet_parameters", []) if isinstance(d, dict) else []:
        per = e.get("per")
        try:
            per = float(per)
        except (TypeError, ValueError):
            continue
        if per > 0:
            ups.append(f"{e.get('name') or e.get('toi') or ''} P={per:.4f} ({e.get('puser')}, {str(e.get('pdate'))[:10]})")
    out["exofop_periods"] = "; ".join(ups)
    return out


def review(r, nasa, ticinfo):
    from dossier import all_sectors, best_lcs
    from duo import alias_table, flatten_sectors
    from ffi_new import new_lcs
    tic = int(r["tic"])
    D = (r.get("dur_fit_h") or r.get("dur_h") or 6.0) / 24
    lcs = best_lcs(tic, all_sectors(tic))
    st = ticinfo.loc[tic]
    have = {lc["meta"]["sector"] for lc in lcs}
    lcs += [lc for lc in new_lcs(tic, float(st.ra), float(st.dec), build_missing=False) if lc["meta"]["sector"] not in have]
    ref = r.get("ref_epoch") or r["epoch"]
    m_ref = measure(lc_for(lcs, ref), ref, D) if lc_for(lcs, ref) else None
    d_ref = m_ref["depth"] if m_ref else (r.get("depth_fit_ppm") or 0) * 1e-6
    evs = []
    for o in r.get("others", []):
        lc = lc_for(lcs, o["t"])
        m = measure(lc, o["t"], D) if lc else None
        ok, why = credible(m, d_ref)
        evs.append({"t": o["t"], "sector": lc["meta"]["sector"] if lc else None,
                    "product": lc["meta"].get("product") if lc else None, "ok": ok, "why": why, **(m or {})})
    good = sorted([ref] + [e["t"] for e in evs if e["ok"]])
    out = {"tic": tic, "id": r["id"], "disp": r.get("disp"), "ref_epoch": ref, "dur_h": D * 24,
           "depth_ref_ppm": d_ref * 1e6, "ref": m_ref, "events": evs, "n_credible": len(good) - 1}
    rho = st.get("rho")
    if rho == rho and rho:
        out["p_min_circular"] = p_min_circular(D * 24, float(rho), k=math.sqrt(max(d_ref, 1e-6)))
    if len(good) >= 2 and d_ref > 0:
        flat = flatten_sectors(lcs, window=max(2.5, 4 * D))
        rows = alias_table(good, 0.6 * d_ref, D * 24, flat, p_min=1.0)
        allowed = [a["P"] for a in rows if a["status"] == "allowed"]
        out["aliases"] = [round(p, 5) for p in allowed]
        out["n_aliases"] = len(allowed)
        out["dt"] = good[-1] - good[0]
    out.update(on_record(tic, r["id"], nasa, float(st.ra), float(st.dec)))
    return out


def _review_safe(r, nasa, ti):
    try:
        return review(r, nasa, ti)
    except Exception as ex:
        return {"tic": r["tic"], "id": r["id"], "error": repr(ex)[:200]}


def main():
    import ticdb
    res = []
    for name in ("period_hunt", "period_hunt_ffi", "period_hunt_singles"):
        p, ck = C / f"{name}.json", CATALOGS.parent / f"{name}.jsonl"
        if p.exists():
            res += json.load(open(p))
        elif ck.exists():   # a run still in progress: its checkpoint
            res += [json.loads(line) for line in open(ck) if line.strip()]
    # the FFI pass supersedes the first pass for the same target
    best = {}
    for r in res:
        key = (r["tic"], r["id"])
        if key not in best or r.get("ffi_sectors"):
            best[key] = r
    todo = [r for r in best.values() if r.get("status") == "new events" and r.get("recovered")]
    print(len(best), "targets;", len(todo), "with new events and the catalogued event re-found", flush=True)
    nasa = pd.read_csv(CATALOGS / "nasa_pscomppars_transiting_20261012.csv", low_memory=False)
    tics = sorted({int(r["tic"]) for r in todo})
    ti = ticdb.query(tics, ["ra", "dec", "Tmag", "rho", "rad", "mass"]).set_index("tic")
    import concurrent.futures as cf
    # checkpoint: one JSON line per reviewed target, so an interrupted run resumes where it stopped
    ck = CATALOGS.parent / "period_review.jsonl"
    out = []
    if ck.exists():
        for line in open(ck):
            try:
                out.append(json.loads(line))
            except ValueError:
                pass
    done = {(o["tic"], o["id"]) for o in out}
    todo = [r for r in todo if (r["tic"], r["id"]) not in done]
    print(len(done), "already reviewed;", len(todo), "to go", flush=True)
    with cf.ProcessPoolExecutor(int(os.environ.get("N_PROC", "4"))) as ex, open(ck, "a") as fh:
        futs = {ex.submit(_review_safe, r, nasa, ti): r for r in todo}
        for i, fu in enumerate(cf.as_completed(futs)):
            o = fu.result()
            out.append(o)
            fh.write(json.dumps(o, default=float) + "\n")
            fh.flush()
            if i % 20 == 0:
                print(i, o["tic"], o["id"], o.get("n_credible"), o.get("n_aliases"), flush=True)
    json.dump(out, open(C / "period_review.json", "w"), indent=1, default=float)
    flat = pd.DataFrame([{k: (v if not isinstance(v, (list, dict)) else json.dumps(v, default=float)[:300])
                          for k, v in o.items() if k not in ("events", "ref")} for o in out])
    flat.to_csv(C / "period_review.csv", index=False)
    print("credible new events:", int((flat.n_credible > 0).sum()), "targets")


if __name__ == "__main__":
    main()
