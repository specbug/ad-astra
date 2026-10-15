"""Per-signal review of the period hunt's new periods (report.period_hunt_summary()["new"]), before any of them is
uploaded as an additional parameter set (independent audits, 15 October 2026: "no blanket certification").

For each period candidate (every alias still allowed):
  - every transit that period predicts inside the data is timed with a trapezoid fit (toi_audit_timing.time_one,
    at the candidate's duration); a predicted transit with good coverage measured below a third of the reference
    depth and more than 4 sigma under it counts as missing;
  - a weighted linear ephemeris with full covariance through the timed transits (errors scaled by sqrt(chi2_red));
  - odd and even depths compared;
  - the pixel verdict (period_vet.json) and the stellar-density comparison (period_deep.json) are carried over.
Tier: "ready" = one allowed period, at least 3 timed transits, none missing, chi2_red < 5, odd/even < 3 sigma, and
the pixels on the target (or unresolved from it); "aliases" = 2-3 periods allowed and each passes the light-curve
checks (upload, if at all, with the aliases in the note); otherwise "hold", with the reasons.

  N_PROC=4 python period_signal_review.py
Output: results/candidates/period_signal_review.json and .csv
"""
import ast
import json
import math
import os
from concurrent.futures import ProcessPoolExecutor, as_completed

import env  # noqa: F401
import numpy as np
import pandas as pd

from env import RESULTS
from toi_audit_timing import fit_line, time_one

C = RESULTS / "candidates"
OUT = C / "period_signal_review.json"


def lit(x):
    if isinstance(x, str):
        try:
            return ast.literal_eval(x)
        except Exception:
            return x
    return x


def box(lc, tc, D):
    t, f = lc["time"], lc["flux"]
    w = np.isfinite(f) & (np.abs(t - tc) < D / 2 + 0.4)
    inn = w & (np.abs(t - tc) < 0.4 * D)
    oot = w & (np.abs(t - tc) > 0.6 * D)
    if inn.sum() < 3 or oot.sum() < 15:
        return None
    cad = float(np.median(np.diff(t)))
    c = np.polyfit(t[oot] - tc, f[oot], 2)
    r = f / np.polyval(c, t - tc)
    sig = 1.4826 * np.median(np.abs(np.diff(r[oot]))) / math.sqrt(2)
    return {"coverage": min(1.0, inn.sum() * cad / (0.8 * D)), "depth": float(1 - r[inn].mean()), "err": float(sig / math.sqrt(inn.sum()))}


def check_period(lcs, P, T, D, depth_ref):
    ev, missing = [], []
    for lc in lcs:
        t = lc["time"]
        for k in range(math.ceil((t.min() - T) / P), math.floor((t.max() - T) / P) + 1):
            tc = T + k * P
            m = time_one(lc, tc, D)
            if m and not m.get("edge") and m["snr"] >= 3 and abs(m["tc"] - tc) < 0.5 * D:
                ev.append({"n": k, "sector": lc["meta"].get("sector"), "tc": m["tc"], "err": m["err"], "depth": m["depth"],
                           "depth_err": m["depth_err"], "snr": m["snr"]})
                continue
            b = box(lc, tc, D)
            if b and b["coverage"] >= 0.6 and b["depth"] < depth_ref / 3 and (depth_ref - b["depth"]) > 4 * b["err"]:
                missing.append({"n": k, "sector": lc["meta"].get("sector"), "tc": tc, "depth_ppm": b["depth"] * 1e6,
                                "err_ppm": b["err"] * 1e6})
    out = {"P_trial": P, "n_timed": len(ev), "n_missing": len(missing), "missing": missing[:10], "events": ev}
    if len(ev) >= 3:
        n = np.array([e["n"] for e in ev], float)
        n0, p, cov, chi2, dof, res = fit_line(n, np.array([e["tc"] for e in ev]), np.array([e["err"] for e in ev]))
        out.update(P=float(p[1]), P_err=float(math.sqrt(cov[1, 1])), T0=float(p[0]), T0_err=float(math.sqrt(cov[0, 0])),
                   cov_T0_P=float(cov[0, 1]), chi2_red=chi2 / dof, o_c_rms_min=float(np.sqrt(np.mean(res ** 2)) * 1440))
        dep = np.array([e["depth"] for e in ev])
        de = np.array([e["depth_err"] for e in ev])
        w = 1 / de ** 2
        out["depth_ppm"] = float(np.sum(w * dep) / w.sum() * 1e6)
        out["depth_err_ppm"] = float(1 / math.sqrt(w.sum()) * 1e6)
        odd = (n.astype(int) % 2) == 1
        if odd.sum() >= 1 and (~odd).sum() >= 1:
            wo, we = w[odd], w[~odd]
            out["odd_even_sigma"] = float(abs(np.sum(wo * dep[odd]) / wo.sum() - np.sum(we * dep[~odd]) / we.sum())
                                          / math.hypot(1 / math.sqrt(wo.sum()), 1 / math.sqrt(we.sum())))
    return out


def review(o, pix, dens):
    from lcaccess import lcs_for
    tic = int(o["tic"])
    lcs = [lc for lc in lcs_for(tic) if lc["time"].size > 100]
    D = float(o["dur_h"]) / 24
    T = float(o["ref_epoch"])
    dref = float(o["depth_ref_ppm"]) * 1e-6
    aliases = [float(x) for x in lit(o.get("aliases")) or []]
    key = f"{tic}|{o['id']}"
    res = {"key": key, "tic": tic, "id": o["id"], "disp": o.get("disp"), "dur_h": float(o["dur_h"]),
           "depth_ref_ppm": float(o["depth_ref_ppm"]), "aliases": aliases, "pixels": pix.get(key, {}).get("verdict"),
           "pixel_offset_sigma": pix.get(key, {}).get("offset_sigma"), "density": dens.get(key),
           "trials": [check_period(lcs, P, T, D, dref) for P in aliases]}
    ok = [t for t in res["trials"] if t["n_timed"] >= 3 and t["n_missing"] == 0 and t.get("chi2_red", 99) < 5
          and (t.get("odd_even_sigma") or 0) < 3]
    pix_ok = bool(res["pixels"]) and (res["pixels"] == "on target" or str(res["pixels"]).startswith("target or"))
    why = []
    if not pix_ok:
        why.append(f"pixels: {res['pixels']}")
    for t in res["trials"]:
        tw = []
        if t["n_timed"] < 3:
            tw.append(f"{t['n_timed']} transits timed")
        if t["n_missing"]:
            tw.append(f"{t['n_missing']} predicted transits missing")
        if t.get("chi2_red", 0) >= 5:
            tw.append(f"timing chi2_red {t['chi2_red']:.1f}")
        if (t.get("odd_even_sigma") or 0) >= 3:
            tw.append(f"odd/even {t['odd_even_sigma']:.1f} sigma")
        if tw:
            why.append(f"P={t['P_trial']:.4f}: " + ", ".join(tw))
    if pix_ok and len(aliases) == 1 and len(ok) == 1:
        res["tier"] = "ready"
    elif pix_ok and 2 <= len(aliases) <= 3 and len(ok) == len(aliases):
        res["tier"] = "aliases"
    else:
        res["tier"] = "hold"
    res["why"] = "; ".join(why)
    res["n_alias_pass"] = len(ok)
    return res


def main():
    from report import period_hunt_summary
    S = period_hunt_summary()
    new = S["new"]
    pix = json.load(open(C / "period_vet.json"))
    pix = {k: {a: lit(b) for a, b in v.items()} for k, v in pix.items()}
    dd = json.load(open(C / "period_deep.json")) if (C / "period_deep.json").exists() else {}
    dens = {}
    for k, v in (dd.items() if isinstance(dd, dict) else []):
        if isinstance(v, dict) and v.get("fits"):
            dens[k] = [{"P": f.get("P"), "rho_solar": f.get("rho_solar"), "rho_z": f.get("rho_z")} for f in v["fits"]]
    res = json.load(open(OUT)) if OUT.exists() else {}
    todo = [o for o in new if f"{o['tic']}|{o['id']}" not in res]
    print(len(new), "new periods;", len(todo), "to review", flush=True)
    with ProcessPoolExecutor(int(os.environ.get("N_PROC", "4"))) as ex:
        futs = {ex.submit(review, o, pix, dens): f"{o['tic']}|{o['id']}" for o in todo}
        for i, fu in enumerate(as_completed(futs)):
            k = futs[fu]
            try:
                res[k] = fu.result()
            except Exception as e:
                res[k] = {"key": k, "tier": "error", "error": repr(e)[:200]}
            if i % 10 == 0:
                json.dump(res, open(OUT, "w"), indent=1, default=float)
                print(i, "/", len(todo), flush=True)
    json.dump(res, open(OUT, "w"), indent=1, default=float)
    rows = [{"key": r["key"], "tier": r.get("tier"), "aliases": r.get("aliases"), "pixels": r.get("pixels"), "why": r.get("why", r.get("error"))}
            for r in res.values()]
    pd.DataFrame(rows).to_csv(C / "period_signal_review.csv", index=False)
    print(pd.Series([r.get("tier") for r in res.values()]).value_counts().to_string())


if __name__ == "__main__":
    main()
