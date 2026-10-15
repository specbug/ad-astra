"""Sector 108 (our own full-frame photometry) against the period choices still allowed for HD 11506 and TOI-2493.01.

For every allowed period, each transit it predicts inside sector 108 is measured: a box of the transit's duration
against a quadratic baseline fitted to the 0.6-day stretch around it (outside the box), with the box depth's error
from the out-of-box scatter scaled by the red-noise factor. A period is ruled out when a predicted transit is at
least 60% covered and the measured depth is more than 4 sigma below the known depth (and below half of it).

  python s108_alias.py
Output: results/candidates/s108_alias.json
"""
import json
import math

import env  # noqa: F401
import numpy as np

from env import CATALOGS, RESULTS

OUT = RESULTS / "candidates" / "s108_alias.json"


def load(tic):
    z = np.load(CATALOGS.parent / "cache" / "ffi_lc" / f"{tic}_s108.npz", allow_pickle=True)
    t, f = z["time"], z["flux"]
    ok = np.isfinite(t) & np.isfinite(f)
    return t[ok], f[ok] / np.nanmedian(f[ok])


def box(t, f, tc, D):
    w = np.abs(t - tc) < 0.3 + D / 2
    inn = w & (np.abs(t - tc) < 0.4 * D)
    oot = w & (np.abs(t - tc) > 0.6 * D)
    cad = float(np.median(np.diff(t)))
    cover = inn.sum() * cad / (0.8 * D)
    if inn.sum() < 3 or oot.sum() < 20:
        return {"tc": tc, "coverage": float(min(cover, 1)), "measured": False}
    c = np.polyfit(t[oot] - tc, f[oot], 2)
    r = f / np.polyval(c, t - tc)
    sig = 1.4826 * np.median(np.abs(np.diff(r[oot]))) / math.sqrt(2)
    n = max(int(round(D * 0.8 / cad)), 1)
    k = oot.sum() // n
    beta = 1.0
    if k >= 4:
        rb = (r[oot][: k * n] - 1).reshape(k, n).mean(axis=1)
        beta = max(1.0, float(np.std(rb) / (sig / math.sqrt(n))))
    d = float(1 - r[inn].mean())
    e = float(sig / math.sqrt(inn.sum()) * beta)
    return {"tc": tc, "coverage": float(min(cover, 1)), "measured": True, "depth_ppm": d * 1e6, "err_ppm": e * 1e6,
            "beta": beta}


def main():
    ex = json.load(open(RESULTS / "report" / "extras.json"))
    fit = json.load(open(RESULTS / "report" / "transit_fits.json"))
    out = {}
    # HD 11506: periods gap/n still allowed by sectors 3, 30, 107
    t, f = load(266680498)
    tc0 = fit["HD11506"]["tc_0"][0]
    D = fit["HD11506"]["t14_h"][0] / 24
    dep = 1e6 * fit["HD11506"]["rp_rs"][0] ** 2
    rows = []
    for a in ex["HD11506"]["aliases"]:
        if a["status"] != "allowed":
            continue
        P = a["P"]
        k = np.arange(math.ceil((t.min() - tc0) / P), math.floor((t.max() - tc0) / P) + 1)
        ev = [box(t, f, tc0 + i * P, D) for i in k]
        bad = [e for e in ev if e["measured"] and e["coverage"] >= 0.6 and e["depth_ppm"] < dep - 4 * e["err_ppm"]
               and e["depth_ppm"] < dep / 2]
        seen = [e for e in ev if e["measured"] and e["coverage"] >= 0.6 and e["depth_ppm"] > 3 * e["err_ppm"]]
        rows.append({"n": a["n"], "P": P, "events": ev, "ruled_out_s108": bool(bad), "transit_seen_s108": bool(seen)})
    surv = [r for r in rows if not r["ruled_out_s108"]]
    out["HD11506"] = {"depth_ppm": dep, "dur_h": D * 24, "n_allowed_before": len(rows), "n_allowed_after": len(surv),
                      "shortest_after": min(r["P"] for r in surv) if surv else None,
                      "seen_in_s108": [r["P"] for r in rows if r["transit_seen_s108"]], "rows": rows}
    print("HD 11506:", len(rows), "->", len(surv), "allowed; shortest", out["HD11506"]["shortest_after"],
          "; S108 transit seen for P =", out["HD11506"]["seen_in_s108"][:10])
    # TOI-2493.01: P, P/2, P/3 in sector 108
    t, f = load(123664207)
    f2 = fit["TOI2493"]
    P, T0 = f2["period_d"][0], f2["tref_btjd"][0]
    D = f2["t14_h"][0] / 24
    dep = 1e6 * f2["rp_rs"][0] ** 2
    res = {}
    for name, Pa in (("P", P), ("P/2", P / 2), ("P/3", P / 3), ("P/4", P / 4)):
        k = np.arange(math.ceil((t.min() - T0) / Pa), math.floor((t.max() - T0) / Pa) + 1)
        res[name] = {"P": Pa, "events": [box(t, f, T0 + i * Pa, D) for i in k]}
        print("TOI-2493", name, [(round(e["tc"], 2), round(e["coverage"], 2), round(e.get("depth_ppm", float("nan"))),
                                  round(e.get("err_ppm", float("nan")))) for e in res[name]["events"]])
    out["TOI2493"] = {"depth_ppm": dep, "dur_h": D * 24, "aliases": res}
    json.dump(out, open(OUT, "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
