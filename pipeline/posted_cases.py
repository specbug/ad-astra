"""Case studies for the community TOIs whose single transit does not belong to the period posted for it.

posted_check.py found five candidates where a period posted to ExoFOP for a community TOI (directly, or through
the TOI it was merged into) does not pass through the community TOI's own transit, or is twice the true period.
For each, the competing ephemerides are measured transit by transit (period_review.measure), the light curves
around every predicted transit are drawn, and the transits are tested at the pixel level (centroid_vet.vet).

  python posted_cases.py [tic ...]
Output: results/candidates/posted_cases.json and one figure per star in results/report/cases/
"""
import json
import math
import sys

import env  # noqa: F401
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from centroid_vet import vet
from env import RESULTS
from lcaccess import lcs_for
from period_review import lc_for, measure

C = RESULTS / "candidates"
OUT = C / "posted_cases.json"
FIG = RESULTS / "report" / "cases"
BTJD = 2457000.0

# ephemerides: (name, P [d], T0 [BTJD], duration [h]); P None = list of single events (times in T0)
CASES = {
    438225366: {"label": "CTOI 438225366.01", "eph": [
        ("posted (41.72 d)", 41.7219, 2221.5204, 5.0),
        ("half of posted", 20.86132, 2221.5204, 5.0),
        ("phase 0.5 of 20.86 d", 20.86132, 2221.5204 + 20.86132 / 2, 5.0)]},
    156716001: {"label": "CTOI 156716001.01 / TOI-7602.01", "eph": [
        ("TOI-7602.01", 10.29778812572357, 2460700.427709 - BTJD, 3.80),
        ("TOI-7602.01 through the CTOI transit", 10.29778812572357, 2458485.406 - BTJD, 3.80),
        ("CTOI events", None, [2458485.406 - BTJD, 2224.52, 3702.77], 8.0)]},
    67711670: {"label": "CTOI 67711670.01 / TOI-6033.01", "eph": [
        ("TOI-6033.01", 18.0401494432717, 2459885.981072 - BTJD, 6.04),
        ("CTOI events", None, [2458798.9854 - BTJD, 3622.80], 6.0)]},
    369532319: {"label": "CTOI 369532319.01 / TOI-5967.01", "eph": [
        ("TOI-5967.01", 15.343585, 2460563.656882 - BTJD, 3.90),
        ("CTOI events", None, [2458755.81 - BTJD, 2831.56], 5.6)]},
    386131764: {"label": "CTOI 386131764.01 / TOI-2633.01", "eph": [
        ("TOI-2633.01", 18.1960490983542, 2459229.361575 - BTJD, 4.86),
        ("CTOI series (18.215 d)", 18.2149, 2459253.5816 - BTJD, 4.0)]},
}


def predicted(lcs, T0, P, D):
    out = []
    for lc in lcs:
        t = lc["time"]
        lo, hi = float(np.nanmin(t)), float(np.nanmax(t))
        out += [(n, T0 + n * P) for n in range(math.ceil((lo + D / 2 - T0) / P), math.floor((hi - D / 2 - T0) / P) + 1)]
    return sorted(set(out))


def measure_all(lcs, eph):
    name, P, T0, dur_h = eph
    D = dur_h / 24
    pts = [(i, t) for i, t in enumerate(T0)] if P is None else predicted(lcs, T0, P, D)
    ev = []
    for n, tc in pts:
        lc = lc_for(lcs, tc)
        m = measure(lc, tc, D) if lc is not None else None
        if m is None or m["coverage"] < 0.5:
            continue
        ev.append({"n": n, "tc": round(tc, 4), "sector": lc["meta"].get("sector"), "product": lc["meta"].get("product"),
                   "depth_ppm": round(m["depth"] * 1e6, 1), "err_ppm": round(abs(m["depth"] / m["snr"]) * 1e6, 1) if m["snr"] else None,
                   "snr": round(m["snr"], 2)})
    return ev


def plot(tic, case, lcs, results):
    FIG.mkdir(parents=True, exist_ok=True)
    marks = []
    for (name, P, T0, dur_h), ev in zip(case["eph"], results):
        for e in ev:
            marks.append((e["tc"], name, dur_h))
    marks.sort()
    if not marks:
        return None
    n = len(marks)
    cols = min(4, n)
    rows = math.ceil(n / cols)
    fig, axs = plt.subplots(rows, cols, figsize=(3.2 * cols, 2.4 * rows), squeeze=False)
    names = [e[0] for e in case["eph"]]
    colors = dict(zip(names, ["#c0392b", "#2471a3", "#7d3c98", "#1e8449"]))
    for ax, (tc, name, dur_h) in zip(axs.flat, marks):
        lc = lc_for(lcs, tc)
        t, f = lc["time"], lc["flux"]
        w = np.abs(t - tc) < 0.75
        if w.sum() < 10:
            ax.set_visible(False)
            continue
        tt, ff = t[w], f[w] / np.nanmedian(f[w])
        ax.plot((tt - tc) * 24, ff, ".", ms=1.5, color="#999999")
        b = np.arange(-18, 18.01, 0.5)
        idx = np.digitize((tt - tc) * 24, b)
        bx = [np.mean((tt[idx == k] - tc) * 24) for k in range(1, len(b)) if (idx == k).sum() > 2]
        by = [np.mean(ff[idx == k]) for k in range(1, len(b)) if (idx == k).sum() > 2]
        ax.plot(bx, by, "-", color=colors[name], lw=1.2)
        ax.axvspan(-dur_h / 2, dur_h / 2, color=colors[name], alpha=0.12)
        ax.set_title(f"{name[:26]}\nS{lc['meta'].get('sector')} {tc:.2f}", fontsize=7)
        ax.tick_params(labelsize=6)
    for ax in list(axs.flat)[n:]:
        ax.set_visible(False)
    fig.suptitle(f"TIC {tic} — {case['label']}", fontsize=9)
    fig.tight_layout()
    p = FIG / f"TIC{tic}_cases.png"
    fig.savefig(p, dpi=110)
    plt.close(fig)
    return str(p)


def main(tics):
    res = json.load(open(OUT)) if OUT.exists() else {}
    for tic in tics:
        case = CASES[tic]
        lcs = [lc for lc in lcs_for(tic) if lc["time"].size > 100]
        out = {"label": case["label"], "sectors": sorted({lc["meta"].get("sector") for lc in lcs}), "ephemerides": []}
        results = []
        for eph in case["eph"]:
            ev = measure_all(lcs, eph)
            results.append(ev)
            out["ephemerides"].append({"name": eph[0], "P": eph[1], "T0": eph[2], "dur_h": eph[3], "events": ev})
            print(tic, eph[0])
            for e in ev:
                print(f"    n={e['n']:4d} S{e['sector']} {e['product']:>6s} tc={e['tc']:.3f} depth={e['depth_ppm']:7.0f} +- {e['err_ppm']} snr={e['snr']}")
        out["figure"] = plot(tic, case, lcs, results)
        # pixel test of the strongest events of each ephemeris (four sectors at most)
        for eo in out["ephemerides"]:
            ev = sorted([e for e in eo["events"] if e["snr"] >= 5], key=lambda e: -e["snr"])
            seen, pick = set(), []
            for e in ev:
                if e["sector"] not in seen and len(seen) < 4:
                    seen.add(e["sector"])
                pick += [(e["sector"], e["tc"])] if e["sector"] in seen else []
            if len(pick) < 1:
                continue
            depth = float(np.median([e["depth_ppm"] for e in ev])) * 1e-6
            try:
                r = vet(tic, pick, eo["dur_h"], max(depth, 1e-4))
            except Exception as ex:
                r = {"verdict": "error", "error": repr(ex)[:200]}
            eo["pixels"] = {k: r.get(k) for k in ("verdict", "offset_arcsec", "n_used", "sectors") if k in r}
            print(tic, eo["name"], "pixels:", r.get("verdict"), r.get("offset_arcsec"), flush=True)
        res[str(tic)] = out
        json.dump(res, open(OUT, "w"), indent=1, default=float)


if __name__ == "__main__":
    main([int(a) for a in sys.argv[1:]] or list(CASES))
