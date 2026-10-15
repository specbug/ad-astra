"""Vetting of the Track A signals that beat the flipped-data baseline (tracka_v2.py) and are on no list.

For each signal:
  - pixels: per sector, the difference images (out of transit minus in transit, eb_pixels._diff_image) of every
    transit in the sector are stacked (inverse-variance weighted) and one PSF is fitted to the stack
    (eb_pixels.fit_diff, offsets measured from the target's own image); sectors are combined and compared with the
    target and the Gaia stars bright enough to make the dip (eb_pixels.combine_centroids, centroid_verdict). The
    signals are too shallow for single-transit centroids;
  - per-sector depths and the phase-folded light curve (figure), from the same stitched light curve the search used;
  - Gaia DR3 (gaia_check.match): RUWE, non-single-star flag, radial-velocity scatter;
  - ExoFOP target page, fetched today: TOIs, community TOIs and period uploads after the 12 October catalogue.

  python tracka_vet.py [TIC:period ...]          (default: every signal in tracka_v2.json "unknown_above")
Output: results/track_a/tracka_vet.json, results/report/cases/tracka_<tic>_<P>.png
"""
import json
import math
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed

import env  # noqa: F401
import numpy as np

from env import RESULTS

A = RESULTS / "track_a"
OUT = A / "tracka_vet.json"
FIG = RESULTS / "report" / "cases"
MAX_SECTORS = 6


def stacked_centroid(t, F, w, tcs, D, ra, dec, gf):
    from eb_pixels import _diff_image, fit_diff
    num = den = base = None
    n = 0
    for tc in tcs:
        r = _diff_image(t, F, tc, D)
        if r is None:
            continue
        diff, err, b = r
        wt = 1 / np.where(err > 0, err, np.inf) ** 2
        num = diff * wt if num is None else num + diff * wt
        den = wt if den is None else den + wt
        base = b if base is None else base + b
        n += 1
    if not n:
        return None
    with np.errstate(divide="ignore", invalid="ignore"):
        diff = num / den
        err = 1 / np.sqrt(den)
    c = fit_diff(diff, err, base / n, w, float(np.median(tcs)), ra, dec, min_snr=3.0, gf=gf)
    if c:
        c["n_transits"] = n
    return c


def pixels_for(tic, rec, sig):
    import ticdb
    from centroid_vet import cube
    from eb_pixels import centroid_verdict, combine_centroids, gaia_field, neighbour_offsets
    st = ticdb.query([tic], ["ra", "dec", "Tmag"]).iloc[0]
    ra, dec, tmag = float(st.ra), float(st.dec), float(st.Tmag)
    gf = gaia_field(ra, dec, tmag)
    secs = rec["sectors"]
    allsec = sorted(set(secs.get("spoc", [])) | set(secs.get("tspoc", [])) | set(secs.get("qlp", [])))
    # 2-minute pixels first, then the most recent full-frame sectors (finer cadence)
    order = sorted(allsec, key=lambda s: (s not in secs.get("spoc", []), -s))[:MAX_SECTORS]
    P, t0, D = sig["period"], sig["t0"], sig["duration_hr"] / 24
    per, cents = [], []
    for s in order:
        try:
            cb = cube(tic, s, ra, dec, secs)
        except Exception as ex:
            per.append({"sector": s, "error": repr(ex)[:120]})
            continue
        if cb is None:
            per.append({"sector": s, "error": "no pixels"})
            continue
        t, F, w, kind = cb
        k = np.arange(math.ceil((t.min() - t0) / P), math.floor((t.max() - t0) / P) + 1)
        tcs = t0 + k * P
        c = stacked_centroid(t, F, w, tcs, D, ra, dec, gf)
        row = {"sector": s, "kind": kind, "n_predicted": int(len(tcs))}
        if c:
            row.update(c)
            if "dra" in c:
                cents.append(dict(c, sector=s))
        per.append(row)
    cen = combine_centroids(cents)
    cand = neighbour_offsets(gf, sig["depth_ppm"] * 1e-6, 150.0)
    return {"per_sector": per, "centroid": cen, "n_candidates": len(cand), **centroid_verdict(cen, cand)}


def fold_and_depths(tic, rec, sig, others):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from periodic import load_star, prepare_star
    lc = prepare_star(load_star(tic, rec["sectors"]))
    t, f, e, sec = lc["time"], lc["flux"], lc["flux_err"], lc["sector"]
    keep = np.ones(t.size, bool)
    for o in others:
        ph_o = ((t - o["t0"] + 0.5 * o["period"]) % o["period"]) - 0.5 * o["period"]
        keep &= np.abs(ph_o) > o["duration_hr"] / 24
    t, f, e, sec = t[keep], f[keep], e[keep], sec[keep]
    P, t0, D = sig["period"], sig["t0"], sig["duration_hr"] / 24
    ph = ((t - t0 + 0.5 * P) % P) - 0.5 * P
    inn = np.abs(ph) < 0.4 * D
    near = (np.abs(ph) > 0.75 * D) & (np.abs(ph) < 3 * D)
    depths = []
    for s in np.unique(sec):
        m = sec == s
        if (m & inn).sum() < 3 or (m & near).sum() < 6:
            continue
        w_in = 1 / e[m & inn] ** 2
        d = float(np.median(f[m & near]) - np.sum(f[m & inn] * w_in) / w_in.sum())
        de = float(1.2533 * np.std(f[m & near]) / math.sqrt((m & inn).sum()))
        depths.append({"sector": int(s), "depth_ppm": d * 1e6, "err_ppm": de * 1e6, "n_in": int((m & inn).sum())})
    # half-phase and odd/even on the same footing
    ph2 = ((t - t0) % P) - 0.5 * P
    sec_in = np.abs(ph2) < 0.4 * D
    ep = np.round((t - t0) / P).astype(int)

    def box(sel):
        if sel.sum() < 3:
            return None, None
        ww = 1 / e[sel] ** 2
        base = np.median(f[near]) if near.sum() else 1.0
        return float((base - np.sum(f[sel] * ww) / ww.sum()) * 1e6), float(1 / math.sqrt(ww.sum()) * 1e6)
    sec_d = box(sec_in)
    odd, even = box(inn & (ep % 2 == 1)), box(inn & (ep % 2 == 0))
    # figure: folded, binned
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.4))
    x = ph * 24
    win = np.abs(x) < max(4 * sig["duration_hr"], 6)
    ax[0].plot(x[win], (f[win] - 1) * 1e6, ".", ms=1.5, color="0.75", rasterized=True)
    bw = max(sig["duration_hr"] / 6, 0.25)
    edges = np.arange(x[win].min(), x[win].max() + bw, bw)
    idx = np.digitize(x[win], edges)
    xb = [x[win][idx == i].mean() for i in np.unique(idx) if (idx == i).sum() > 3]
    yb = [(np.median(f[win][idx == i]) - 1) * 1e6 for i in np.unique(idx) if (idx == i).sum() > 3]
    ax[0].plot(xb, yb, "o", ms=3.5, color="#b5402a")
    lim = 3 * sig["depth_ppm"]
    ax[0].set_ylim(-lim, lim * 0.8)
    ax[0].set_xlabel("hours from mid-transit")
    ax[0].set_ylabel("ppm")
    ax[0].set_title(f"TIC {tic}  P = {P:.5f} d  ({len(np.unique(sec))} sectors)", fontsize=9)
    if depths:
        ax[1].errorbar([d["sector"] for d in depths], [d["depth_ppm"] for d in depths],
                       [d["err_ppm"] for d in depths], fmt="o", ms=3, color="#2a5db5")
        ax[1].axhline(sig["depth_ppm"], color="0.5", lw=0.8)
        ax[1].axhline(0, color="0.8", lw=0.8)
        ax[1].set_xlabel("sector")
        ax[1].set_ylabel("depth (ppm)")
    fig.tight_layout()
    p = FIG / f"tracka_{tic}_{P:.3f}.png"
    fig.savefig(p, dpi=110)
    plt.close(fig)
    return {"sector_depths": depths, "secondary_ppm": sec_d, "odd_ppm": odd, "even_ppm": even, "figure": str(p.name)}


def gaia_for(tic):
    import gaia_check
    import ticdb
    st = ticdb.query([tic], ["ra", "dec", "Tmag"]).iloc[0]
    g = gaia_check.match(float(st.ra), float(st.dec), float(st.Tmag))
    if g is None:
        return {"gaia": "no match"}
    keep = ("source_id", "ruwe", "non_single_star", "radial_velocity_error", "rv_nb_transits", "phot_g_mean_mag",
            "parallax", "ipd_frac_multi_peak", "ipd_gof_harmonic_amplitude")
    return {k: (g[k].item() if hasattr(g[k], "item") else g[k]) for k in keep if k in g.index}


def exofop_for(tic):
    import exofop_check
    from env import CATALOGS
    p = CATALOGS / "exofop" / "targets" / f"{tic}.json"
    if p.exists():
        p.unlink()          # fetch today's page, not a cached one
    try:
        d = exofop_check.fetch(tic)
        return exofop_check.summarise(tic, d)
    except Exception as ex:
        return {"error": repr(ex)[:150]}


def vet_one(tic, P):
    rec = json.load(open(A / "stars" / f"{tic}.json"))
    sigs = rec["signals"]
    sig = min(sigs, key=lambda s: abs(s["period"] / P - 1))
    others = [s for s in sigs if s is not sig]
    out = {"tic": tic, "period": sig["period"], "t0": sig["t0"], "duration_hr": sig["duration_hr"],
           "depth_ppm": sig["depth_ppm"], "sde": sig["sde"]}
    for name, fn in (("lc", lambda: fold_and_depths(tic, rec, sig, others)), ("pixels", lambda: pixels_for(tic, rec, sig)),
                     ("gaia", lambda: gaia_for(tic)), ("exofop", lambda: exofop_for(tic))):
        try:
            out[name] = fn()
        except Exception as ex:
            out[name] = {"error": repr(ex)[:200]}
    return out


def main(args):
    if args:
        todo = [(int(a.split(":")[0]), float(a.split(":")[1])) for a in args]
    else:
        v2 = json.load(open(A / "tracka_v2.json"))
        todo = [(int(r["tic"]), float(r["period"])) for r in v2["unknown_above"]]
    res = json.load(open(OUT)) if OUT.exists() else {}
    todo = [x for x in todo if f"{x[0]}|{x[1]:.4f}" not in res]
    print(len(todo), "to vet", flush=True)
    with ProcessPoolExecutor(6) as ex:
        futs = {ex.submit(vet_one, tic, P): (tic, P) for tic, P in todo}
        for fu in as_completed(futs):
            tic, P = futs[fu]
            try:
                r = fu.result()
            except Exception as e:
                r = {"tic": tic, "period": P, "error": repr(e)[:200]}
            res[f"{tic}|{P:.4f}"] = r
            json.dump(res, open(OUT, "w"), indent=1, default=float)
            pv = r.get("pixels", {})
            print(tic, f"P={P:.4f}", pv.get("verdict"), "RUWE", r.get("gaia", {}).get("ruwe"), flush=True)


if __name__ == "__main__":
    main(sys.argv[1:])
