"""Draft ExoFOP bulk planet-parameter upload (params_planet_YYYYMMDD_001.txt format, fetched from ExoFOP's template
on 15 October 2026) for the results that may be uploaded without a publication: additional parameter sets
(flag = newparams) on existing TOIs and community TOIs. ExoFOP requires a paper only for flag = newctoi, so TIC
402710280.01 is not here.

Field definitions used (ExoFOP's units):
  depth     transit depth in ppm as seen in the light curve: the fitted limb-darkened model's depth at mid-transit
            (not (Rp/R*)^2, which goes in r_planet as Rp/R*); for the TOI ephemeris rows, the mean in-transit depth
            of a trapezoid fitted to each transit (toi_audit_timing.py), weighted over transits;
  epoch     BJD_TDB of an observed transit; period and epoch errors from the fits (scaled by sqrt(chi2_red) when
            above 1);
  duration  hours, first to fourth contact, from the transit fit (left blank where only timed at the listed
            duration);
  dens      fitted stellar density in g/cm^3 (solar density 1.408 g/cm^3);
  radius    R_Earth from Rp/R* and the TIC stellar radius (its error included).
Every row carries the same placeholder data tag (replace with your own before uploading), prop_period 0, no paper,
and a note of at most 120 characters. A row is written only when its checks pass; the companion CSV lists every
candidate row with the reason it was or was not included, and toi_review/ holds an O-C plot per TOI row.

  python exofop_upload.py
Output: results/candidates/exofop_upload/params_planet_DRAFT.txt, upload_review.csv, toi_review/*.png
"""
import json
import math

import env  # noqa: F401
import numpy as np
import pandas as pd

import transit_fit as tf
from env import RESULTS

C = RESULTS / "candidates"
D = C / "exofop_upload"
RHO_SUN = 1.408
TAG = "YYYYMMDD_username_adastra_nnnnn"
COLS = ("target|flag|disp|discovery|candname|period|period_unc|epoch|epoch_unc|depth|depth_unc|duration|duration_unc|"
        "inc|inc_unc|imp|imp_unc|r_planet|r_planet_err|ar_star|ar_star_err|radius|radius_unc|mass|mass_unc|temp|temp_unc|"
        "insol|insol_unc|dens|dens_unc|sma|sma_unc|ecc|ecc_unc|arg_peri|arg_peri_err|time_peri|time_peri_unc|vsa|vsa_unc|"
        "tag|group|prop_period|paper|notes").split("|")


def err(v):
    return max(v[1], v[2])


def ar_from_duration(k, b, t14_h, P):
    """a/R* that gives this total duration at period P for a circular orbit (shape-only fits have no period)."""
    x = math.sqrt(max((1 + k) ** 2 - b ** 2, 1e-9))
    s = math.sin(math.pi * t14_h / 24 / P)
    return math.sqrt((x / s) ** 2 + b ** 2)


def model_depth(f, P=None):
    """Depth at mid-transit of the fitted limb-darkened model, and its error from Rp/R*'s."""
    P = P or (f["period_d"][0] if "period_d" in f else 300.0)
    k, b = f["rp_rs"][0], f["b"][0]
    aR = f["a_rs"][0] if "a_rs" in f else ar_from_duration(k, b, f["t14_h"][0], P)
    inc = math.degrees(math.acos(min(b / aR, 1)))
    d = 1 - float(tf._model(np.array([0.0]), 0.0, P, k, aR, inc, f["u"])[0])
    return d * 1e6, d * 1e6 * 2 * err(f["rp_rs"]) / k


def row(**kw):
    r = {c: "" for c in COLS}
    r.update({k: ("" if v is None else v) for k, v in kw.items()})
    r.update(tag=TAG, prop_period=0)
    assert len(str(r["notes"])) <= 120, (r["target"], len(r["notes"]))
    return r


def fmt(x, n):
    return "" if x is None or x != x else f"{x:.{n}f}"


def periodic_row(target, f, ex, rstar, rstar_e, note, disp="PC"):
    dep, dep_e = model_depth(f)
    k = f["rp_rs"][0]
    rp = k * rstar * 109.08
    rp_e = rp * math.hypot(err(f["rp_rs"]) / k, rstar_e / rstar)
    return row(target=target, flag="newparams", disp=disp, period=fmt(f["period_d"][0], 6), period_unc=fmt(err(f["period_d"]), 6),
               epoch=fmt(f["tref_btjd"][0] + 2457000, 5), epoch_unc=fmt(err(f["tref_btjd"]), 5), depth=fmt(dep, 0),
               depth_unc=fmt(dep_e, 0), duration=fmt(f["t14_h"][0], 3), duration_unc=fmt(err(f["t14_h"]), 3),
               imp=fmt(f["b"][0], 3), imp_unc=fmt(err(f["b"]), 3), r_planet=fmt(k, 5), r_planet_err=fmt(err(f["rp_rs"]), 5),
               ar_star=fmt(f["a_rs"][0], 2), ar_star_err=fmt(err(f["a_rs"]), 2), radius=fmt(rp, 2), radius_unc=fmt(rp_e, 2),
               temp=fmt(ex.get("teq"), 0), insol=fmt(ex.get("s_earth"), 2),
               dens=fmt(f["rho_solar"][0] * RHO_SUN, 3), dens_unc=fmt(err(f["rho_solar"]) * RHO_SUN, 3),
               sma=fmt(ex.get("a_au"), 4), notes=note)


def oc_plot(r, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ev = r["events"]
    n = np.array([e["n"] for e in ev])
    oc = np.array([e["o_c_min"] for e in ev])
    er = np.array([e["err"] * 1440 for e in ev])
    fig, ax = plt.subplots(1, 2, figsize=(9, 3))
    ax[0].errorbar(n, oc, er, fmt="o", ms=3, color="#2a5db5")
    ax[0].axhline(0, color="0.6", lw=0.8)
    # where the listed ephemeris would put each transit, relative to the corrected one
    if not r.get("no_listed"):
        Pl, Tl = r["P_listed"], r["T0_listed"]
        tt = np.array([e["tc"] for e in ev])
        kl = np.round((tt - Tl) / Pl)
        ax[0].plot(n, (Tl + kl * Pl - (tt - oc / 1440)) * 1440, "x", color="#b5402a", ms=4, label="listed ephemeris")
    ax[0].set_xlabel("cycle")
    ax[0].set_ylabel("O-C (min)")
    if not r.get("no_listed"):
        ax[0].legend(fontsize=7)
    ax[0].set_title(f"{r.get('label') or 'TOI-' + str(r['toi'])}  P = {r['P']:.6f} ± {r['P_err']:.6f} d  (χ²_red {r['chi2_red']:.1f})", fontsize=8)
    odd = n % 2 == 1
    for sel, col, lab in ((~odd, "#2a5db5", "even"), (odd, "#b5402a", "odd")):
        ax[1].errorbar(n[sel], np.array([e["depth"] for e in ev])[sel] * 1e6, np.array([e["depth_err"] for e in ev])[sel] * 1e6,
                       fmt="o", ms=3, color=col, label=lab)
    ax[1].axhline(r["depth_ppm_listed"], color="0.6", lw=0.8, label="listed depth")
    ax[1].set_xlabel("cycle (odd/even)")
    ax[1].set_ylabel("depth (ppm)")
    ax[1].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(path, dpi=100)
    plt.close(fig)


def main():
    D.mkdir(parents=True, exist_ok=True)
    (D / "toi_review").mkdir(exist_ok=True)
    fits = json.load(open(RESULTS / "report" / "transit_fits.json"))
    ex = json.load(open(RESULTS / "report" / "extras.json"))
    cf = json.load(open(C / "corrections_fit.json"))
    rows, review = [], []

    def add(r, why):
        rows.append(r)
        review.append({"target": r["target"], "included": True, "why": why, "notes": r["notes"]})

    def skip(target, why):
        review.append({"target": target, "included": False, "why": why})

    # TOI-2493.01: period for a single-transit TOI
    f = fits["TOI2493"]
    ok = f.get("chain_summary", {}).get("converged")
    r = periodic_row("TOI2493.01", f, ex["TOI2493"], 0.652, 0.059,
                     "Period from 3 transits S33,S98,S107 (2-min); shorter aliases excluded incl. S108. Ad Astra survey")
    add(r, "3 transits, aliases excluded; MCMC " + ("converged" if ok else "convergence not confirmed"))
    # HIP 56966 = CTOI 266376971.01: period
    f = fits["HIP56966"]
    r = periodic_row("TIC266376971.01", f, ex["HIP56966"], 1.069, 0.05,
                     "Period from 4 transits S10,S64,S99,S100; centroid offset marginal (2.7 sigma). Ad Astra survey")
    add(r, "4 transits; centroid caveat in note")
    # HD 11506 = CTOI 266680498.01: second transit, shape; period unknown
    f = fits["HD11506"]
    k = f["rp_rs"][0]
    dep = (k ** 2) * 1e6
    dep, dep_e = model_depth(f, 300.0)
    rp = k * 1.35 * 109.08
    add(row(target="TIC266680498.01", flag="newparams", disp="PC", epoch=fmt(f["tc_1"][0] + 2457000, 5),
            epoch_unc=fmt(err(f["tc_1"]), 5), depth=fmt(dep, 0), depth_unc=fmt(dep_e, 0), duration=fmt(f["t14_h"][0], 3),
            duration_unc=fmt(err(f["t14_h"]), 3), imp=fmt(f["b"][0], 3), imp_unc=fmt(err(f["b"]), 3),
            r_planet=fmt(k, 5), r_planet_err=fmt(err(f["rp_rs"]), 5), radius=fmt(rp, 2),
            radius_unc=fmt(rp * math.hypot(err(f["rp_rs"]) / k, 0.06 / 1.35), 2),
            notes="2nd transit S107 (1st S30); P=2149.389/n d, 32 values >=25.9 d allowed. Ad Astra survey"),
        "second transit; period left blank (aliases)")
    # CTOI 438225366.01: recurrence at half the posted period, both interpretations in the note
    a = cf["CTOI 438225366.01"]
    r = row(target="TIC438225366.01", flag="newparams", disp="PC", period=fmt(a["period_d"][0], 6),
            period_unc=fmt(err(a["period_d"]), 6), epoch=fmt(a["tref_btjd"][0] + 2457000, 5), epoch_unc=fmt(err(a["tref_btjd"]), 5),
            depth=fmt(a["depth_ppm"][0], 0), depth_unc=fmt(err(a["depth_ppm"]), 0), duration=fmt(a["t14_h"][0], 3),
            duration_unc=fmt(err(a["t14_h"]), 3), imp=fmt(a["b"][0], 3), imp_unc=fmt(err(a["b"]), 3),
            r_planet=fmt(a["rp_rs"][0], 5), r_planet_err=fmt(err(a["rp_rs"]), 5),
            notes="8 equal dips every 20.8615 d. Orbit 20.86 d if transits, or ~41.72 d if twin-EB eclipses. Ad Astra")
    r["depth"] = fmt(model_depth(a)[0], 0)
    add(r, "recurrence; EB alternative in note")
    # TOI-2633.01: ephemeris
    a = cf["TOI-2633.01"]
    r = row(target="TOI2633.01", flag="newparams", disp="PC", period=fmt(a["period_d"][0], 6), period_unc=fmt(err(a["period_d"]), 6),
            epoch=fmt(a["tref_btjd"][0] + 2457000, 5), epoch_unc=fmt(err(a["tref_btjd"]), 5), depth=fmt(model_depth(a)[0], 0),
            duration=fmt(a["t14_h"][0], 3), duration_unc=fmt(err(a["t14_h"]), 3), r_planet=fmt(a["rp_rs"][0], 5),
            r_planet_err=fmt(err(a["rp_rs"]), 5),
            notes="Ephemeris from 6 transits S7-S88; listed T0 is off-phase. Density 3x TIC: blend/ecc possible. Ad Astra")
    add(r, "6 transits; density caveat in note")
    # TOI ephemeris corrections (model timing)
    tm = json.load(open(C / "toi_audit_timing.json"))
    for key, t in sorted(tm.items(), key=lambda kv: -abs(kv[1].get("sigma_listed") or 0)):
        tgt = f"TOI{t.get('toi')}"
        if t.get("status") != "corrected":
            continue
        why = []
        if not t.get("solid"):
            why.append("too few or too sparse transits, or chi2_red >= 5")
        if t.get("already_posted"):
            why.append("matching period already posted: " + "; ".join(t["already_posted"])[:80])
        oe = t.get("odd_even_sigma")
        if oe and oe > 3:
            why.append(f"odd/even depths differ {oe:.1f} sigma (EB at 2P?)")
        if t.get("chi2_red", 0) >= 5:
            why.append(f"timing chi2_red {t['chi2_red']:.1f}")
        if why:
            skip(tgt, "; ".join(why))
            if t.get("events") and len(t["events"]) >= 3 and "P" in t:
                (D / "toi_review_held").mkdir(exist_ok=True)
                oc_plot(t, D / "toi_review_held" / f"TOI{t['toi']}.png")
            continue
        secs = sorted({e["sector"] for e in t["events"] if e.get("sector")})
        srange = f"S{secs[0]}-S{secs[-1]}" if secs else ""
        note = (f"Alt. ephemeris from {t['n_timed']} transit times {srange}; listed off {abs(t['offset_h']):.0f} h on "
                f"2026-12-01. Ad Astra")
        r = row(target=tgt, flag="newparams", disp="PC", period=fmt(t["P"], 7), period_unc=fmt(t["P_err"], 7),
                epoch=fmt(t["T0"] + 2457000, 5), epoch_unc=fmt(t["T0_err"], 5), depth=fmt(t["depth_ppm"], 0), notes=note)
        add(r, f"{t['n_timed']} transits, chi2_red {t['chi2_red']:.1f}, offset {t['offset_h']:+.1f} h = {t['sigma_listed']:.0f}x listed error")
        oc_plot(t, D / "toi_review" / f"TOI{t['toi']}.png")
    # periods for period-less TOIs and community TOIs (period_signal_review.py, tier "ready")
    psr = json.load(open(C / "period_signal_review.json"))
    (D / "period_review").mkdir(exist_ok=True)
    for key, x in sorted(psr.items()):
        name = str(x.get("id", ""))
        tgt = ("TOI" + name.split("TOI-")[1]) if name.startswith("TOI-") else ("TIC" + name.split("CTOI ")[1]) if name.startswith("CTOI ") else name
        if x.get("tier") != "ready":
            skip(tgt, f"period review: {x.get('tier')}: {x.get('why') or x.get('error', '')}"[:200])
            continue
        t = x["trials"][0]
        secs = sorted({e["sector"] for e in t["events"] if e.get("sector")})
        note = f"Period from {t['n_timed']} transit times S{secs[0]}-S{secs[-1]}; no other period fits. Ad Astra survey"
        add(row(target=tgt, flag="newparams", disp="", period=fmt(t["P"], 6), period_unc=fmt(t["P_err"], 6),
                epoch=fmt(t["T0"] + 2457000, 5), epoch_unc=fmt(t["T0_err"], 5), depth=fmt(t["depth_ppm"], 0),
                depth_unc=fmt(t["depth_err_ppm"], 0), notes=note),
            f"period review ready: {t['n_timed']} transits, chi2_red {t['chi2_red']:.1f}, pixels {x.get('pixels')}")
        ev = [dict(e, o_c_min=0.0) for e in t["events"]]
        n = np.array([e["n"] for e in ev], float)
        tt = np.array([e["tc"] for e in ev])
        k = np.round((tt - t["T0"]) / t["P"])
        for e, kk in zip(ev, k):
            e["o_c_min"] = float((e["tc"] - t["T0"] - kk * t["P"]) * 1440)
        oc_plot({"events": ev, "P_listed": t["P"], "T0_listed": t["T0"], "toi": name, "P": t["P"], "P_err": t["P_err"],
                 "chi2_red": t["chi2_red"], "depth_ppm_listed": x["depth_ref_ppm"], "label": name, "no_listed": True},
                D / "period_review" / f"{tgt}.png")
    with open(D / "params_planet_DRAFT.txt", "w") as fh:
        fh.write("\\ DRAFT, not for upload as is: replace the tag, review each row (upload_review.csv, toi_review/).\n")
        fh.write("|".join(COLS) + "\n")
        for r in rows:
            fh.write("|".join(str(r[c]) for c in COLS) + "\n")
    pd.DataFrame(review).to_csv(D / "upload_review.csv", index=False)
    inc = sum(1 for x in review if x["included"])
    print(inc, "rows written;", len(review) - inc, "held back")
    print(pd.DataFrame(review).groupby("included").size())


if __name__ == "__main__":
    main()
