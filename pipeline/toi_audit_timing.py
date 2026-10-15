"""Transit-model timing for the TOI ephemeris audit (v2 of toi_audit_review.py's timing step).

The first pass timed each transit by maximising a box statistic on 41 grid points across one duration, with an
error D/(2 SNR) floored at half a cadence: good enough to flag a wrong ephemeris, not to publish one (independent
audit, 15 October 2026). Here every transit the chosen ephemeris predicts inside the data is fitted with a
trapezoid (total duration D, ingress and egress D/8 each) whose mid-time is scanned on a 1-minute grid over +-0.6 D,
with the depth amplitude and a quadratic baseline refitted at every step. The 1-sigma error is the delta-chi2 = 1
half-width of a parabola through the minimum, multiplied by the red-noise factor beta (residual rms in bins of
D/8 to D/2 against white noise). A transit is kept when its depth amplitude is at least 3 sigma, its minimum is
not at the grid edge and both sides of the transit have data. The ephemeris is a weighted straight line with its
full 2x2 covariance, errors scaled by sqrt(chi2_red) when above 1, and every transit's O-C is kept for review.
Odd and even depths are compared from the same fits.

  N_PROC=6 python toi_audit_timing.py            (every TOI the first pass flagged and fitted)
Output: results/candidates/toi_audit_timing.json and .csv
"""
import json
import math
import os
from concurrent.futures import ProcessPoolExecutor, as_completed

import env  # noqa: F401
import numpy as np
import pandas as pd

from env import RESULTS
from toi_audit_review import OUT as V1, REF, _check_ref, classify

C = RESULTS / "candidates"
OUT = C / "toi_audit_timing.json"


def trapezoid(t, tc, D):
    """Unit-depth trapezoid: 1 inside, 0 outside, linear ingress/egress of D/8."""
    x = np.abs(t - tc)
    tau = D / 8
    return np.clip((D / 2 - x) / tau, 0, 1)


def beta_factor(t, r, sig, D):
    cad = float(np.median(np.diff(t)))
    bs = []
    for w in (D / 8, D / 4, D / 2):
        n = max(int(round(w / cad)), 1)
        k = len(r) // n
        if n < 2 or k < 4:
            continue
        rb = r[: k * n].reshape(k, n).mean(axis=1)
        bs.append(float(np.std(rb) / (sig / math.sqrt(n)) * math.sqrt(k / (k - 1))))
    return max(1.0, float(np.median(bs))) if bs else 1.0


def time_one(lc, tc0, D):
    t, f = lc["time"], lc["flux"]
    w = np.isfinite(t) & np.isfinite(f) & (np.abs(t - tc0) < 1.5 * D + 0.3)
    t, f = t[w], f[w]
    if t.size < 30:
        return None
    if (t < tc0 - 0.6 * D).sum() < 5 or (t > tc0 + 0.6 * D).sum() < 5 or (np.abs(t - tc0) < 0.4 * D).sum() < 3:
        return None
    f = f / np.median(f)
    oot = np.abs(t - tc0) > 0.75 * D
    sig = 1.4826 * np.median(np.abs(np.diff(f[oot]))) / math.sqrt(2)
    if not sig > 0:
        return None
    grid = tc0 + np.arange(-0.6 * D, 0.6 * D + 1e-9, 1 / 1440)
    chi, amps = [], []
    for tc in grid:
        dt = t - tc
        A = np.vstack([np.ones_like(t), dt, dt ** 2, -trapezoid(t, tc, D)]).T
        c, *_ = np.linalg.lstsq(A, f, rcond=None)
        chi.append(float(np.sum(((f - A @ c) / sig) ** 2)))
        amps.append(c[3])
    chi = np.array(chi)
    i = int(np.argmin(chi))
    if i < 3 or i > len(grid) - 4:
        return {"edge": True}
    lo, hi = max(i - 8, 0), min(i + 9, len(grid))
    a2, a1, _ = np.polyfit(grid[lo:hi] - grid[i], chi[lo:hi], 2)
    if a2 <= 0:
        return None
    tb = float(grid[i] - a1 / (2 * a2))
    dt = t - tb
    A = np.vstack([np.ones_like(t), dt, dt ** 2, -trapezoid(t, tb, D)]).T
    c, *_ = np.linalg.lstsq(A, f, rcond=None)
    r = f - A @ c
    beta = beta_factor(t, r, sig, D)
    cov = np.linalg.inv(A.T @ A) * (sig * beta) ** 2
    amp, amp_e = float(c[3]), float(math.sqrt(cov[3, 3]))
    return {"tc": tb, "err": float(beta / math.sqrt(a2)), "beta": beta, "depth": amp, "depth_err": amp_e,
            "snr": amp / amp_e if amp_e > 0 else 0.0, "edge": False,
            "dchi2_flat": float(np.sum(((f - A[:, :3] @ np.linalg.lstsq(A[:, :3], f, rcond=None)[0]) / sig) ** 2) - chi[i])}


def fit_line(n, tt, st):
    n0 = float(np.round(np.average(n, weights=1 / st ** 2)))
    A = np.vstack([np.ones_like(n), n - n0]).T
    W = np.diag(1 / st ** 2)
    cov = np.linalg.inv(A.T @ W @ A)
    p = cov @ A.T @ W @ tt
    res = tt - A @ p
    chi2 = float(np.sum((res / st) ** 2))
    dof = max(len(tt) - 2, 1)
    s2 = max(chi2 / dof, 1.0)
    return n0, p, cov * s2, chi2, dof, res


def review(key, r1):
    import toi_audit as ta
    from dossier import all_sectors, best_lcs
    from period_review import lc_for
    tic = int(r1["tic"])
    D = r1["dur_h"] / 24
    lcs = [lc for lc in best_lcs(tic, all_sectors(tic)) if lc["time"].size > 100]
    P1, T1 = r1["P"], r1["T0"]
    preds = ta.predicted(lcs, T1, P1, D)
    if len(preds) > 400:     # short periods over many sectors: an even subsample keeps the run bounded
        idx = np.linspace(0, len(preds) - 1, 400).astype(int)
        preds = [preds[i] for i in idx]
    ev = []
    for n, tc in preds:
        lc = lc_for(lcs, tc)
        if lc is None:
            continue
        m = time_one(lc, tc, D)
        if m and not m.get("edge") and m["snr"] >= 3:
            ev.append({"n": int(n), "sector": lc["meta"].get("sector"), **{k: m[k] for k in ("tc", "err", "beta", "depth", "depth_err", "snr")}})
    out = {"key": key, "tic": tic, "toi": r1["toi"], "P_listed": r1["P_listed"], "T0_listed": r1["T0_listed"],
           "dur_h": r1["dur_h"], "depth_ppm_listed": r1["depth_ppm"], "n_predicted": len(preds), "n_timed": len(ev),
           "events": ev, "P_v1": P1, "T0_v1": T1}
    if len(ev) < 3:
        out["status"] = "too few timed transits"
        return out
    n = np.array([e["n"] for e in ev], float)
    tt = np.array([e["tc"] for e in ev])
    st = np.array([e["err"] for e in ev])
    n0, p, cov, chi2, dof, res = fit_line(n, tt, st)
    # the epoch at the weighted-mean cycle: period-epoch covariance there is near zero, but it is kept anyway
    T0, P = float(p[0]), float(p[1])      # T0: the transit at cycle n0 (v1 numbering), where covariance is smallest
    out.update(P=P, P_err=float(math.sqrt(cov[1, 1])), T0=T0, T0_err=float(math.sqrt(cov[0, 0])),
               cov_T0_P=float(cov[0, 1]), ref_cycle=n0, chi2=chi2, dof=dof, chi2_red=chi2 / dof,
               o_c_min=[float(x * 1440) for x in res], o_c_rms_min=float(np.sqrt(np.mean(res ** 2)) * 1440))
    for e, x in zip(ev, res):
        e["o_c_min"] = float(x * 1440)
    dep = np.array([e["depth"] for e in ev])
    de = np.array([e["depth_err"] for e in ev])
    odd = (n.astype(int) % 2) == 1
    if odd.sum() >= 2 and (~odd).sum() >= 2:
        wo, we = 1 / de[odd] ** 2, 1 / de[~odd] ** 2
        mo, me = np.sum(wo * dep[odd]) / wo.sum(), np.sum(we * dep[~odd]) / we.sum()
        out["odd_even_sigma"] = float(abs(mo - me) / math.hypot(1 / math.sqrt(wo.sum()), 1 / math.sqrt(we.sum())))
        out["odd_even_ppm"] = [float(mo * 1e6), float(me * 1e6)]
    w = 1 / de ** 2
    out["depth_ppm"] = float(np.sum(w * dep) / w.sum() * 1e6)
    # 1 Dec 2026 comparison with the full covariance for the corrected ephemeris
    k = round((REF - T0) / P)
    out.update(classify(T0, P, out["T0_err"], out["P_err"], r1["T0_listed"], r1["P_listed"], r1.get("_Pe_l"),
                        r1.get("_T0e_l"), D))
    out["corrected_err_h"] = float(math.sqrt(max(cov[0, 0] + k * k * cov[1, 1] + 2 * k * cov[0, 1], 0)) * 24)
    out["solid"] = bool(len(ev) >= 5 and len(ev) >= 0.7 * r1["chosen"]["n_data"] and out["chi2_red"] < 5)
    out["P_shift_sigma_vs_v1"] = float(abs(P - P1) / max(out["P_err"], 1e-9))
    return out


def main():
    import toi_audit as ta
    _check_ref()
    v1 = json.load(open(V1))
    t = ta.tois()
    cat = {f"{float(x):.2f}": (float(a) if a == a else None, float(b) if b == b else None)
           for x, a, b in zip(t["TOI"], t["Period (days) err"], t["Epoch (BJD) err"])}
    todo = {}
    for k, r in v1.items():
        if r.get("status") in ("corrected", "refined", "within a transit duration") and "P" in r:
            r = dict(r)
            r["_Pe_l"], r["_T0e_l"] = cat.get(r["toi"], (None, None))
            todo[k] = r
    res = json.load(open(OUT)) if OUT.exists() else {}
    todo = {k: r for k, r in todo.items() if k not in res}
    print(len(todo), "to time", flush=True)
    with ProcessPoolExecutor(int(os.environ.get("N_PROC", "6"))) as ex:
        futs = {ex.submit(review, k, r): k for k, r in todo.items()}
        for i, fu in enumerate(as_completed(futs)):
            k = futs[fu]
            try:
                res[k] = fu.result()
            except Exception as e:
                res[k] = {"key": k, "status": "error", "error": repr(e)[:200]}
            if i % 20 == 0:
                json.dump(res, open(OUT, "w"), indent=1, default=float)
                print(i, "/", len(todo), flush=True)
    for k, r in res.items():
        if k in v1:
            r["already_posted"] = v1[k].get("already_posted", [])
    json.dump(res, open(OUT, "w"), indent=1, default=float)
    rows = [{a: b for a, b in r.items() if not isinstance(b, (list, dict))} for r in res.values()]
    pd.DataFrame(rows).to_csv(C / "toi_audit_timing.csv", index=False)
    print(pd.Series([r.get("status") for r in res.values()]).value_counts().to_string())


if __name__ == "__main__":
    main()
