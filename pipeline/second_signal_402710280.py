"""The 3.78-day transit signal of TIC 402710280 (the second, weaker periodicity found once the 36.34-day transits of
TIC 402710280.01 are masked).

Steps:
  1. light curves of all eight sectors (TESS-SPOC 12, 39, 66; QLP 93; our FFI photometry 100-103), the 36.34-d
     transits masked (+-6 h), each divided by a 0.75-d running median that ignores the 3.78-d transit phases;
  2. independent box-least-squares searches of sectors 12-93 and 100-103, and a refined ephemeris from all data;
  3. depth in every sector, and in the epoch groups n mod k (k = 2..5: odd/even and the P*k aliases), with
     errors from boxes placed at 400 random phases (so they include correlated noise);
  4. secondary eclipse at phase 0.5 and the deepest box over phases 0.3-0.7;
  5. transit fit (batman/emcee, fit_periodic_windows) with the stellar density free, compared with the star's;
  6. the dilution required if the 3.1-3.6" background neighbour (G = 16.15) were the source.

  python second_signal_402710280.py [nstep]
  python second_signal_402710280.py --refit nstep     (step 5 only, on the stored ephemeris)
Output: results/candidates/second_signal_402710280.json, results/report/cases/TIC402710280_3.78d.png
"""
import json
import sys

import env  # noqa: F401
import numpy as np
from astropy.timeseries import BoxLeastSquares

import transit_fit as tf
from env import RESULTS
from lcaccess import lcs_for

TIC = 402710280
P1, T1 = 36.33989590709609, 2460838.20190249 - 2457000      # TIC 402710280.01 (BTJD)
TC1 = [2384.6094, 3838.1847, 4092.5851, 4128.9449, 3111.4]    # measured mid-times; S66 predicted
OUT = RESULTS / "candidates" / "second_signal_402710280.json"
FIG = RESULTS / "report" / "cases" / "TIC402710280_3.78d.png"
TEFF = 5182.0
RHO_STAR = (1.61, 0.38)       # TIC (M, R)
NEIGH_FRAC = 10 ** (-6.37 / 2.5)


def phase_dist(t, P, t0):
    return np.abs(((t - t0 + P / 2) % P) - P / 2)


def running_median(t, f, mask, w=0.75):
    o = np.argsort(t)
    t, f, mask = t[o], f[o], mask[o]
    lo = np.searchsorted(t, t - w / 2)
    hi = np.searchsorted(t, t + w / 2)
    tr = np.full(t.size, np.nan)
    for i in range(t.size):
        seg = f[lo[i]:hi[i]][~mask[lo[i]:hi[i]]]
        if seg.size > 10:
            tr[i] = np.median(seg)
    back = np.empty_like(o)
    back[o] = np.arange(o.size)
    return tr[back]


def load(P2, T2):
    lcs = lcs_for(TIC)
    rows = []
    for lc in lcs:
        t, f, e = (np.asarray(lc[k], float) for k in ("time", "flux", "flux_err"))
        ok = np.isfinite(t) & np.isfinite(f) & np.isfinite(e)
        t, f, e = t[ok], f[ok], e[ok]
        keep = np.min(np.abs(t[:, None] - np.array(TC1)[None, :]), axis=1) > 0.25
        keep &= phase_dist(t, P1, T1) > 0.25
        t, f, e = t[keep], f[keep], e[keep]
        tr = running_median(t, f, phase_dist(t, P2, T2) < 0.06)
        ok = np.isfinite(tr)
        s = lc["meta"]["sector"]
        rows.append((s, t[ok], f[ok] / tr[ok], e[ok] / np.median(f), float(lc["meta"]["cadence_s"])))
    return lcs, rows


def box_depth(t, f, P, t0, dur):
    inn = phase_dist(t, P, t0) < dur / 2
    out = (phase_dist(t, P, t0) > dur) & (phase_dist(t, P, t0) < 0.5)
    if inn.sum() < 3:
        return np.nan, 0
    return float((np.median(f[out]) - np.mean(f[inn])) * 1e6), int(inn.sum())


def depth_with_error(t, f, P, t0, dur, rng, n=400):
    d, nin = box_depth(t, f, P, t0, dur)
    shifts = rng.uniform(0.1, 0.9, n) * P
    null = np.array([box_depth(t, f, P, t0 + s, dur)[0] for s in shifts])
    null = null[np.isfinite(null)]
    return d, float(1.4826 * np.median(np.abs(null - np.median(null)))), nin


def main(nstep=12000):
    rng = np.random.default_rng(3)
    P2, T2 = 3.7837134787718334, 1628.3234410499908          # tracka_vet ephemeris (sectors 12-93)
    lcs, rows = load(P2, T2)
    T = np.concatenate([r[1] for r in rows])
    F = np.concatenate([r[2] for r in rows])
    E = np.concatenate([r[3] for r in rows])
    S = np.concatenate([[r[0]] * len(r[1]) for r in rows])
    out = {"tic": TIC, "sectors": [int(r[0]) for r in rows]}

    # 2. independent searches and the joint ephemeris
    per = np.exp(np.linspace(np.log(0.7), np.log(20), 200000))
    searches = {}
    for name, sel in (("sectors 12-93", S <= 93), ("sectors 100-103", S >= 100)):
        r = BoxLeastSquares(T[sel], F[sel], E[sel]).power(per, [0.03, 0.05, 0.08, 0.12], objective="snr")
        sde = (r.power - np.median(r.power)) / np.std(r.power)
        top = []
        for i in np.argsort(r.power)[::-1]:
            if all(abs(per[i] / q["P"] - 1) > 0.01 for q in top):
                top.append({"P": float(per[i]), "snr": float(r.power[i]), "sde": float(sde[i]),
                            "depth_ppm": float(r.depth[i] * 1e6), "dur_h": float(r.duration[i] * 24),
                            "t0": float(r.transit_time[i])})
            if len(top) == 5:
                break
        searches[name] = top
        print(name, [(round(q["P"], 5), round(q["snr"], 1)) for q in top], flush=True)
    out["searches"] = searches
    fine = np.linspace(P2 - 0.002, P2 + 0.002, 8001)
    r = BoxLeastSquares(T, F, E).power(fine, [0.03, 0.04, 0.05, 0.06, 0.08], objective="snr")
    i = int(np.argmax(r.power))
    P2, T2, dur = float(fine[i]), float(r.transit_time[i]), float(r.duration[i])
    n_mid = round((np.median(T) - T2) / P2)
    T2 = T2 + n_mid * P2                                      # reference epoch near the middle of the data
    out["bls_joint"] = {"P": P2, "t0_btjd": T2, "dur_h": dur * 24, "depth_ppm": float(r.depth[i] * 1e6),
                        "depth_err_ppm": float(r.depth_err[i] * 1e6), "snr": float(r.power[i])}
    print("joint", out["bls_joint"], flush=True)

    # 3. per-sector and per-epoch-group depths
    sec = []
    for s, t, f, e, cad in rows:
        d, de, nin = depth_with_error(t, f, P2, T2, dur, rng)
        n_tr = len(np.unique(np.round((t[phase_dist(t, P2, T2) < dur / 2] - T2) / P2)))
        sec.append({"sector": int(s), "depth_ppm": d, "err_ppm": de, "n_in": nin, "n_transits": n_tr})
    out["sectors_depth"] = sec
    d_all, e_all, _ = depth_with_error(T, F, P2, T2, dur, rng, 1000)
    out["depth_all"] = [d_all, e_all]
    groups = {}
    ep = np.round((T - T2) / P2).astype(int)
    for k in (2, 3, 4, 5):
        groups[k] = []
        for j in range(k):
            sel = (ep % k) == j
            d, de, _ = depth_with_error(T[sel], F[sel], P2, T2, dur, rng)
            groups[k].append([d, de])
    out["epoch_groups"] = groups
    # 4. secondary eclipse
    d05, e05, _ = depth_with_error(T, F, P2, T2 + P2 / 2, dur, rng, 1000)
    scan = [(ph, box_depth(T, F, P2, T2 + ph * P2, dur)[0]) for ph in np.arange(0.3, 0.7001, 0.002)]
    best = max(scan, key=lambda x: x[1])
    out["secondary"] = {"phase0.5_ppm": [d05, e05], "deepest_0.3_0.7": {"phase": best[0], "depth_ppm": best[1],
                                                                       "sigma": best[1] / e05}}
    print("sectors", [(x["sector"], round(x["depth_ppm"]), round(x["err_ppm"])) for x in sec], flush=True)
    print("all", d_all, e_all, "groups", groups, "secondary", out["secondary"], flush=True)

    # 5. transit fit on the folded light curve (period fixed at the joint BLS value; a free period lets the walkers
    #    hop between ephemerides one cycle apart over the ~650 cycles), stellar density free
    fit, ch = fit_folded(rows, P2, T2, nstep)
    out["fit_free_rho"] = fit
    rho = 10 ** ch[:, 3]
    out["rho_compare"] = {"p_rho_ge_catalogue_minus_1sigma": float(np.mean(rho > RHO_STAR[0] - RHO_STAR[1])),
                          "rho_star_catalogue": RHO_STAR}
    k = ch[:, 1]
    out["radius_earth_if_on_target"] = [float(np.median(k) * 0.817 * 109.1), float(np.std(k) * 0.817 * 109.1)]
    # 6. neighbour
    out["neighbour"] = {"flux_fraction_G": NEIGH_FRAC,
                        "required_eclipse_depth": d_all * 1e-6 / NEIGH_FRAC}
    json.dump(out, open(OUT, "w"), indent=1, default=float)
    print(json.dumps({k: v for k, v in out.items() if k not in ("searches",)}, indent=1, default=float)[:4000])
    plot(T, F, P2, T2, dur, sec, d_all)


def fit_folded(rows, P, t0, nstep, nwalk=32, seed=2):
    """Folded-light-curve transit fit: mid-time offset (min), Rp/R*, b, log10 rho*; per-sector white-noise errors
    times the sector's red-noise factor in 1-hour bins; each cadence's exposure time integrated."""
    import emcee
    from scipy.optimize import minimize
    u = tf.tess_ld(TEFF)
    D = []
    for s, t, f, e, cad in rows:
        x = ((t - t0 + P / 2) % P) - P / 2
        m = np.abs(x) < 0.25
        if m.sum() < 30:
            continue
        r = f - 1
        sig = 1.4826 * np.median(np.abs(np.diff(r))) / np.sqrt(2)
        oot = np.abs(x) > 0.08
        nb = max(int(round((1 / 24) / (cad / 86400))), 1)
        rr = r[oot][: r[oot].size // nb * nb].reshape(-1, nb).mean(axis=1)
        beta = max(1.0, float(1.4826 * np.median(np.abs(rr - np.median(rr))) / (sig / np.sqrt(nb))))
        D.append((x[m], f[m], sig * beta, cad))

    def lnp(th):
        dt, k, b, lrho = th
        if not (-60 < dt < 60 and 0.002 < k < 0.3 and 0 <= b < 1 + k and -2 < lrho < 2):
            return -np.inf
        aR = ((10 ** lrho) * tf.RHO_SUN * tf.G * (P * 86400) ** 2 / (3 * np.pi)) ** (1 / 3)
        if aR <= 1 + k or b >= aR:
            return -np.inf
        inc = np.degrees(np.arccos(b / aR))
        return sum(-0.5 * np.sum(((f - tf._model(x, dt / 1440, P, k, aR, inc, u, exp_s=cad)) / sg) ** 2)
                   for x, f, sg, cad in D)

    best = None
    for x0 in ([0, 0.02, 0.3, 0.2], [0, 0.022, 0.85, 0.2], [0, 0.03, 0.95, 0.0], [0, 0.02, 0.6, -0.5]):
        r = minimize(lambda p: -lnp(p), x0, method="Nelder-Mead", options={"maxiter": 8000, "xatol": 1e-6, "fatol": 1e-3})
        if best is None or r.fun < best.fun:
            best = r
    rng = np.random.default_rng(seed)
    pos = best.x + np.array([1.0, 5e-4, 0.02, 0.02]) * rng.normal(size=(nwalk, 4))
    pos[:, 2] = np.abs(pos[:, 2])
    sampler = emcee.EnsembleSampler(nwalk, 4, lnp)
    sampler.run_mcmc(pos, nstep, progress=False)
    ch = sampler.get_chain(discard=nstep // 3, thin=5, flat=True)
    dt, k, b, lrho = ch.T
    aR = ((10 ** lrho) * tf.RHO_SUN * tf.G * (P * 86400) ** 2 / (3 * np.pi)) ** (1 / 3)
    t14 = tf.t14_hours(P, aR, k, b)
    out = {"P_fixed": P, "t0_btjd": t0, "u": u, "n_points": int(sum(len(d[0]) for d in D)),
           "chi2_min": float(2 * best.fun), "sigma_per_cadence_scaled": [float(d[2]) for d in D],
           "chain_summary": tf.convergence(sampler, nstep, ["dt", "rp_rs", "b", "log10_rho"])}
    for name, v in (("dt_min", dt), ("rp_rs", k), ("b", b), ("log10_rho", lrho), ("rho_solar", 10 ** lrho),
                    ("a_rs", aR), ("t14_h", t14)):
        q = np.percentile(v, [16, 50, 84])
        out[name] = [float(q[1]), float(q[2] - q[1]), float(q[1] - q[0])]
    out["b_5"] = float(np.percentile(b, 5))
    out["rho_95"] = float(np.percentile(10 ** lrho, 95))
    med = np.median(ch, axis=0)
    aRm = ((10 ** med[3]) * tf.RHO_SUN * tf.G * (P * 86400) ** 2 / (3 * np.pi)) ** (1 / 3)
    out["model_depth_ppm"] = float((1 - tf._model(np.array([0.0]), 0.0, P, med[1], aRm,
                                                   np.degrees(np.arccos(med[2] / aRm)), u, exp_s=120.0)[0]) * 1e6)
    print("folded fit", {k_: v for k_, v in out.items() if k_ != "chain_summary"}, out["chain_summary"]["steps_over_tau"], flush=True)
    return out, ch


def plot(T, F, P, t0, dur, sec, d_all):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ph = ((T - t0 + P / 2) % P) - P / 2
    m = np.abs(ph) < 0.25
    x, y = ph[m] * 24, (F[m] - 1) * 1e6
    edges = np.arange(-6, 6.01, 0.25)
    idx = np.digitize(x, edges)
    xb = np.array([x[idx == j].mean() for j in range(1, len(edges)) if (idx == j).sum() > 5])
    yb = np.array([y[idx == j].mean() for j in range(1, len(edges)) if (idx == j).sum() > 5])
    eb = np.array([y[idx == j].std() / np.sqrt((idx == j).sum()) for j in range(1, len(edges)) if (idx == j).sum() > 5])
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.4))
    ax[0].errorbar(xb, yb, eb, fmt="o", ms=3, color="k")
    ax[0].axhline(0, color="0.6", lw=0.8)
    ax[0].set_xlabel("hours from mid-transit")
    ax[0].set_ylabel("relative flux (ppm)")
    ax[0].set_title(f"TIC 402710280, P = {P:.5f} d, 15-min bins")
    s = [x["sector"] for x in sec]
    ax[1].errorbar(range(len(s)), [x["depth_ppm"] for x in sec], [x["err_ppm"] for x in sec], fmt="o", color="k")
    ax[1].axhline(d_all, color="C3", lw=1)
    ax[1].axhline(0, color="0.6", lw=0.8)
    ax[1].set_xticks(range(len(s)))
    ax[1].set_xticklabels([f"S{v}" for v in s])
    ax[1].set_ylabel("box depth (ppm)")
    fig.tight_layout()
    fig.savefig(FIG, dpi=130)


def refit(nstep):
    """Re-run only the transit fit (step 5) on the stored ephemeris, e.g. with a longer chain."""
    out = json.load(open(OUT))
    P2, T2 = out["bls_joint"]["P"], out["bls_joint"]["t0_btjd"]
    _, rows = load(P2, T2)
    fit, ch = fit_folded(rows, P2, T2, nstep)
    out["fit_free_rho"] = fit
    rho = 10 ** ch[:, 3]
    out["rho_compare"] = {"p_rho_below_catalogue_plus_2sigma": float(np.mean(rho < RHO_STAR[0] + 2 * RHO_STAR[1])),
                          "rho_star_catalogue": RHO_STAR}
    k = ch[:, 1]
    out["radius_earth_if_on_target"] = [float(np.median(k) * 0.817 * 109.1), float(np.std(k) * 0.817 * 109.1)]
    np.save(RESULTS / "candidates" / "second_signal_chain.npy", ch.astype(np.float32))
    json.dump(out, open(OUT, "w"), indent=1, default=float)


if __name__ == "__main__":
    if "--refit" in sys.argv:
        refit(int(sys.argv[-1]))
    else:
        main(int(sys.argv[1]) if len(sys.argv) > 1 else 12000)
