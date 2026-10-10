"""Joint transit fits of a few individual transits in SPOC 2-minute data (batman + emcee).

Each transit window (+-`win` days) is normalised by a quadratic fitted to its out-of-transit part, then all
windows are fitted together with one transit shape. Two parameterisations:
  - periodic: t0 (reference transit), P, Rp/R*, b, log10(rho*/rho_sun); a/R* follows from rho* and P
  - shape-only (period unknown, e.g. a duo): Rp/R*, b, T14 and one mid-time per transit; the orbit is
    a long circular one whose a/R* reproduces T14, so only the shape is constrained
Limb darkening is fixed at quadratic coefficients for the TESS band (Claret 2017, interpolated coarsely
in Teff). Output: posterior summaries and binned data + model curves for plotting.
"""
import io

import env  # noqa: F401
import batman
import emcee
import numpy as np
from astropy.io import fits
from scipy.optimize import minimize

import s3io

G = 6.674e-11
RHO_SUN = 1408.0  # kg m^-3


def tess_ld(teff):
    """Rough TESS-band quadratic limb darkening (Claret 2017, logg 4.5, solar metallicity)."""
    grid = np.array([3500, 4500, 5500, 6000, 6500, 7000])
    u1 = np.array([0.20, 0.42, 0.38, 0.32, 0.28, 0.25])
    u2 = np.array([0.38, 0.20, 0.21, 0.22, 0.24, 0.26])
    return float(np.interp(teff, grid, u1)), float(np.interp(teff, grid, u2))


def load_windows(tic, events, win, t14_h):
    """events: list of (sector, tc_guess). Returns list of (t, f, sigma) normalised windows."""
    out = []
    for s, tc in events:
        tm = s3io.spoc_templates()[str(s)]
        key = f"tess/public/tid/s{s:04d}/{s3io.tic_dir(tic)}/{tm[0]}-s{s:04d}-{tic:016d}-{tm[1]}-s_lc.fits"
        d = fits.open(io.BytesIO(s3io.fetch(key)))[1].data
        t = np.asarray(d["TIME"], float)
        f = np.asarray(d["PDCSAP_FLUX"], float)
        q = np.asarray(d["QUALITY"])
        ok = np.isfinite(t) & np.isfinite(f) & (q == 0) & (np.abs(t - tc) < win)
        t, f = t[ok], f[ok]
        oot = np.abs(t - tc) > 0.75 * t14_h / 24
        c = np.polyfit(t[oot] - tc, f[oot], 2)
        f = f / np.polyval(c, t - tc)
        r = f[oot] - 1
        sig = 1.4826 * np.median(np.abs(np.diff(r))) / np.sqrt(2)
        out.append((t, f, sig))
    return out


def _model(t, tc, P, rp, aR, inc, u, exp_s=120.0):
    pm = batman.TransitParams()
    pm.t0, pm.per, pm.rp, pm.a, pm.inc = tc, P, rp, aR, inc
    pm.ecc, pm.w, pm.u, pm.limb_dark = 0.0, 90.0, list(u), "quadratic"
    return batman.TransitModel(pm, t, supersample_factor=5, exp_time=exp_s / 86400).light_curve(pm)


def t14_hours(P, aR, rp, b):
    inc = np.arccos(b / aR)
    arg = np.sqrt(np.clip((1 + rp) ** 2 - b ** 2, 0, None)) / (aR * np.sin(inc))
    return P / np.pi * np.arcsin(np.clip(arg, 0, 1)) * 24


def fit_periodic(tic, events, epochs, P0, t14_h, teff, nwalk=40, nstep=4000, win=1.0, seed=1):
    """Transits at epochs (integers) of one planet; events[i] is (sector, tc_guess) for epochs[i]."""
    data = load_windows(tic, events, win, t14_h)
    u = tess_ld(teff)
    epochs = np.asarray(epochs)
    ref = int(np.round(np.median(epochs)))
    tcs0 = np.array([tc for _, tc in events])
    if len(epochs) >= 2:
        P0 = float(np.polyfit(epochs, tcs0, 1)[0])
    tc_ref0 = float(np.mean(tcs0 - (epochs - ref) * P0))

    def unpack(th):
        tref, P, rp, b, lrho = th
        aR = ((10 ** lrho) * RHO_SUN * G * (P * 86400) ** 2 / (3 * np.pi)) ** (1 / 3)
        return tref, P, rp, b, aR

    def lnp(th):
        tref, P, rp, b, lrho = th
        if not (0.002 < rp < 0.3 and 0 <= b < 1 + rp and -3 < lrho < 2 and abs(P - P0) < 0.05):
            return -np.inf
        _, _, _, _, aR = unpack(th)
        if aR <= 1 + rp or b >= aR:
            return -np.inf
        inc = np.degrees(np.arccos(b / aR))
        ll = 0.0
        for (t, f, sig), ep in zip(data, epochs):
            m = _model(t, tref + (ep - ref) * P, P, rp, aR, inc, u)
            ll += -0.5 * np.sum(((f - m) / sig) ** 2)
        return ll

    rho_guess = np.log10(0.5)
    x0 = np.array([tc_ref0, P0, 0.025, 0.4, rho_guess])
    res = minimize(lambda p: -lnp(p), x0, method="Nelder-Mead",
                   options={"maxiter": 40000, "xatol": 1e-8, "fatol": 1e-4})
    rng = np.random.default_rng(seed)
    scale = np.array([5e-4, 5e-6, 5e-4, 0.02, 0.02])
    pos = res.x + scale * rng.normal(size=(nwalk, 5))
    pos[:, 3] = np.abs(pos[:, 3])
    sampler = emcee.EnsembleSampler(nwalk, 5, lnp)
    sampler.run_mcmc(pos, nstep, progress=False)
    ch = sampler.get_chain(discard=nstep // 3, thin=5, flat=True)
    tref, P, rp, b, lrho = ch.T
    aR = unpack(ch.T)[4]
    t14 = t14_hours(P, aR, rp, b)
    out = {"ref_epoch": ref, "u": u, "n_points": int(sum(len(t) for t, _, _ in data)),
           "chi2_min": float(-2 * lnp(np.median(ch, axis=0)))}
    for name, v in (("tref_btjd", tref), ("period_d", P), ("rp_rs", rp), ("b", b), ("log10_rho", lrho),
                    ("rho_solar", 10 ** lrho), ("a_rs", aR), ("t14_h", t14), ("depth_ppm", rp ** 2 * 1e6)):
        q = np.percentile(v, [16, 50, 84])
        out[name] = [float(q[1]), float(q[2] - q[1]), float(q[1] - q[0])]
    med = np.median(ch, axis=0)
    out["plot"] = plot_data(data, [med[0] + (ep - ref) * med[1] for ep in epochs], med[1], med[2],
                            unpack(med)[4], np.degrees(np.arccos(med[3] / unpack(med)[4])), u)
    out["chain_summary"] = {"acceptance": float(np.mean(sampler.acceptance_fraction))}
    return out, ch


def fit_shape(tic, events, t14_h, teff, nwalk=40, nstep=4000, win=0.6, seed=1, P_dummy=300.0):
    """Shape-only fit for transits whose period is unknown: Rp/R*, b, T14, and each mid-time."""
    data = load_windows(tic, events, win, t14_h)
    u = tess_ld(teff)
    n = len(data)

    def aR_from(rp, b, T14):
        return np.sqrt(max((1 + rp) ** 2 - b ** 2, 1e-12)) / np.sin(np.pi * T14 / 24 / P_dummy)

    def lnp(th):
        rp, b, T14 = th[:3]
        if not (0.002 < rp < 0.3 and 0 <= b < 1 + rp and 0.5 < T14 < 30):
            return -np.inf
        aR = aR_from(rp, b, T14)
        if b >= aR:
            return -np.inf
        inc = np.degrees(np.arccos(b / aR))
        ll = 0.0
        for i, (t, f, sig) in enumerate(data):
            ll += -0.5 * np.sum(((f - _model(t, th[3 + i], P_dummy, rp, aR, inc, u)) / sig) ** 2)
        return ll

    x0 = np.array([0.025, 0.4, t14_h] + [tc for _, tc in events])
    res = minimize(lambda p: -lnp(p), x0, method="Nelder-Mead",
                   options={"maxiter": 40000, "xatol": 1e-8, "fatol": 1e-4})
    rng = np.random.default_rng(seed)
    ndim = 3 + n
    scale = np.array([5e-4, 0.02, 0.02] + [5e-4] * n)
    pos = res.x + scale * rng.normal(size=(nwalk, ndim))
    pos[:, 1] = np.abs(pos[:, 1])
    sampler = emcee.EnsembleSampler(nwalk, ndim, lnp)
    sampler.run_mcmc(pos, nstep, progress=False)
    ch = sampler.get_chain(discard=nstep // 3, thin=5, flat=True)
    out = {"u": u, "n_points": int(sum(len(t) for t, _, _ in data))}
    names = ["rp_rs", "b", "t14_h"] + [f"tc_{i}" for i in range(n)]
    for i, name in enumerate(names):
        q = np.percentile(ch[:, i], [16, 50, 84])
        out[name] = [float(q[1]), float(q[2] - q[1]), float(q[1] - q[0])]
    q = np.percentile(ch[:, 0] ** 2 * 1e6, [16, 50, 84])
    out["depth_ppm"] = [float(q[1]), float(q[2] - q[1]), float(q[1] - q[0])]
    med = np.median(ch, axis=0)
    aR = aR_from(med[0], med[1], med[2])
    out["plot"] = plot_data(data, list(med[3:]), P_dummy, med[0], aR, np.degrees(np.arccos(med[1] / aR)), u)
    # per-transit depth with the shape fixed at the median
    per = []
    for i, (t, f, sig) in enumerate(data):
        def chi(r):
            return np.sum(((f - _model(t, med[3 + i], P_dummy, r[0], aR, np.degrees(np.arccos(med[1] / aR)), u)) / sig) ** 2)
        r0 = minimize(chi, [med[0]], method="Nelder-Mead").x[0]
        grid = np.linspace(r0 * 0.85, r0 * 1.15, 121)
        c = np.array([chi([g]) for g in grid])
        okk = grid[c < c.min() + 1]
        per.append({"rp_rs": float(r0), "rp_rs_lo": float(okk.min()), "rp_rs_hi": float(okk.max()),
                    "depth_ppm": float(r0 ** 2 * 1e6)})
    out["per_transit"] = per
    return out, ch


def plot_data(data, tcs, P, rp, aR, inc, u, bin_min=20.0):
    """Binned points (hours from mid-transit, ppm) and model curves for each window."""
    panels = []
    for (t, f, sig), tc in zip(data, tcs):
        x = (t - tc) * 24
        edges = np.arange(x.min(), x.max() + bin_min / 60, bin_min / 60)
        k = np.digitize(x, edges)
        xb, yb, eb = [], [], []
        for j in np.unique(k):
            m = k == j
            if m.sum() < 3:
                continue
            xb.append(float(x[m].mean()))
            yb.append(float((f[m].mean() - 1) * 1e6))
            eb.append(float(sig / np.sqrt(m.sum()) * 1e6))
        tm = np.linspace(x.min(), x.max(), 400)
        mm = _model(tc + tm / 24, tc, P, rp, aR, inc, u)
        panels.append({"tc": float(tc), "x": np.round(xb, 3).tolist(), "y": np.round(yb, 1).tolist(),
                       "e": np.round(eb, 1).tolist(), "mx": np.round(tm, 3).tolist(),
                       "my": np.round((mm - 1) * 1e6, 1).tolist(), "sigma_ppm_2min": float(sig * 1e6)})
    return panels
