"""Track A, triage v2: which flag-free signals are already known, and how many unknown ones beat the flipped-data
false-alarm baseline.

The flipped light curves (flux mirrored about its median) went through the same search and the same flags, star for
star, so a flag-free flipped signal is a false alarm of the kind the real search also makes: its count above a
threshold estimates how many real flag-free signals above that threshold are noise.

Known: a TOI (any disposition), a community TOI or a NASA confirmed planet on the same star (TIC, or within 10"
for NASA) at the same period or a 1/3-3x multiple of it (0.5%); or a TOI within 2' at the same period (a blend).

  python tracka_v2.py
Output: results/track_a/tracka_v2.json, results/track_a/tracka_v2_unknown.csv
"""
import json

import env  # noqa: F401
import numpy as np
import pandas as pd

from env import CATALOGS, RESULTS

A = RESULTS / "track_a"
RATIOS = [1, 2, 3, 1 / 2, 1 / 3, 2 / 3, 3 / 2]


def same_period(p, q, tol=5e-3):
    return any(abs(p / (q * r) - 1) < tol for r in RATIOS)


def catalogs():
    toi = pd.read_csv(CATALOGS / "exofop" / "toi_20261012.csv", low_memory=False)
    toi = toi.assign(P=pd.to_numeric(toi["Period (days)"], errors="coerce"), tic=toi["TIC ID"].astype(int),
                     disp=toi["TFOPWG Disposition"].fillna(toi["TESS Disposition"]).fillna(""))
    ctoi = pd.read_csv(CATALOGS / "exofop" / "ctoi_20261012.csv", low_memory=False)
    ctoi = ctoi.assign(P=pd.to_numeric(ctoi["Period (days)"], errors="coerce"), tic=ctoi["TIC ID"].astype(int))
    ps = pd.read_csv(CATALOGS / "nasa_pscomppars_transiting_20261012.csv")
    return toi, ctoi, ps


def sep_arcsec(ra1, de1, ra2, de2):
    ra1, de1, ra2, de2 = map(np.radians, (ra1, de1, ra2, de2))
    c = np.sin(de1) * np.sin(de2) + np.cos(de1) * np.cos(de2) * np.cos(ra1 - ra2)
    return np.degrees(np.arccos(np.clip(c, -1, 1))) * 3600


def known(r, pos, toi, ctoi, ps):
    out = []
    t = toi[toi.tic == r.tic]
    for x in t.itertuples():
        if x.P == x.P and same_period(r.period, x.P):
            out.append(f"TOI-{x.TOI} ({x.disp}, P={x.P:.4f})")
    for x in ctoi[ctoi.tic == r.tic].itertuples():
        if x.P == x.P and same_period(r.period, x.P):
            out.append(f"CTOI {x.CTOI} (P={x.P:.4f})")
    ra, de = pos
    s = sep_arcsec(ra, de, ps.ra.values, ps.dec.values)
    for x in ps[s < 10].itertuples():
        if x.pl_orbper == x.pl_orbper and same_period(r.period, x.pl_orbper):
            out.append(f"{x.pl_name} (P={x.pl_orbper:.4f})")
    if not out:
        s = sep_arcsec(ra, de, toi.RA.values, toi.Dec.values) if toi.RA.dtype != object else None
        if s is not None:
            for x, d in zip(toi[s < 120].itertuples(), s[s < 120]):
                if x.tic != r.tic and x.P == x.P and same_period(r.period, x.P, 2e-3):
                    out.append(f"blend: TOI-{x.TOI} ({x.disp}) {d:.0f}\" away")
    other = [f"TOI-{x.TOI} (P={x.P:.3f})" for x in t.itertuples() if x.P == x.P and not same_period(r.period, x.P)]
    return out, other


def radec_toi(toi):
    """ExoFOP gives RA/Dec as sexagesimal strings; convert to degrees."""
    from astropy.coordinates import Angle
    import astropy.units as u
    if toi.RA.dtype == object:
        toi = toi.copy()
        toi["RA"] = Angle(toi.RA.astype(str).values, unit=u.hourangle).deg
        toi["Dec"] = Angle(toi.Dec.astype(str).values, unit=u.deg).deg
    return toi


def main():
    tri = pd.read_csv(A / "triage_v2.csv", low_memory=False)
    tg = pd.read_csv(A / "targets_search.csv").set_index("tic")
    toi, ctoi, ps = catalogs()
    toi = radec_toi(toi)
    free = tri[tri["hard_flags"].fillna("") == ""].copy()
    rows = []
    for r in free.itertuples():
        pos = (tg.at[r.tic, "ra"], tg.at[r.tic, "dec"]) if r.tic in tg.index else (np.nan, np.nan)
        k, other = known(r, pos, toi, ctoi, ps)
        st = tg.loc[r.tic] if r.tic in tg.index else None
        rp = float(np.sqrt(max(r.depth_ppm, 0) * 1e-6) * (st.rad if st is not None else 0.5) * 109.1)
        rows.append({"tic": r.tic, "kind": r.kind, "period": r.period, "t0": r.t0, "duration_hr": r.duration_hr,
                     "depth_ppm": r.depth_ppm, "rp_rearth": rp, "sde": r.sde, "snr": r.snr, "snr_emp": r.snr_emp,
                     "n_tr_data": r.n_tr_data, "soft_flags": r._asdict().get("flags"),
                     "Tmag": st.Tmag if st is not None else np.nan, "Teff": st.Teff if st is not None else np.nan,
                     "rad": st.rad if st is not None else np.nan, "known": "; ".join(k), "other_tois": "; ".join(other)})
    df = pd.DataFrame(rows)
    real = df[df.kind == "signals"]
    flip = df[df.kind == "signals_inverted"]
    unk = real[real.known == ""]
    # false-alarm calibration on SNR against red noise, and on SDE
    cal = {}
    for col in ("snr_emp", "sde"):
        rows_c = []
        for thr in sorted(set(np.round(np.quantile(flip[col], np.linspace(0, 1, 41)), 2))):
            rows_c.append({"thr": float(thr), "flipped": int((flip[col] >= thr).sum()),
                           "real_unknown": int((unk[col] >= thr).sum()), "real_known": int(((real.known != "") & (real[col] >= thr)).sum())})
        cal[col] = rows_c
    # the threshold that leaves at most one flipped signal, jointly on both
    def n_flip(a, b):
        return int(((flip.snr_emp >= a) & (flip.sde >= b)).sum())
    best = None
    for a in np.arange(6, 40, 0.5):
        for b in np.arange(7, 40, 0.5):
            if n_flip(a, b) <= 1:
                n = int(((unk.snr_emp >= a) & (unk.sde >= b)).sum())
                if best is None or n > best[2]:
                    best = (float(a), float(b), n)
                break
    a, b, _ = best
    sel = unk[(unk.snr_emp >= a) & (unk.sde >= b)].sort_values("sde", ascending=False)
    out = {"n_free_real": len(real), "n_free_flipped": len(flip), "n_known": int((real.known != "").sum()),
           "n_unknown": len(unk), "threshold": {"snr_emp": a, "sde": b, "flipped_above": n_flip(a, b),
                                                "unknown_above": len(sel)},
           "calibration": cal, "known": real[real.known != ""][["tic", "period", "sde", "known"]].to_dict("records"),
           "unknown_above": sel.to_dict("records")}
    json.dump(out, open(A / "tracka_v2.json", "w"), indent=1, default=float)
    unk.sort_values("sde", ascending=False).to_csv(A / "tracka_v2_unknown.csv", index=False)
    print("flag-free real", len(real), "flipped", len(flip), "known", out["n_known"], "unknown", len(unk))
    print("threshold", out["threshold"])
    print(sel[["tic", "period", "depth_ppm", "rp_rearth", "sde", "snr_emp", "n_tr_data", "Tmag", "Teff", "soft_flags", "other_tois"]].to_string())


if __name__ == "__main__":
    main()
