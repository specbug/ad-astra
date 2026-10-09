"""Track A target selection: K/M dwarfs within 150 pc whose TESS data were never jointly searched.

A star qualifies when it has >= 3 sectors spanning >= 300 d in total, at most two of them at 2-minute
cadence (so SPOC's 2-min multi-sector search saw at most two sectors, and nobody's search combined the
2-minute and full-frame sectors), low contamination, and is not already a TOI host. Stars are ranked by
the expected signal-to-noise of a 1.5 Earth-radius planet on a 10-day orbit, using the measured QLP noise.
"""
import numpy as np
import pandas as pd

import env  # noqa: F401
from env import CATALOGS, RESULTS

# 1-hour box-average noise (ppm) of QLP 200-s light curves vs Tmag, measured in the sector-80 pilot
NOISE_T = np.array([7.25, 8.25, 8.75, 9.25, 9.75, 10.25, 10.75, 11.25, 11.75, 12.25, 12.75, 13.25])
NOISE_PPM = np.array([218, 196, 244, 273, 335, 440, 610, 881, 1058, 1623, 2543, 4014])


def noise_1h(tmag):
    return np.interp(tmag, NOISE_T, NOISE_PPM) * 1e-6


def main():
    tic = pd.read_parquet(CATALOGS / "tic_master.parquet")
    s = tic[(tic.objType == "STAR") & (tic.d < 150) & (tic.Tmag < 13.5) & (tic.Teff.between(2700, 5400))
            & (tic.rad < 1.0) & (tic.logg > 4.2) & (tic.contratio.fillna(0) < 0.2)
            & ~tic.disposition.isin(["DUPLICATE", "ARTIFACT", "SPLIT"])].copy()
    qlp = pd.read_parquet(CATALOGS / "qlp_manifest.parquet")
    tsp = pd.read_parquet(CATALOGS / "tspoc_manifest.parquet")
    spc = pd.read_parquet(CATALOGS / "spoc_manifest.parquet")
    allsec = pd.concat([qlp[["tic", "sector"]], tsp[["tic", "sector"]], spc[["tic", "sector"]]]).drop_duplicates()
    allsec = allsec[allsec.tic.isin(s.tic)]
    g = allsec.groupby("tic").sector.agg(n_sec="nunique", s_min="min", s_max="max")
    n2 = spc[spc.tic.isin(s.tic)].groupby("tic").sector.nunique().rename("n_2min")
    s = s.join(g, on="tic").join(n2, on="tic")
    s["n_2min"] = s.n_2min.fillna(0).astype(int)
    s = s[s.n_sec >= 3]
    # approximate span: sectors are ~27.4 d apart
    s["span_d"] = (s.s_max - s.s_min) * 27.4
    s = s[s.span_d >= 300]
    toi = pd.read_csv(CATALOGS / "known" / "exofop_toi_autottv_2025.csv", usecols=["TIC ID"])
    s["toi_host"] = s.tic.isin(toi["TIC ID"])
    # expected SNR of a 1.5 Re planet, P = 10 d, central transit
    rs = s.rad.clip(0.08, 1.0)
    depth = (1.5 / 109.1 / rs) ** 2 / (1 + s.contratio.fillna(0))
    rho = (s.mass / rs ** 3).fillna(1.0)
    dur_h = 13.0 * (10 / 365.25) ** (1 / 3) * rho ** (-1 / 3)
    n_tr = s.n_sec * 27.4 * 0.85 / 10
    s["exp_snr"] = depth / (noise_1h(s.Tmag) / np.sqrt(dur_h)) * np.sqrt(n_tr)
    s["tier"] = np.where(s.n_2min == 0, "ffi_only", "few_2min")
    s = s[(s.n_2min <= 2) & ~s.toi_host].sort_values("exp_snr", ascending=False)
    out = RESULTS / "track_a"
    out.mkdir(parents=True, exist_ok=True)
    cols = ["tic", "ra", "dec", "Tmag", "Teff", "rad", "mass", "d", "contratio", "n_sec", "n_2min", "span_d",
            "exp_snr", "tier"]
    s[cols].to_csv(out / "targets_all.csv", index=False, float_format="%.5g")
    print("candidates:", len(s), s.tier.value_counts().to_dict())
    for thr in (5, 7, 10, 15):
        print(f"exp_snr >= {thr}:", int((s.exp_snr >= thr).sum()))
    print(s[cols].head(10).to_string())


if __name__ == "__main__":
    main()
