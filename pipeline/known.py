"""Known-object cross-matching: is a star already a TOI / known planet host / in another list?

Sources available offline in this environment:
  - ExoFOP TOI table export (Feb 2026; TIC IDs, periods, epochs, dispositions), from the public
    shporer/autottv repository
  - NASA Exoplanet Archive TOI snapshot of 5 Oct 2026 (positions, periods; copied from the
    juliensimon/tess-toi-candidates Hugging Face dataset), when present
Matching is by TIC ID and, for the snapshot without TIC IDs, by position (10 arcsec).
"""
import numpy as np
import pandas as pd

from env import CATALOGS

_cache = {}


def exofop_toi():
    if "exofop" not in _cache:
        d = pd.read_csv(CATALOGS / "known" / "exofop_toi_autottv_2025.csv")
        d = d.rename(columns={"TIC ID": "tic", "Period (days)": "period", "Epoch (BJD)": "epoch_bjd",
                              "Depth (ppm)": "depth_ppm", "Duration (hours)": "dur_h",
                              "TFOPWG Disposition": "tfop", "TESS Disposition": "tess_disp", "Comments": "comments"})
        _cache["exofop"] = d
    return _cache["exofop"]


def hf_toi():
    if "hf" not in _cache:
        p = CATALOGS / "hf_toi_20261005.csv"
        _cache["hf"] = pd.read_csv(p) if p.exists() else None
    return _cache["hf"]


def sep_arcsec(ra1, dec1, ra2, dec2):
    ra1, dec1, ra2, dec2 = map(np.radians, (ra1, dec1, ra2, dec2))
    c = np.sin(dec1) * np.sin(dec2) + np.cos(dec1) * np.cos(dec2) * np.cos(ra1 - ra2)
    return np.degrees(np.arccos(np.clip(c, -1, 1))) * 3600


def lookup(tic, ra=None, dec=None, radius_arcsec=10.0):
    """All known entries for one star: list of dicts."""
    out = []
    e = exofop_toi()
    for _, r in e[e.tic == int(tic)].iterrows():
        out.append({"source": "ExoFOP TOI (Feb 2026)", "id": f"TOI-{r['TOI']}", "period": r["period"],
                    "depth_ppm": r["depth_ppm"], "dur_h": r["dur_h"], "disp": f"{r['tfop']}/{r['tess_disp']}",
                    "comments": r["comments"]})
    h = hf_toi()
    if h is not None and ra is not None:
        s = sep_arcsec(ra, dec, h.ra_deg.values, h.dec_deg.values)
        for _, r in h[s < radius_arcsec].iterrows():
            out.append({"source": "NASA TOI table (5 Oct 2026)", "id": f"TOI-{r['toi_prefix']}",
                        "period": r.get("period_days"), "depth_ppm": r.get("transit_depth_ppm"),
                        "disp": r.get("disposition")})
    return out


def flag_known(df, ra_col="ra", dec_col="dec"):
    """Add is_toi / toi_ids columns to a DataFrame with tic (and optionally ra/dec)."""
    e = exofop_toi()
    ids = e.groupby("tic").TOI.apply(lambda x: ",".join(str(v) for v in x))
    df = df.copy()
    df["toi_exofop"] = df.tic.map(ids)
    h = hf_toi()
    if h is not None and ra_col in df:
        from scipy.spatial import cKDTree
        def xyz(ra, dec):
            ra, dec = np.radians(ra), np.radians(dec)
            return np.c_[np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)]
        tree = cKDTree(xyz(h.ra_deg.values, h.dec_deg.values))
        d, i = tree.query(xyz(df[ra_col].values, df[dec_col].values), k=1)
        sep = 2 * np.degrees(np.arcsin(np.clip(d / 2, 0, 1))) * 3600
        df["toi_hf"] = np.where(sep < 10, h.toi_prefix.values[np.clip(i, 0, len(h) - 1)].astype(str), None)
    else:
        df["toi_hf"] = None
    df["is_toi"] = df.toi_exofop.notna() | df.toi_hf.notna()
    return df
