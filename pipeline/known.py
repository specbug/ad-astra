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


def confirmed_planets():
    """NASA Exoplanet Archive confirmed planets (ps table, default_flag=1) of 5 Oct 2026, copied page by page from
    the juliensimon/nasa-exoplanets Hugging Face mirror into catalogs/hf_ps/; None if absent."""
    if "ps" not in _cache:
        pages = sorted((CATALOGS / "hf_ps").glob("page_*.csv"))
        _cache["ps"] = (pd.concat([pd.read_csv(f) for f in pages]).drop_duplicates("row").sort_values("row")
                        .reset_index(drop=True) if pages else None)
    return _cache["ps"]


def flag_confirmed(df, ra_col="ra", dec_col="dec", radius_arcsec=10.0):
    """Add known_planet / known_planet_P / known_planet_sep columns: confirmed planets whose host lies within
    radius_arcsec of the star, or whose host is named 'TIC <tic>'. Positions are compared as catalogued (the TIC
    at epoch 2000, the archive mostly at the Gaia epoch), so fast-moving stars rely on the name match."""
    ps = confirmed_planets()
    df = df.copy()
    df["known_planet"] = None
    df["known_planet_P"] = None
    df["known_planet_sep"] = np.nan
    if ps is None:
        return df
    from scipy.spatial import cKDTree

    def xyz(ra, dec):
        ra, dec = np.radians(ra), np.radians(dec)
        return np.c_[np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)]
    tree = cKDTree(xyz(ps.ra.values, ps.dec.values))
    ok = df[ra_col].notna() & df[dec_col].notna()
    chord = 2 * np.sin(np.radians(radius_arcsec / 3600) / 2)
    pos = xyz(df[ra_col].fillna(0).values, df[dec_col].fillna(0).values)
    host_tic = ps.hostname.astype(str).str.extract(r"^TIC[ -]?(\d+)", expand=False)
    by_name = {}
    for j, t in host_tic.dropna().items():
        by_name.setdefault(int(t), []).append(j)
    for i, (p, good, tic) in enumerate(zip(pos, ok.values, df.tic.values)):
        hits = set(tree.query_ball_point(p, chord)) if good else set()
        idx = sorted(hits | set(by_name.get(int(tic), [])))
        if not idx:
            continue
        rows = ps.loc[idx]
        df.iat[i, df.columns.get_loc("known_planet")] = "; ".join(rows.pl_name.astype(str))
        df.iat[i, df.columns.get_loc("known_planet_P")] = "; ".join(f"{v:.4g}" for v in rows.pl_orbper.fillna(np.nan))
        if good:
            df.iat[i, df.columns.get_loc("known_planet_sep")] = float(sep_arcsec(df[ra_col].values[i], df[dec_col].values[i],
                                                                                rows.ra.values, rows.dec.values).min())
    return df


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
        pos = xyz(df[ra_col].fillna(0).values, df[dec_col].fillna(0).values)
        d, i = tree.query(pos, k=1)
        sep = 2 * np.degrees(np.arcsin(np.clip(d / 2, 0, 1))) * 3600
        sep = np.where(df[ra_col].notna() & df[dec_col].notna(), sep, np.inf)
        df["toi_hf"] = np.where(sep < 10, h.toi_prefix.values[np.clip(i, 0, len(h) - 1)].astype(str), None)
        if "tic" in h:  # the complete snapshot also carries TIC IDs (column 'toi_id' in the source)
            by_tic = h.groupby("tic").toi_prefix.first().astype(str)
            df["toi_hf"] = df["toi_hf"].where(df["toi_hf"].notna(), df.tic.map(by_tic))
    else:
        df["toi_hf"] = None
    df["is_toi"] = df.toi_exofop.notna() | df.toi_hf.notna()
    return df


def spoc_tces(tic, sectors=None):
    """SPOC threshold-crossing events for a 2-min target, parsed from DV XML reports on the S3 bucket.

    Every sector directory holding the star's 2-min data (and multi-sector runs filed under their last
    sector) is listed; returns one dict per planetCandidate entry (period, epoch, depth, MES, EB flag).
    """
    import re
    import s3io
    from env import CATALOGS
    if sectors is None:
        m = pd.read_parquet(CATALOGS / "spoc_manifest.parquet", filters=[("tic", "==", int(tic))])
        sectors = sorted(set(m.sector.astype(int)))
        # multi-sector runs are filed under later sectors; also look at the standard multi-sector end sectors
        sectors = sorted(set(sectors) | {13, 26, 36, 39, 41, 55, 60, 65, 69, 78, 83, 92, 96})
    out = []
    for s in sectors:
        try:
            keys = [k for k, _ in s3io.list_keys(f"tess/public/tid/s{s:04d}/{s3io.tic_dir(tic)}/") if k.endswith("_dvr.xml")]
        except Exception:
            continue
        for k in keys:
            x = s3io.fetch(k).decode(errors="ignore")
            run = re.search(r"-(s\d{4}-s\d{4})-", k)
            for m in re.finditer(r"<dv:planetCandidate ([^>]*)>", x):
                a = dict(re.findall(r'(\w+)="([^"]*)"', m.group(1)))
                out.append({"run": run.group(1) if run else str(s), "planet": int(a.get("planetNumber", 0)),
                            "period": float(a.get("orbitalPeriodInDays", "nan")), "epoch_tjd": float(a.get("epochTjd", "nan")),
                            "mes": float(a.get("maxMultipleEventSigma", "nan")),
                            "n_tr": int(a.get("observedTransitCount", 0)),
                            "suspected_eb": a.get("suspectedEclipsingBinary") == "true"})
            # fitted depth per planet (all-transits fit)
            for m in re.finditer(r'planetNumber="(\d+)" transitModelName[^>]*>.*?name="transitDepthPpm" uncertainty="[^"]*" value="([^"]+)"', x, re.S):
                for o in out:
                    if o["run"] == (run.group(1) if run else str(s)) and o["planet"] == int(m.group(1)) and "depth_ppm" not in o:
                        o["depth_ppm"] = float(m.group(2))
    return out


def nearby_tois(ra, dec, radius_arcsec=150.0, period=None, ptol=0.01):
    """TOIs (ExoFOP Feb 2026 + NASA Oct 2026 snapshot) within radius; optionally only those whose period
    matches `period` or a small multiple/fraction of it (a signal often gets assigned to a neighbour)."""
    out = []
    e = exofop_toi()
    s = sep_arcsec(ra, dec, e.RA.values, e.Dec.values)
    for (_, r), sp in zip(e[s < radius_arcsec].iterrows(), s[s < radius_arcsec]):
        out.append({"source": "ExoFOP", "toi": r["TOI"], "tic": int(r["tic"]), "sep": float(sp), "period": r["period"],
                    "depth_ppm": r["depth_ppm"], "disp": f"{r['tfop']}/{r['tess_disp']}"})
    h = hf_toi()
    if h is not None:
        s = sep_arcsec(ra, dec, h.ra_deg.values, h.dec_deg.values)
        for (_, r), sp in zip(h[s < radius_arcsec].iterrows(), s[s < radius_arcsec]):
            out.append({"source": "NASA-Oct2026", "toi": r["toi_prefix"], "tic": None, "sep": float(sp),
                        "period": r.get("period_days"), "depth_ppm": r.get("transit_depth_ppm"), "disp": r.get("disposition")})
    if period is not None:
        def match(p):
            if p is None or not np.isfinite(p) or p <= 0:
                return False
            r_ = period / p
            return any(abs(r_ - k) < ptol * k for k in (1, 2, 3, 0.5, 1 / 3))
        out = [o for o in out if match(o["period"])]
    return out
