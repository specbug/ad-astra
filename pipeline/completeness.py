"""Blind-recovery test of the survey on known TOIs (Track B completeness).

For every TOI with a period > 15 d (at most ~2 transits per sector), host in the survey sample, and a
transit predicted inside a surveyed sector, check whether the survey catalogued a dip within 0.3 d of
the predicted time (any Track B event), and whether it survived stage-2 as planet-like.
"""
import glob
import json

import env  # noqa: F401
import numpy as np
import pandas as pd

import known
import s3io
from env import RESULTS, CATALOGS

BTJD = 2457000.0


def sector_spans(sectors):
    """Time span (BTJD) of each sector from one QLP light curve header (cached)."""
    path = CATALOGS / "sector_spans.json"
    spans = json.loads(path.read_text()) if path.exists() else {}
    for s in sectors:
        if str(s) in spans:
            continue
        tl = s3io.read_target_list("qlp", s)
        for tic in tl.tic.astype("int64").values[:50]:
            lc = s3io.load("qlp", int(tic), s)
            if lc is not None and lc["time"].size > 100:
                spans[str(s)] = [float(lc["time"].min()), float(lc["time"].max())]
                break
    path.write_text(json.dumps(spans))
    return {int(k): v for k, v in spans.items()}


def run():
    toi = known.exofop_toi()
    toi = toi[(toi.period > 15) & (toi.depth_ppm > 1500)].copy()
    sample = pd.read_parquet(RESULTS / "samples" / "survey_sample.parquet").tic
    toi = toi[toi.tic.isin(sample)]
    qlp = pd.read_parquet(CATALOGS / "qlp_manifest.parquet", filters=[("tic", "in", list(toi.tic.unique()))])
    done = sorted(int(p.split("/s")[-1][:4]) for p in glob.glob(str(RESULTS / "survey" / "s[0-9][0-9][0-9][0-9]"))
                  if (RESULTS / "survey" / p.split("/")[-1] / "summary.json").exists())
    spans = sector_spans(done)
    ev = []
    for s in done:
        f = RESULTS / "survey" / f"s{s:04d}" / "trackb_events.parquet"
        if f.exists():
            d = pd.read_parquet(f)
            ev.append(d[~d.invert.astype(bool)])
    ev = pd.concat(ev, ignore_index=True) if ev else pd.DataFrame(columns=["tic", "sector", "t0"])
    st2 = RESULTS / "candidates" / "trackb_stage2.parquet"
    st2 = pd.read_parquet(st2) if st2.exists() else None
    rows = []
    for _, r in toi.iterrows():
        secs = set(qlp[qlp.tic == r.tic].sector.astype(int)) & set(done)
        for s in secs:
            a, b = spans[s]
            ep = r.epoch_bjd - BTJD
            k = np.arange(np.ceil((a - ep) / r.period), np.floor((b - ep) / r.period) + 1)
            for tc in ep + k * r.period:
                e = ev[(ev.tic == r.tic) & (ev.sector == s) & (np.abs(ev.t0 - tc) < 0.3)]
                pl = None
                if st2 is not None and len(e):
                    m = st2[(st2.tic == r.tic) & (st2.sector == s) & (np.abs(st2.t0 - tc) < 0.3) & ~st2.invert.astype(bool)]
                    pl = bool(m.planet_like.any()) if len(m) else False
                rows.append({"toi": r.TOI, "tic": int(r.tic), "sector": s, "t_pred": tc, "period": r.period,
                             "depth_ppm": r.depth_ppm, "tmag": r["TESS Mag"], "found": len(e) > 0, "planet_like": pl,
                             "disp": f"{r.tfop}/{r.tess_disp}"})
    df = pd.DataFrame(rows)
    out = RESULTS / "candidates"
    out.mkdir(parents=True, exist_ok=True)
    df.to_csv(out / "completeness_toi.csv", index=False)
    if len(df):
        print(f"predicted TOI transits in surveyed sectors: {len(df)}; found as Track B events: {df.found.mean():.2%}")
        for lo, hi in ((1500, 3000), (3000, 6000), (6000, 12000), (12000, 1e6)):
            x = df[(df.depth_ppm >= lo) & (df.depth_ppm < hi)]
            if len(x):
                print(f"  depth {lo:>6.0f}-{hi:<7.0f} ppm: {len(x):4d} transits, found {x.found.mean():.2%}, "
                      f"planet-like among found {x[x.found].planet_like.mean() if x.found.any() else float('nan'):.2%}")
    return df


if __name__ == "__main__":
    run()
