#!/usr/bin/env python3
"""Reconstruct legacy QLP ingestion with only the known pi Men c control."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import tempfile

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "adastra-matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(Path(tempfile.gettempdir()) / "adastra-cache"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy.timeseries import BoxLeastSquares
from scipy.ndimage import median_filter

from adastra.archive import ArchiveClient, query_mast
from adastra.qlp import read_qlp


TIC = 261136679
SECTORS = (1, 27)
KNOWN_PERIOD = 6.2678399
KNOWN_EPOCH_BTJD = 1425.789204
KNOWN_DURATION = 2.952024 / 24


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def flatten(time, flux, width=.5):
    trend = np.empty_like(flux)
    cuts = np.r_[0, np.flatnonzero(np.diff(time) > .5) + 1, len(time)]
    for start, end in zip(cuts[:-1], cuts[1:]):
        cadence = np.median(np.diff(time[start:end])) if end-start > 1 else width
        size = max(1, int(round(width / cadence)))
        size += (size % 2 == 0)
        trend[start:end] = median_filter(flux[start:end], size=size, mode="nearest")
    return flux/trend, trend


def fixed_ephemeris_depths(time, flux):
    """Local median differences at the published phase; no phase optimization."""
    first = int(np.ceil((time.min()-KNOWN_EPOCH_BTJD)/KNOWN_PERIOD))
    last = int(np.floor((time.max()-KNOWN_EPOCH_BTJD)/KNOWN_PERIOD))
    events = []
    for cycle in range(first, last+1):
        center = KNOWN_EPOCH_BTJD + cycle*KNOWN_PERIOD
        distance = np.abs(time-center)
        inside = distance < KNOWN_DURATION/2
        outside = (distance > KNOWN_DURATION) & (distance < .4)
        measured = int(inside.sum()) >= 2 and int(outside.sum()) >= 4
        events.append({"cycle": cycle, "center_btjd": center,
                       "in_transit_samples": int(inside.sum()), "local_out_samples": int(outside.sum()),
                       "median_depth_ppm": float((np.median(flux[outside])-np.median(flux[inside]))*1e6)
                       if measured else None})
    return {"period_days": KNOWN_PERIOD, "epoch_btjd": KNOWN_EPOCH_BTJD,
            "duration_days": KNOWN_DURATION,
            "source": "NASA Exoplanet Archive ps default row: Kunovac Hodzic et al. 2021",
            "depth_definition": "local out-of-transit median minus in-transit median on flattened flux",
            "events": events, "interpretation": "descriptive depths; no estimated detection significance"}


def run(out, cache):
    out.mkdir(parents=True, exist_ok=True)
    declaration = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "purpose": "fresh reconstruction: legacy QLP ingestion and known-planet control",
        "tic_id": TIC, "sectors": list(SECTORS),
        "known_period_days": KNOWN_PERIOD,
        "known_epoch_btjd": KNOWN_EPOCH_BTJD,
        "known_duration_days": KNOWN_DURATION,
        "known_period_source": "NASA Exoplanet Archive pscomppars, pi Men c / TIC 261136679, retrieved 2026-10-10 UTC",
        "flux": "KSPSAP_FLUX", "quality": "QUALITY == 0; finite positive flux and errors",
        "detrending_window_days": .5, "gap_threshold_days": .5,
        "search_period_days": [.5, 13.], "period_limit": "also no longer than half measured baseline",
        "durations_days": [.04, .065, .1, .15, .2],
        "frequency_grid": "nearest-grid accumulated timing drift <= minimum duration / 4",
        "control_recovery_criterion": f"strongest BLS period within 1% of {KNOWN_PERIOD} days",
        "source_sha256": digest(__file__),
        "limitations": ["Known control, not an independent discovery or held-out test",
                        "Conditional BLS depth S/N is not a survey false-alarm probability",
                        "No discovery target or holdout flux is requested"],
    }
    (out/"declared-control.json").write_text(json.dumps(declaration, indent=2)+"\n")
    request = {"service": "Mast.Caom.Filtered", "params": {
        "columns": "*", "filters": [
            {"paramName": "target_name", "values": [str(TIC)]},
            {"paramName": "provenance_name", "values": ["QLP"]}]},
        "format": "json", "pagesize": 1000, "page": 1}
    response = query_mast(request)
    (out/"mast-observations.json").write_text(json.dumps(response, indent=2)+"\n")
    if response.get("status") != "COMPLETE":
        raise RuntimeError("MAST query did not complete")
    rows = response.get("data", [])
    if len(rows) >= 1000:
        raise RuntimeError("Unexpectedly paginated control response; refusing incomplete selection")
    client = ArchiveClient(cache)
    results = []
    for sector in SECTORS:
        selected = [r for r in rows if int(r.get("sequence_number", -1)) == sector
                    and str(r.get("dataURL", "")).endswith("_llc.fits")]
        if len(selected) != 1:
            raise RuntimeError(f"Expected one QLP sector-{sector} control; got {len(selected)}")
        path = client.download(selected[0]["dataURL"])
        lc = read_qlp(path)
        sap = read_qlp(path, "SAP_FLUX")
        if lc.metadata["tic_id"] != TIC or lc.metadata["sector"] != sector:
            raise ValueError("Downloaded identity differs from the control declaration")
        use = lc.good_mask()
        t, flux, errors = lc.time_btjd[use], lc.flux[use], lc.flux_error[use]
        order = np.argsort(t, kind="stable")
        t, flux, errors = t[order], flux[order], errors[order]
        unique = np.r_[True, np.diff(t) > 0]
        duplicates = int(np.sum(~unique))
        t, flux, errors = t[unique], flux[unique], errors[unique]
        norm = np.median(flux)
        f, trend = flatten(t, flux/norm)
        e = errors/norm/trend
        baseline = float(np.ptp(t))
        max_period = min(13., baseline/2)
        fmin, fmax = 1/max_period, 1/.5
        n = int(np.ceil((fmax-fmin)*2*baseline*max_period/.04)) + 1
        periods = (1/np.linspace(fmin, fmax, n))[::-1]
        bls = BoxLeastSquares(t, f, e)
        power = bls.power(periods, [.04, .065, .1, .15, .2])
        k = int(np.argmax(power.power))
        fixed = bls.power([KNOWN_PERIOD], [.04, .065, .1, .15, .2])
        recovered = bool(abs(float(power.period[k])/KNOWN_PERIOD-1) < .01)
        row = {
            "sector": sector, "tic_id": TIC, "uri": selected[0]["dataURL"],
            "local_path": str(path), "sha256": digest(path),
            "ingestion": lc.metadata, "sap_error_available": sap.flux_error is not None,
            "quality_zero_finite_positive_samples": int(np.sum(use)),
            "duplicate_times_removed": duplicates, "analysis_samples": len(t),
            "baseline_days": baseline, "median_cadence_seconds": float(np.median(np.diff(t))*86400),
            "grid_period_count": n, "strongest_period_days": float(power.period[k]),
            "strongest_duration_days": float(power.duration[k]),
            "strongest_epoch_btjd": float(power.transit_time[k]),
            "strongest_depth_ppm": float(power.depth[k]*1e6),
            "strongest_depth_snr_conditional": float(power.depth_snr[k]),
            "known_period_recovered_as_strongest": recovered,
            "fixed_published_ephemeris": fixed_ephemeris_depths(t, f),
            "fixed_known_period_diagnostic": {
                "period_days": KNOWN_PERIOD,
                "fitted_epoch_btjd": float(fixed.transit_time[0]),
                "fitted_duration_days": float(fixed.duration[0]),
                "depth_ppm": float(fixed.depth[0]*1e6),
                "depth_snr_conditional": float(fixed.depth_snr[0]),
                "interpretation": "period fixed; duration and phase fitted to this same curve"},
        }
        results.append(row)
        fig, axes = plt.subplots(3, 1, figsize=(10, 9), constrained_layout=True)
        axes[0].plot(t, (flux/norm-1)*1e6, ".", ms=2, alpha=.5)
        axes[0].set(xlabel="BTJD (BJD TDB − 2457000)", ylabel="QLP flux (ppm)")
        axes[1].plot(periods, power.power, lw=.8)
        axes[1].axvline(KNOWN_PERIOD, color="tab:red", ls="--", label="Known pi Men c period")
        axes[1].set(xlabel="Trial period (days)", ylabel="BLS power")
        axes[1].legend()
        phase = ((t-float(fixed.transit_time[0])+KNOWN_PERIOD/2)%KNOWN_PERIOD-KNOWN_PERIOD/2)
        axes[2].plot(phase, (f-1)*1e6, ".", ms=2, alpha=.35)
        edges = np.linspace(-.4, .4, 41)
        ys = [np.median((f[(phase>=a)&(phase<b)]-1)*1e6)
              if np.any((phase>=a)&(phase<b)) else np.nan for a,b in zip(edges[:-1],edges[1:])]
        axes[2].plot((edges[:-1]+edges[1:])/2, ys, "o-", ms=3)
        axes[2].set(xlim=(-.4,.4), xlabel="Days from fitted transit epoch (known period)", ylabel="Detrended flux (ppm)")
        fig.suptitle(f"Known control: pi Men c, QLP sector {sector}; strongest period {power.period[k]:.6f} d")
        fig.savefig(out/f"sector-{sector}.png", dpi=150)
        fig.savefig(out/f"sector-{sector}.pdf")
        plt.close(fig)
    payload = {"declaration": declaration, "completed_utc": datetime.now(timezone.utc).isoformat(),
               "results": results, "scope": "Only these two known-control light curves were downloaded"}
    (out/"results.json").write_text(json.dumps(payload, indent=2)+"\n")
    print(json.dumps({"output": str(out), "results": [{k:r[k] for k in
        ["sector", "analysis_samples", "strongest_period_days", "known_period_recovered_as_strongest"]}
        for r in results]}, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path("runs/reconstruction/qlp-control"))
    parser.add_argument("--cache", type=Path, default=Path("data/cache"))
    args = parser.parse_args()
    run(args.output, args.cache)
