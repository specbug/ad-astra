"""Exploratory QLP screening, reconstructed after loss of the old workspace.

Scores and event-window counts are descriptive, not detection probabilities.
QLP KSPSAP is already processed by its authors; ``unflattened`` below means no
additional running median in this implementation, not untouched detector flux.
"""
from dataclasses import asdict, dataclass
import math

import numpy as np
from astropy.timeseries import BoxLeastSquares, LombScargle
from scipy.ndimage import median_filter


@dataclass(frozen=True)
class ScreenConfig:
    minimum_samples: int = 100
    gap_days: float = 0.5
    median_window_days: float = 0.5
    minimum_period_days: float = 0.5
    maximum_period_days: float = 13.0
    durations_days: tuple = (0.04, 0.065, 0.1, 0.15, 0.2)
    maximum_grid_points: int = 200000
    bls_oversample: int = 10
    maximum_peaks: int = 3
    masking_duration_multiplier: float = 1.5
    maximum_masked_fraction: float = 0.5
    ls_minimum_period_days: float = 0.05
    ls_maximum_period_days: float = 20.0
    ls_samples_per_peak: int = 5
    period_agreement_fraction: float = 0.01
    minimum_sampled_windows: int = 3
    minimum_positive_windows: int = 3

    def __post_init__(self):
        if self.minimum_samples < 100 or self.maximum_grid_points < 2:
            raise ValueError("At least 100 samples and two grid points are required")
        if not (0 < self.gap_days and 0 < self.median_window_days):
            raise ValueError("Gap and median window must be positive")
        if not 0 < self.minimum_period_days < self.maximum_period_days:
            raise ValueError("Invalid BLS period range")
        if not self.durations_days or any(not 0 < x < self.minimum_period_days for x in self.durations_days):
            raise ValueError("BLS durations must be positive and below the shortest period")
        if not 1 <= self.maximum_peaks <= 3:
            raise ValueError("Only one to three deterministic peaks are supported")
        if not (self.masking_duration_multiplier >= 1 and 0 < self.maximum_masked_fraction <= 0.5):
            raise ValueError("Invalid masking configuration")
        if not 0 < self.ls_minimum_period_days < self.ls_maximum_period_days:
            raise ValueError("Invalid LS period range")
        if self.bls_oversample < 1 or self.ls_samples_per_peak < 1:
            raise ValueError("Periodogram oversampling must be positive")
        if not 0 < self.period_agreement_fraction < 0.1:
            raise ValueError("Period agreement fraction is invalid")
        if self.minimum_sampled_windows < 3 or self.minimum_positive_windows < 3:
            raise ValueError("At least three sampled and positive windows are required")


def clean_curve(time, flux, error, quality):
    """Quality-zero, finite positive values; stable sort and first-time dedup.

    No sigma clipping occurs. Errors are mandatory for this weighted screen.
    Arrays returned are independent copies, with flux/error median normalized.
    """
    if error is None:
        raise ValueError("This weighted screen requires explicit flux errors")
    arrays = [np.asarray(x) for x in (time, flux, error, quality)]
    if any(x.ndim != 1 for x in arrays) or len({len(x) for x in arrays}) != 1:
        raise ValueError("Expected equally sized one-dimensional arrays")
    t, f, e = (np.asarray(x, dtype=float) for x in arrays[:3])
    q = arrays[3]
    good = (q == 0) & np.isfinite(t) & np.isfinite(f) & (f > 0) & np.isfinite(e) & (e > 0)
    t, f, e = t[good], f[good], e[good]
    order = np.argsort(t, kind="stable")
    t, f, e = t[order], f[order], e[order]
    unique = np.r_[True, np.diff(t) > 0] if len(t) else np.empty(0, dtype=bool)
    duplicates = int(np.count_nonzero(~unique))
    t, f, e = t[unique], f[unique], e[unique]
    if len(t) == 0:
        raise ValueError("No usable samples")
    normalizer = float(np.median(f))
    return t, f / normalizer, e / normalizer, {
        "input_samples": len(arrays[0]), "quality_finite_positive_samples": int(good.sum()),
        "duplicate_times_removed": duplicates, "analysis_samples": len(t),
        "median_flux_normalization": normalizer, "baseline_days": float(np.ptp(t)),
        "filter": "quality == 0; finite time, positive finite flux/error; no clipping",
        "duplicate_policy": "first qualifying input row at each exact time, after stable sort",
    }


def flatten_curve(time, flux, error, config):
    """Divide by an odd-width running median independently across >gap gaps."""
    trend = np.empty_like(flux)
    bounds = np.r_[0, np.flatnonzero(np.diff(time) > config.gap_days) + 1, len(time)]
    segments = []
    for start, end in zip(bounds[:-1], bounds[1:]):
        cadence = float(np.median(np.diff(time[start:end]))) if end - start > 1 else config.median_window_days
        size = max(1, int(round(config.median_window_days / cadence)))
        size += size % 2 == 0
        trend[start:end] = median_filter(flux[start:end], size=size, mode="nearest")
        segments.append({"start_index": int(start), "end_index_exclusive": int(end),
                         "samples": int(end-start), "median_cadence_days": cadence,
                         "median_filter_samples": size})
    if not np.all(np.isfinite(trend) & (trend > 0)):
        raise ValueError("Nonpositive or nonfinite running-median trend")
    return flux / trend, error / trend, trend, segments


def dense_period_grid(baseline, config):
    upper = min(config.maximum_period_days, baseline / 2)
    if upper <= config.minimum_period_days:
        raise ValueError("Baseline cannot cover two allowed periods")
    low_frequency, high_frequency = 1 / upper, 1 / config.minimum_period_days
    # Half a frequency step at the longest period, over the full baseline,
    # accumulates no more than a quarter of the shortest trial duration.
    n = math.ceil((high_frequency-low_frequency) * 2 * baseline * upper / min(config.durations_days)) + 1
    if n > config.maximum_grid_points:
        raise ValueError(f"Required {n} grid points exceeds cap {config.maximum_grid_points}; refusing coarsening")
    frequency = np.linspace(low_frequency, high_frequency, n)
    periods = (1 / frequency)[::-1]
    df = float(frequency[1]-frequency[0])
    return periods, {
        "points": n, "minimum_period_days": float(periods[0]), "maximum_period_days": float(periods[-1]),
        "uniform_frequency_step_per_day": df,
        "nearest_grid_worst_case_timing_drift_days": 0.5 * df * baseline * upper,
        "required_maximum_drift_days": min(config.durations_days) / 4,
    }


def event_windows(time, flux, period, epoch, duration):
    """Count sampled predicted windows; positive local median depth is not S/N.

    In-window samples lie within duration/2; local comparison samples lie between
    duration/2 and min(1.5*duration, period/2). One sample suffices for coverage.
    A depth needs >=2 in-window and >=3 comparison samples. No statistical
    significance or number of detected transits can be inferred from counts.
    """
    first = math.ceil((time[0]-duration/2-epoch) / period)
    last = math.floor((time[-1]+duration/2-epoch) / period)
    events = []
    for index in range(first, last + 1):
        center = epoch + index * period
        distance = np.abs(time-center)
        inside = distance <= duration/2
        outside = (distance > duration/2) & (distance <= min(1.5*duration, period/2))
        ni, no = int(inside.sum()), int(outside.sum())
        if not ni:
            continue
        depth = float(np.median(flux[outside]) - np.median(flux[inside])) if ni >= 2 and no >= 3 else None
        events.append({"index": index, "center_btjd": float(center), "in_window_samples": ni,
                       "local_comparison_samples": no, "local_median_depth_ppm": None if depth is None else depth*1e6})
    return {"sampled_windows": len(events),
            "positive_local_median_depth_windows": sum(x["local_median_depth_ppm"] is not None and x["local_median_depth_ppm"] > 0 for x in events),
            "events": events,
            "interpretation": "Coverage and signs of local sample medians, NOT statistically detected transit events"}


def _finite(value):
    value = float(value)
    return value if np.isfinite(value) else None


def screen_curve(time, flux, error, quality, config=None):
    """Return up to three masked BLS peaks and a descriptive weighted LS peak."""
    config = config or ScreenConfig()
    if isinstance(config, dict):
        config = ScreenConfig(**config)
    t, f, e, clean = clean_curve(time, flux, error, quality)
    if len(t) < config.minimum_samples:
        raise ValueError(f"Only {len(t)} samples; need {config.minimum_samples}")
    y, dy, trend, segments = flatten_curve(t, f, e, config)
    periods, grid = dense_period_grid(float(np.ptp(t)), config)
    selected = np.ones(len(t), dtype=bool)
    peaks = []
    stop_reason = "maximum_peaks_reached"
    for rank in range(1, config.maximum_peaks+1):
        if np.count_nonzero(selected) < config.minimum_samples:
            stop_reason = "fewer_than_minimum_selected_samples"
            break
        if 1 - selected.mean() > config.maximum_masked_fraction:
            stop_reason = "cumulative_mask_fraction_exceeded"
            break
        bls = BoxLeastSquares(t[selected], y[selected], dy[selected])
        power = bls.power(periods, config.durations_days, objective="likelihood", method="fast", oversample=config.bls_oversample)
        if not np.any(np.isfinite(power.power)):
            stop_reason = "no_finite_bls_power"
            break
        k = int(np.nanargmax(power.power))
        period, epoch, duration = (float(x[k]) for x in (power.period, power.transit_time, power.duration))
        original_events = event_windows(t, y, period, epoch, duration)
        selection_events = event_windows(t[selected], y[selected], period, epoch, duration)
        # The exact upper endpoint lies within half of the final grid spacing.
        boundary = period >= periods[-1] - (periods[-1]-periods[-2])/2
        peak = {"rank": rank, "period_days": period, "epoch_btjd": epoch, "duration_days": duration,
                "depth_ppm": _finite(power.depth[k]*1e6), "depth_error_ppm": _finite(power.depth_err[k]*1e6),
                "depth_snr_conditional": _finite(power.depth_snr[k]), "bls_power": _finite(power.power[k]),
                "selected_samples": int(selected.sum()), "masked_fraction_at_search": float(1-selected.mean()),
                "at_upper_period_boundary": bool(boundary), "original_events": original_events,
                "selection_events": selection_events}
        peak["passes_event_window_triage"] = bool(not boundary and all(
            events["sampled_windows"] >= config.minimum_sampled_windows and
            events["positive_local_median_depth_windows"] >= config.minimum_positive_windows
            for events in (original_events, selection_events)))
        peaks.append(peak)
        phase = (t-epoch+period/2) % period - period/2
        selected &= np.abs(phase) > config.masking_duration_multiplier * duration / 2
    upper_ls = min(config.ls_maximum_period_days, np.ptp(t)/2)
    ls_result = {"input": "median-normalized KSPSAP flux without additional running-median division; upstream QLP processing remains",
                 "interpretation": "Descriptive strongest weighted LS sinusoid; no survey false-alarm probability"}
    if np.ptp(f) == 0:
        ls_result.update(status="constant_flux", period_days=None, power=None)
    elif upper_ls <= config.ls_minimum_period_days:
        ls_result.update(status="insufficient_baseline", period_days=None, power=None)
    else:
        ls = LombScargle(t, f, e, fit_mean=True, center_data=True)
        frequencies = ls.autofrequency(minimum_frequency=1/upper_ls, maximum_frequency=1/config.ls_minimum_period_days,
                                       samples_per_peak=config.ls_samples_per_peak)
        powers = ls.power(frequencies, normalization="standard")
        if not np.any(np.isfinite(powers)):
            ls_result.update(status="no_finite_power", period_days=None, power=None)
        else:
            k = int(np.nanargmax(powers))
            params = ls.model_parameters(frequencies[k])
            ls_result.update(status="complete", period_days=float(1/frequencies[k]), power=_finite(powers[k]),
                             amplitude_ppm=float(np.hypot(params[1], params[2])*1e6), trial_frequencies=len(frequencies))
    return {"status": "complete", "config": asdict(config), "cleaning": clean, "segments": segments,
            "grid": grid, "bls_peaks": peaks, "bls_stop_reason": stop_reason, "lomb_scargle": ls_result,
            "limitations": ["Exploratory ranking only; multiplicity, correlated noise and selection effects are uncalibrated",
                            "Positive local median-depth windows are not transit detections",
                            "Masked peak selection is adaptive within each curve",
                            "Median edge handling and pre-existing QLP corrections can alter astrophysical shapes"]}


def compare_sectors(first, second, config=None):
    """Evaluate every up-to-nine peak pair, including harmonic coincidences."""
    config = config or ScreenConfig()
    if isinstance(config, dict):
        config = ScreenConfig(**config)
    pairs = []
    for a in first["bls_peaks"]:
        for b in second["bls_peaks"]:
            p1, p2 = a["period_days"], b["period_days"]
            relative = abs(p1-p2) / min(p1, p2)
            ratio = max(p1, p2) / min(p1, p2)
            harmonic = min((2, 3), key=lambda h: abs(ratio/h-1))
            harmonic_delta = abs(ratio/harmonic-1)
            qualified = a["passes_event_window_triage"] and b["passes_event_window_triage"]
            pairs.append({"first_rank": a["rank"], "second_rank": b["rank"], "first_period_days": p1,
                          "second_period_days": p2, "fractional_direct_difference": relative,
                          "nearest_harmonic": harmonic, "fractional_harmonic_difference": harmonic_delta,
                          "both_pass_event_window_triage": bool(qualified),
                          "direct_inspection_priority": bool(qualified and relative <= config.period_agreement_fraction),
                          "harmonic_inspection_priority": bool(qualified and harmonic_delta <= config.period_agreement_fraction)})
    la, lb = first["lomb_scargle"].get("period_days"), second["lomb_scargle"].get("period_days")
    ls_agreement = bool(la is not None and lb is not None and abs(la-lb)/min(la, lb) <= config.period_agreement_fraction)
    return {"pairs": pairs, "any_direct_inspection_priority": any(p["direct_inspection_priority"] for p in pairs),
            "any_harmonic_inspection_priority": any(p["harmonic_inspection_priority"] for p in pairs),
            "ls_period_agreement": ls_agreement, "interpretation": "Inspection priorities, not validated candidate classifications"}
