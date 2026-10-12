"""Behavioral screening checks on synthetic arrays, never protected science data."""
import numpy as np
import pytest

from adastra.screen import (
    ScreenConfig, clean_curve, compare_sectors, dense_period_grid,
    event_windows, flatten_curve, screen_curve,
)


def test_quality_filter_deduplicates_stably_without_clipping_or_mutating():
    time = np.array([2., 1., 1., 3., 4., 5., 6.])
    flux = np.array([100., 1., 2., .1, 1., np.nan, 1.])
    error = np.array([.1, .1, .1, .1, .1, .1, 0.])
    quality = np.array([0, 0, 0, 0, 1, 0, 0])
    original = flux.copy()
    t, f, e, details = clean_curve(time, flux, error, quality)
    np.testing.assert_array_equal(t, [1., 2., 3.])
    np.testing.assert_array_equal(f, [1., 100., .1])
    np.testing.assert_allclose(e, [.1, .1, .1])
    np.testing.assert_array_equal(flux, original)
    assert details["duplicate_times_removed"] == 1
    assert details["quality_finite_positive_samples"] == 4
    assert not np.shares_memory(t, time)


def test_weighted_screen_rejects_missing_errors_or_insufficient_usable_samples():
    with pytest.raises(ValueError, match="explicit flux errors"):
        clean_curve([1, 2], [1, 1], None, [0, 0])
    time = np.arange(100) * .02
    quality = np.zeros(100)
    quality[0] = 1
    with pytest.raises(ValueError, match="Only 99 samples"):
        screen_curve(time, np.ones(100), np.ones(100), quality)


def test_running_median_keeps_separate_gap_baselines_and_scales_errors():
    time = np.r_[np.arange(7)*.1, 10+np.arange(7)*.1]
    flux = np.r_[np.ones(7), 5*np.ones(7)]
    error = np.full(14, .1)
    flattened, scaled_error, trend, segments = flatten_curve(time, flux, error, ScreenConfig())
    np.testing.assert_allclose(flattened, 1)
    np.testing.assert_allclose(trend, flux)
    np.testing.assert_allclose(scaled_error, error/flux)
    assert [row["samples"] for row in segments] == [7, 7]


def test_period_grid_resolves_declared_phase_drift_and_refuses_coarsening():
    config = ScreenConfig()
    periods, diagnostics = dense_period_grid(27., config)
    assert np.all(np.diff(periods) > 0)
    assert periods[0] == pytest.approx(.5)
    assert periods[-1] == pytest.approx(13)
    assert diagnostics["nearest_grid_worst_case_timing_drift_days"] <= min(config.durations_days)/4+1e-14
    frequencies = 1/periods
    np.testing.assert_allclose(np.diff(frequencies), np.diff(frequencies)[0], rtol=1e-10)
    with pytest.raises(ValueError, match="refusing coarsening"):
        dense_period_grid(27., ScreenConfig(maximum_grid_points=100))
    with pytest.raises(ValueError, match="Baseline cannot cover"):
        dense_period_grid(.9, config)


def test_window_coverage_is_distinct_from_positive_depth_or_significance():
    # Three sampled predicted events; only the first has enough local coverage.
    time = np.array([-.14, -.12, -.10, -.02, .02, .10, .12, .14, 1., 1.98, 2.02])
    flux = np.ones(len(time))
    flux[[3, 4]] -= .01
    result = event_windows(time, flux, period=1., epoch=0., duration=.1)
    assert result["sampled_windows"] == 3
    assert result["positive_local_median_depth_windows"] == 1
    assert result["events"][0]["local_median_depth_ppm"] == pytest.approx(10000.)
    assert result["events"][1]["local_median_depth_ppm"] is None
    assert result["events"][2]["local_median_depth_ppm"] is None


def test_synthetic_transit_recovers_period_and_phase_after_quality_filter():
    rng = np.random.default_rng(483)
    time = np.arange(0, 12, .02)
    period, epoch, duration = 2., .41, .08
    phase = (time-epoch+period/2) % period-period/2
    flux = 1+rng.normal(0, .0001, len(time))
    flux[np.abs(phase) < duration/2] -= .01
    quality = np.zeros(len(time), dtype=int)
    # A severe flagged feature must not influence period selection.
    quality[100:105] = 1
    flux[100:105] = .1
    result = screen_curve(time, flux, np.full(len(time), .0001), quality,
                          ScreenConfig(maximum_period_days=5))
    peak = result["bls_peaks"][0]
    assert abs(peak["period_days"]/period-1) < .005
    phase_error = abs((peak["epoch_btjd"]-epoch+period/2) % period-period/2)
    assert phase_error < duration/2
    assert peak["passes_event_window_triage"]
    assert peak["original_events"]["sampled_windows"] == 6
    assert peak["depth_ppm"] == pytest.approx(10000, rel=.08)
    assert result["cleaning"]["analysis_samples"] == len(time)-5


def test_synthetic_sinusoid_recovers_descriptive_period_and_amplitude():
    time = np.arange(0, 12, .02)
    period, amplitude = .71, .006
    flux = 1+amplitude*np.sin(2*np.pi*time/period)
    result = screen_curve(time, flux, np.full(len(time), .0001), np.zeros(len(time)),
                          ScreenConfig(maximum_period_days=5, maximum_peaks=1))
    ls = result["lomb_scargle"]
    assert abs(ls["period_days"]/period-1) < .01
    assert ls["amplitude_ppm"] == pytest.approx(amplitude*1e6, rel=.03)
    assert "no survey false-alarm probability" in ls["interpretation"]


def test_constant_flux_has_no_positive_depth_triage_or_ls_period():
    time = np.arange(0, 6, .02)
    result = screen_curve(time, np.ones(len(time)), np.full(len(time), .001), np.zeros(len(time)))
    assert result["lomb_scargle"]["status"] == "constant_flux"
    assert result["lomb_scargle"]["period_days"] is None
    assert all(not peak["passes_event_window_triage"] for peak in result["bls_peaks"])


def test_sector_comparison_checks_all_pairs_and_requires_event_triage():
    def curve(periods, accepted):
        return {"bls_peaks": [{"rank": i+1, "period_days": period,
                               "passes_event_window_triage": i in accepted}
                              for i, period in enumerate(periods)],
                "lomb_scargle": {"period_days": None}}
    result = compare_sectors(curve([1., 2., 4.], {0, 1}), curve([6., 2.01, 4.01], {0, 1, 2}))
    assert len(result["pairs"]) == 9
    direct = [pair for pair in result["pairs"] if pair["direct_inspection_priority"]]
    assert [(p["first_rank"], p["second_rank"]) for p in direct] == [(2, 2)]
    assert result["any_harmonic_inspection_priority"]
    assert not result["ls_period_agreement"]
