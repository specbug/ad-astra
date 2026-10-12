"""Science-facing diagnostics checks using synthetic arrays, never archive flux."""
import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest

SPEC = importlib.util.spec_from_file_location("qlp_diagnostics", Path(__file__).parents[1] / "scripts" / "qlp_diagnostics.py")
diagnostics = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(diagnostics)


def test_close_periods_can_fail_independent_phase_coherence():
    first = {"period_days": 2.0, "epoch_btjd": 0.25, "duration_days": 0.2}
    second = {**first, "period_days": 2.012}
    result = diagnostics.pair_phase(first, second, np.arange(30, 55, 0.02))
    assert abs(first["period_days"] - second["period_days"]) / first["period_days"] < 0.01
    assert result["median_absolute_center_offset_days"] > result["window_overlap_tolerance_days"]
    assert not result["median_centers_overlap"]


def test_integer_epoch_cycle_shift_does_not_create_phase_mismatch():
    first = {"period_days": 2.0, "epoch_btjd": 0.25, "duration_days": 0.2}
    result = diagnostics.pair_phase(first, {**first, "epoch_btjd": 100.25}, np.arange(30, 55, 0.02))
    assert result["median_centers_overlap"]
    assert result["maximum_absolute_center_offset_days"] == 0


def test_odd_even_depths_and_single_event_domination_are_visible():
    time = np.arange(0.0, 10.0, 0.01)
    flux = np.ones_like(time)
    for cycle in range(5):
        depth = 0.006 if cycle % 2 else 0.003
        flux[np.abs(time - (cycle * 2.0 + 0.25)) <= 0.1] -= depth
    events = diagnostics.event_windows(time, flux, 2.0, 0.25, 0.2)
    summary = diagnostics.depth_robustness(events)
    assert summary["measurable_windows"] == summary["positive_windows"] == 5
    assert summary["odd_minus_even_median_depth_ppm"] == pytest.approx(3000)
    one_deep = {"sampled_windows": 4, "events": [
        {"index": index, "local_median_depth_ppm": depth}
        for index, depth in enumerate([10000.0, 0.0, -100.0, 0.0])
    ]}
    summary = diagnostics.depth_robustness(one_deep)
    assert summary["positive_windows"] == 1
    assert summary["largest_positive_depth_fraction_of_positive_depth_sum"] == 1
    assert summary["leave_one_event_out_mean_depth_range_ppm"][0] < 0


def test_all_recorded_sectors_checked_before_any_cached_flux_open(tmp_path, monkeypatch):
    root = Path(__file__).parents[1]
    payload = {
        "status": "complete", "original_sources_unchanged": True,
        "all_requested_targets_accounted_for": True, "target_count": 1,
        "declaration": {"allowed_discovery_sectors": [1, 2], "protected_holdouts_accessed": False,
                        "source_sha256": {name: diagnostics.sha256_file(root / name) for name in
                                          ("src/adastra/archive.py", "src/adastra/qlp.py", "src/adastra/screen.py")}},
        "results": [{"tic": 999, "status": "complete", "sectors": [
            {"sector": 1, "status": "complete"}, {"sector": 27, "status": "complete"}]}],
    }
    path = tmp_path / "results.json"
    path.write_text(json.dumps(payload))
    def forbidden(*args, **kwargs):
        pytest.fail("No cached flux may be accessed when any result lists a protected sector")
    monkeypatch.setattr(diagnostics, "load_cached_sector", forbidden)
    with pytest.raises(ValueError, match="outside literal S1/S2"):
        diagnostics.run(path, tmp_path / "output", root)
    assert not (tmp_path / "output").exists()


def test_uri_sector_mismatch_rejected_before_any_file_read(tmp_path):
    row = {"sector": 1, "status": "complete", "input_path": "missing.fits",
           "uri": "mast:HLSP/qlp/s0027/0000/0000/0000/0999/hlsp_qlp_tess_ffi_s0027-0000000000000999_tess_v01_llc.fits"}
    with pytest.raises(ValueError, match="exact S1/S2 product"):
        diagnostics.load_cached_sector(row, 999, tmp_path)
