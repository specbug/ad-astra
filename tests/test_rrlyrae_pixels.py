"""Meaningful arithmetic and holdout guards for the contamination diagnostic."""
import importlib.util
from pathlib import Path

import numpy as np
import pytest

spec = importlib.util.spec_from_file_location(
    "rrlyrae_pixels", Path(__file__).resolve().parents[1] / "scripts/qlp_rrlyrae_pixels.py")
rr = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rr)


def test_fixed_harmonics_recovered_with_segment_baselines_and_background():
    time = np.r_[np.linspace(0, 12, 600), np.linspace(14, 26, 600)]
    background = np.cos(time / 4) * 11
    phase = 2 * np.pi * (time - time[0]) / rr.PERIOD
    cube = np.ones((len(time), 3, 3)) * 1000
    cube += background[:, None, None]
    cube[600:] += 100
    cube += .3 * time[:, None, None]
    cube[:, 1, 2] += 5 * np.cos(phase) - 3 * np.sin(phase) + 2 * np.cos(2 * phase)
    result = rr.fit_pixels(time, cube, np.ones_like(cube), background)
    assert result["segment_count"] == 2
    injected = next(p for p in result["pixels"] if p["x"] == 2 and p["y"] == 1)
    np.testing.assert_allclose(injected["harmonic_coefficients_electrons_per_s"],
                               [[5, -3], [2, 0], [0, 0]], atol=1e-9)
    for pixel in result["pixels"]:
        if pixel is not injected:
            assert pixel["fundamental_semiamplitude_electrons_per_s"] < 1e-9


def test_sector_guard_rejects_protected_or_nonliteral_sector_before_network():
    for sector in (3, 27, 95, "1", 1., True):
        with pytest.raises(ValueError, match="Literal sector 1/2"):
            rr.fetch(sector)


def test_zero_pixel_uncertainties_are_not_treated_as_measurements():
    time = np.linspace(0, 26, 1000)
    result = rr.fit_pixels(time, np.ones((1000, 1, 1)), np.zeros((1000, 1, 1)), None)
    assert result["pixels"][0]["status"] == "insufficient_positive_uncertainties"
