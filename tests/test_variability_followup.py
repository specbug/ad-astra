import importlib.util
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from qlp_variability_followup import fit, EPOCH


def test_sinusoid_recovery_with_an_isolated_sample_and_segment_trends():
    period = .477
    t = np.r_[np.arange(1325., 1335., .02), [1337.], np.arange(1340., 1350., .02)]
    baseline = np.where(t < 1336, 23*(t-1330)+180, -12*(t-1345)-50)
    baseline[t == 1337] = 900
    y = 283*np.sin(2*np.pi*(t-EPOCH)/period+.3) + baseline
    result, _ = fit(t, y, period)
    assert abs(result['amplitude']-283) < 1e-7
    assert abs(result['phase_cycles_at_btjd1340']-.3/(2*np.pi)) < 1e-9
    assert result['samples'] == len(t)


def test_robust_fit_reduces_single_extreme_artifact_without_dropping_samples():
    t = np.arange(1325, 1352, .02)
    period = 1.194
    y = 120*np.sin(2*np.pi*(t-EPOCH)/period) + np.random.default_rng(827).normal(0, 30, len(t))
    y[35] += 100000
    ordinary, _ = fit(t, y, period)
    robust, _ = fit(t, y, period, robust=True)
    assert abs(robust['amplitude']-120) < 5
    assert abs(ordinary['amplitude']-120) > 20
    assert robust['samples'] == len(t)
    assert robust['converged']
