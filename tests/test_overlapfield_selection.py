"""Guard consequential metadata-selection failure and coverage semantics."""
import importlib.util
import json
from pathlib import Path
import sys

import pytest


scripts = Path(__file__).parents[1] / "scripts"
sys.path.insert(0, str(scripts))
try:
    spec = importlib.util.spec_from_file_location("overlapfield_selector_test", scripts / "select_qlp_overlapfield.py")
    selection = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(selection)
finally:
    sys.path.remove(str(scripts))


def test_recorded_query_failure_is_not_retried_on_resume(tmp_path, monkeypatch):
    def unavailable(*args, **kwargs):
        raise TimeoutError("archive unavailable")

    monkeypatch.setattr(selection.base, "mast_pages", unavailable)
    meta = {"ID": "123", "ra": 7.2, "dec": -60.0}
    first = selection.candidate(1, meta, tmp_path, "declaration", {}, {})
    assert not first["eligible"]
    assert first["excluded_reasons"] == ["observation_query_failed"]

    def would_retry(*args, **kwargs):
        pytest.fail("A completed failed query must remain an exclusion on resume")

    monkeypatch.setattr(selection.base, "mast_pages", would_retry)
    assert selection.candidate(1, meta, tmp_path, "declaration", {}, {}) == first


def test_resume_refuses_changed_raw_metadata(tmp_path, monkeypatch):
    def empty_archive(query, directory, **kwargs):
        selection.write(directory / "response.json", {"data": []})
        return []

    monkeypatch.setattr(selection.base, "mast_pages", empty_archive)
    meta = {"ID": "123", "ra": 7.2, "dec": -60.0}
    result = selection.candidate(1, meta, tmp_path, "declaration", {}, {})
    assert result["excluded_reasons"] == ["missing_sector_1_or_2"]
    (tmp_path / "observations/123/response.json").write_text("changed")
    with pytest.raises(ValueError, match="metadata hash differs"):
        selection.candidate(1, meta, tmp_path, "declaration", {}, {})


def test_empty_gaia_cone_is_unchecked_coverage(tmp_path, monkeypatch):
    def query(adql, path):
        return [{"n": 0}] if "COUNT" in adql else []

    monkeypatch.setattr(selection, "gaia_query", query)
    result = selection.audit_gaia(123, {"ra": 7.2, "dec": -60.0}, tmp_path)
    assert result["target_coverage"] == "unchecked_no_close_gaia_source"
    assert result["excluded_reasons"] == []


def test_gaia_truncation_fails_and_variable_neighbor_excludes(tmp_path, monkeypatch):
    row = {"source_id": 1, "ra": 7.2, "dec": -60.02, "phot_g_mean_mag": 15,
           "phot_variable_flag": "VARIABLE"}
    monkeypatch.setattr(selection, "gaia_query", lambda query, path: [{"n": 2}] if "COUNT" in query else [row.copy()])
    with pytest.raises(ValueError, match="independent count"):
        selection.audit_gaia(123, {"ra": 7.2, "dec": -60.0}, tmp_path)
    monkeypatch.setattr(selection, "gaia_query", lambda query, path: [{"n": 1}] if "COUNT" in query else [row.copy()])
    result = selection.audit_gaia(123, {"ra": 7.2, "dec": -60.0}, tmp_path)
    assert result["excluded_reasons"] == ["gaia_dr3_variable_within_90_arcsec"]
    assert result["variable_neighbors"][0]["separation_arcsec"] == pytest.approx(72)
