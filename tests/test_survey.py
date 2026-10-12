"""Frozen-manifest barriers; these tests have no network or science FITS access."""
from copy import deepcopy
from dataclasses import asdict
import importlib.util
import json
from pathlib import Path

import pytest


SPEC = importlib.util.spec_from_file_location("qlp_survey", Path(__file__).parents[1] / "scripts" / "qlp_survey.py")
survey = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(survey)


def product(sector, start, tic=123456789):
    target = f"{tic:016d}"
    path = "/".join(target[i:i+4] for i in range(0, 16, 4))
    return {"sector": sector, "uri": f"mast:HLSP/qlp/s{sector:04d}/{path}/hlsp_qlp_tess_ffi_s{sector:04d}-{target}_tess_v01_llc.fits",
            "obsid": sector+1000, "t_min": start, "t_max": start+27, "provenance_name": "QLP"}


@pytest.fixture
def frozen():
    target = {"tic": 123456789, "discovery": [product(1, 0), product(2, 28)],
              "holdout": [product(27, 700), product(28, 728)]}
    manifest = {"status": survey.FROZEN_STATUS, "review_status": "approved", "selected_count": 1, "targets": [target]}
    config = {"status": survey.FROZEN_STATUS, "review_status": "approved", "manifest_sha256": "manifest_hash",
              "source_sha256": {"source": "hash"}, "science": asdict(survey.ScreenConfig())}
    return manifest, config


def validate(frozen):
    return survey.validate_frozen(*frozen, "manifest_hash", {"source": "hash"})


def test_reviewed_frozen_identity_is_accepted_without_data_access(frozen):
    assert validate(frozen) == survey.ScreenConfig()


@pytest.mark.parametrize("which,field,value", [(0, "status", "provisional"), (1, "review_status", "pending"),
    (1, "manifest_sha256", "wrong_hash"), (1, "source_sha256", {"source": "changed"}), (0, "selected_count", 2)])
def test_unreviewed_or_changed_run_is_rejected(frozen, which, field, value):
    frozen[which][field] = value
    with pytest.raises(ValueError):
        validate(frozen)


@pytest.mark.parametrize("sector", [3, 27, 74, 80, 94, True, "1"])
def test_science_sector_guard_is_literal(sector):
    row = product(1, 0)
    row["sector"] = sector
    with pytest.raises(ValueError, match="limited literally"):
        survey.check_product(123456789, row, {27, 28})


def test_uri_path_and_filename_cannot_point_at_another_identity():
    row = product(1, 0)
    row["uri"] = row["uri"].replace("s0001-", "s0027-")
    with pytest.raises(ValueError, match="sector differs"):
        survey.check_product(123456789, row, {27, 28})
    row = product(1, 0, tic=234567891)
    with pytest.raises(ValueError, match="TIC identity"):
        survey.check_product(123456789, row, {27, 28})
    row = product(1, 0)
    row["uri"] = row["uri"].replace("/0000/0001/", "/0001/0001/")
    with pytest.raises(ValueError, match="TIC identity"):
        survey.check_product(123456789, row, {27, 28})


def test_qlp_provenance_and_finite_span_are_required():
    row = product(1, 0)
    row.pop("provenance_name")
    with pytest.raises(ValueError, match="QLP provenance"):
        survey.check_product(123456789, row, {27, 28})
    row = product(1, 0)
    row["t_max"] = float("inf")
    with pytest.raises(ValueError, match="metadata span"):
        survey.check_product(123456789, row, {27, 28})


def test_repeated_holdout_does_not_satisfy_two_sector_requirement(frozen):
    rows = frozen[0]["targets"][0]["holdout"]
    rows[1] = deepcopy(rows[0])
    with pytest.raises(ValueError, match="must be distinct"):
        validate(frozen)


@pytest.mark.parametrize("sector", [3, 26, 74, 79, 94])
def test_holdout_release_must_be_a_later_validated_sector(frozen, sector):
    frozen[0]["targets"][0]["holdout"][0] = product(sector, 700)
    with pytest.raises(ValueError, match="later validated legacy"):
        validate(frozen)


def test_holdouts_require_identity_and_separated_single_season(frozen):
    rows = frozen[0]["targets"][0]["holdout"]
    rows[0] = product(27, 700, tic=987654321)
    with pytest.raises(ValueError, match="TIC identity"):
        validate(frozen)
    rows[:] = [product(27, 200), product(28, 228)]
    with pytest.raises(ValueError, match="at least 180 days"):
        validate(frozen)
    rows[:] = [product(27, 700), product(28, 1000)]
    with pytest.raises(ValueError, match="one season"):
        validate(frozen)


def test_output_overwrite_barrier_precedes_downloads(tmp_path, frozen, monkeypatch):
    manifest, config = frozen
    config["source_sha256"] = survey.source_hashes()
    manifest_path, config_path = tmp_path/"manifest.json", tmp_path/"config.json"
    manifest_path.write_text(json.dumps(manifest))
    config["manifest_sha256"] = survey.sha256_file(manifest_path)
    config_path.write_text(json.dumps(config))
    output = tmp_path/"results"
    output.mkdir()
    prior = output/"prior.txt"
    prior.write_text("keep this partial run")
    monkeypatch.setattr(survey.ArchiveClient, "download", lambda *args: pytest.fail("download barrier failed"))
    with pytest.raises(ValueError, match="refusing overwrite"):
        survey.run(manifest_path, config_path, output, tmp_path/"cache", workers=1)
    assert prior.read_text() == "keep this partial run"


def test_source_change_prevents_any_target_download(frozen, monkeypatch):
    monkeypatch.setattr(survey.ArchiveClient, "download", lambda *args: pytest.fail("download barrier failed"))
    result = survey.run_target(frozen[0]["targets"][0], asdict(survey.ScreenConfig()), "unused", {"source": "old"})
    assert result["status"] == "failed"
    assert result["sectors"] == []
    assert "Source files changed" in result["error"]
