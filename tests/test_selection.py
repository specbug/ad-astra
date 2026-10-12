import importlib.util
from pathlib import Path

SPEC = importlib.util.spec_from_file_location("select_qlp_pilot", Path(__file__).parents[1] / "scripts" / "select_qlp_pilot.py")
selection = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(selection)


def obs(sector, start):
    return {"sequence_number": sector, "dataURL": f"mast:HLSP/s{sector:04d}/a_llc.fits", "obsid": sector, "obs_id": f"a_{sector}", "provenance_name": "QLP", "t_min": start, "t_max": start + 27}


def test_discovery_remains_sectors_one_two_with_later_season_reserved():
    observations = [obs(1, 0), obs(2, 28), obs(27, 700), obs(28, 728), obs(80, 1800), obs(81, 1828), obs(94, 2200), obs(95, 2228)]
    result, reason = selection.eligible_observations(observations)
    assert reason is None
    assert [row["sector"] for row in result["discovery"]] == [1, 2]
    assert [row["sector"] for row in result["holdout"]] == [80, 81]
    assert 94 not in result["available_safe_sectors"]


def test_latest_season_singleton_is_rejected_not_cherry_picked():
    result, reason = selection.eligible_observations([obs(1, 0), obs(2, 28), obs(27, 700), obs(28, 728), obs(80, 1800)])
    assert result is None
    assert reason == "latest_eligible_season_has_fewer_than_two_sectors"


def test_missing_early_discovery_is_not_substituted():
    result, reason = selection.eligible_observations([obs(2, 28), obs(3, 56), obs(27, 700), obs(28, 728)])
    assert result is None
    assert reason == "missing_sector_1_or_2"


def test_catalog_coordinate_parsing_and_ra_wrap():
    assert selection.degrees("12:00:00", hours=True) == 180
    assert selection.degrees("-66:30:00") == -66.5
    assert selection.angular_separation(359.999, 0, 0.001, 0) < 7.20001
    assert selection.tic_number("TIC 123") == 123
    assert selection.tic_number("123.01") is None
