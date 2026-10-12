#!/usr/bin/env python3
"""Fixed product-level injection grid on two existing pi Men QLP controls."""
import argparse
from dataclasses import asdict
from pathlib import Path
import json

import numpy as np

from adastra.archive import sha256_file, utc_now, write_json
from adastra.qlp import read_qlp
from adastra.screen import ScreenConfig, clean_curve, screen_curve
import adastra.qlp
import adastra.screen


def compare_truth(peaks, period, epoch, duration):
    """Exact-period recovery also requires the predicted event phase to agree."""
    comparisons = []
    for peak in peaks:
        difference = abs(peak["period_days"] / period - 1)
        phase_offset = abs((peak["epoch_btjd"] - epoch + period/2) % period - period/2)
        exact = difference < .01 and phase_offset < duration
        ratio = peak["period_days"] / period
        alias = min((.3333333333333333, .5, 2., 3.), key=lambda x: abs(ratio/x-1))
        comparisons.append({
            "rank": peak["rank"], "fractional_period_error": difference,
            "phase_offset_days": phase_offset, "recovered_period_and_phase": bool(exact),
            "passes_event_window_triage": peak["passes_event_window_triage"],
            "nearest_alias_ratio": alias, "alias_period_agreement": bool(abs(ratio/alias-1) < .01),
            "alias_interpretation": "period coincidence only; not counted as exact recovery",
        })
    return {"recovered": any(x["recovered_period_and_phase"] for x in comparisons),
            "recovered_and_passes_event_window_triage": any(
                x["recovered_period_and_phase"] and x["passes_event_window_triage"] for x in comparisons),
            "comparisons": comparisons}


def main(args):
    out = args.output
    if (out/"results.json").exists() or (out/"declared-injections.json").exists():
        raise FileExistsError("Choose a new output path; existing declared experiments are immutable")
    config = ScreenConfig()
    controls = json.loads(args.controls.read_text())["results"]
    if {r["sector"] for r in controls} != {1, 27} or len(controls) != 2:
        raise ValueError("This declaration is limited to pi Men QLP sectors 1 and 27")
    sources = {str(Path(p)): sha256_file(p) for p in
               [__file__, adastra.qlp.__file__, adastra.screen.__file__]}
    cases = [{"period_days": period, "depth_fraction": depth,
              "duration_days": .08*period**(1/3), "epoch_phase_from_first_clean_time": .371}
             for period in [1., 3., 7.] for depth in [.001, .003, .01]]
    declaration = {
        "declared_utc": utc_now(), "experiment": "fresh fixed-grid QLP product injections",
        "config": asdict(config), "cases_per_control": cases,
        "input_controls_manifest": str(args.controls), "input_manifest_sha256": sha256_file(args.controls),
        "sources_sha256": sources,
        "injection": "subtract depth * original quality-zero median KSPSAP_FLUX in a box before this code's median division; retain original errors and noise",
        "phase_definition": "epoch = earliest qualifying original sample + 0.371 * period",
        "recovery": "any of the three deterministic peaks within 1% of the injected period and phase offset less than injected full duration",
        "matched_controls": "same uninjected real curve, one per sector, with every input and setting shared",
        "limitations": ["Upstream QLP extraction and detrending are not replayed",
                        "Known pi Men c and other real artifacts remain in these data",
                        "Two sectors of one bright star, one phase per period, and boxes do not establish general completeness",
                        "Uninjected curves are not noise-only null simulations",
                        "Period aliases are recorded but do not count as exact recovery"],
    }
    write_json(out/"declared-injections.json", declaration)
    results = []
    for record in sorted(controls, key=lambda x:x["sector"]):
        if record["tic_id"] != 261136679:
            raise ValueError("Unexpected control identity")
        path = Path(record["local_path"])
        if sha256_file(path) != record["sha256"]:
            raise ValueError("Control input hash differs from original archive verification")
        lc = read_qlp(path)
        t, _, _, clean = clean_curve(lc.time_btjd, lc.flux, lc.flux_error, lc.quality)
        original = screen_curve(lc.time_btjd, lc.flux, lc.flux_error, lc.quality, config)
        known = compare_truth(original["bls_peaks"], 6.2678399, 1425.789204, 2.952024/24)
        row = {"sector": record["sector"], "tic_id": 261136679,
               "input_path": str(path), "input_sha256": record["sha256"],
               "uninjected": original, "known_planet_recovery": known, "cases": []}
        for case in cases:
            for source, expected in sources.items():
                if sha256_file(source) != expected:
                    raise RuntimeError(f"Source changed after experiment declaration: {source}")
            period, duration, depth = case["period_days"], case["duration_days"], case["depth_fraction"]
            epoch = float(t[0] + .371*period)
            phase = (lc.time_btjd-epoch+period/2)%period-period/2
            inside = np.abs(phase) < duration/2
            injected_flux = lc.flux.copy()
            injected_flux[inside] -= depth*clean["median_flux_normalization"]
            injected = screen_curve(lc.time_btjd, injected_flux, lc.flux_error, lc.quality, config)
            if injected["cleaning"]["analysis_samples"] != clean["analysis_samples"]:
                raise ValueError("Injection changed sample eligibility; this experiment forbids selection changes")
            recovery = compare_truth(injected["bls_peaks"], period, epoch, duration)
            row["cases"].append({**case, "epoch_btjd": epoch, "result": injected, **recovery})
            write_json(out/"partial-results.json", {"declaration": declaration, "completed_sectors": results, "current_sector": row})
            print(json.dumps({"sector":record["sector"],"period":period,"depth":depth,
                              "recovered":recovery["recovered"],
                              "event_window_triage":recovery["recovered_and_passes_event_window_triage"]}), flush=True)
        results.append(row)
    for source, expected in sources.items():
        if sha256_file(source) != expected:
            raise RuntimeError(f"Source changed after experiment declaration: {source}")
    cases_all = [c for r in results for c in r["cases"]]
    summary = {"cases": len(cases_all), "exact_period_and_phase_recoveries": sum(c["recovered"] for c in cases_all),
               "recoveries_passing_event_window_triage": sum(c["recovered_and_passes_event_window_triage"] for c in cases_all),
               "by_depth": [{"depth_fraction":d,"cases":sum(c["depth_fraction"]==d for c in cases_all),
                             "recovered":sum(c["recovered"] and c["depth_fraction"]==d for c in cases_all)}
                            for d in [.001,.003,.01]]}
    write_json(out/"results.json", {"declaration":declaration,"completed_utc":utc_now(),"results":results,"summary":summary})
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--controls", type=Path, default=Path("runs/reconstruction/qlp-control/results.json"))
    parser.add_argument("--output", type=Path, default=Path("runs/reconstruction/qlp-injections"))
    main(parser.parse_args())
