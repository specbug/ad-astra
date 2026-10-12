#!/usr/bin/env python3
"""Bind the reviewed metadata pilot, validated controls and exact source files.

This performs no network requests and opens no science light curves. Review
approval records the research workflow decision, not a discovery certification.
"""
from copy import deepcopy
from dataclasses import asdict
import json
from pathlib import Path

from adastra.archive import sha256_file, utc_now, write_json
from adastra.screen import ScreenConfig
from qlp_survey import FROZEN_STATUS, source_hashes, validate_frozen


def enrich_product(tic, product):
    matches = []
    for root in (Path("data/qlp-pilot"), Path("data/qlp-pilot-v2"), Path("data/qlp-pilot-v3")):
        for path in sorted((root / "observations" / str(tic)).glob("response-[0-9][0-9][0-9].json")):
            provenance = json.loads(path.with_name(path.name + ".provenance.json").read_text())
            if sha256_file(path) != provenance["sha256"]:
                raise ValueError(f"Metadata snapshot hash mismatch: {path}")
            response = json.loads(path.read_text())
            if response.get("status") != "COMPLETE":
                raise ValueError("Incomplete metadata snapshot")
            for row in response["data"]:
                if row["obsid"] == product["obsid"] and row["dataURL"] == product["uri"]:
                    if (row["provenance_name"] != "QLP" or int(row["target_name"]) != tic
                            or int(row["sequence_number"]) != product["sector"]
                            or row["t_min"] != product["t_min"] or row["t_max"] != product["t_max"]):
                        raise ValueError("Frozen product identity conflicts with raw metadata")
                    matches.append((row, path, provenance))
    if not matches:
        raise ValueError(f"No verified raw observation metadata for TIC {tic}")
    row, path, provenance = matches[0]
    product.update(provenance_name=row["provenance_name"], metadata_response=str(path),
                   metadata_response_sha256=provenance["sha256"], metadata_retrieved_utc=provenance["retrieved_utc"])


def main():
    provisional = Path("data/qlp-pilot-v3/provisional-manifest.json")
    output = Path("data/qlp-frozen")
    if output.exists() and any(output.iterdir()):
        raise ValueError("Refusing to overwrite an existing freeze")
    manifest = deepcopy(json.loads(provisional.read_text()))
    if manifest["selected_count"] != 24 or manifest["flux_accessed"] is not False or not manifest.get("selection_wrapper_complete"):
        raise ValueError("Expected the complete 24-target metadata-only selection")
    controls_path = Path("runs/reconstruction/qlp-injections/results.json")
    controls = json.loads(controls_path.read_text())
    for path, digest in controls["declaration"]["sources_sha256"].items():
        if sha256_file(Path(path)) != digest:
            raise ValueError(f"Source differs from the completed control experiment: {path}")
    if controls["summary"]["cases"] != 18 or controls["summary"]["exact_period_and_phase_recoveries"] != 18:
        raise ValueError("Control results differ from the reviewed experiment")
    manifest.update(status=FROZEN_STATUS, review_status="approved", frozen_utc=utc_now(),
                    provisional_manifest_sha256=sha256_file(provisional),
                    review_meaning="Approved for exploratory screening after code/control review; no detection or novelty claim")
    for target in manifest["targets"]:
        for product in target["discovery"] + target["holdout"]:
            enrich_product(target["tic"], product)
    config = {
        "status": FROZEN_STATUS, "review_status": "approved", "frozen_utc": utc_now(),
        "science": asdict(ScreenConfig()), "source_sha256": source_hashes(),
        "control_results": str(controls_path), "control_results_sha256": sha256_file(controls_path),
        "reviewed_control_outcomes": {
            "injections": "18/18 recovered in period and phase: 1, 3, 7 days; 1000, 3000, 10000 ppm; two pi Mensae products",
            "uninjected_pi_mensae_c": "Recovered in sector 1, missed in sector 27 even in the three-peak screen",
            "limitation": "Only deep injected signals on one bright star; not population completeness or upstream pipeline recovery",
        },
        "matched_control": "Reflect every cleaned median-normalized curve as 2-flux before further detrending; same screen and pair checks",
        "priority_rules": "All up-to-nine period pairs; direct or 2/3 harmonic agreement within 1 percent; >=3 sampled and positive-depth windows; LS strongest-period agreement within 1 percent",
        "interpretation": "Exploratory priorities only. Positive windows are not significant events. Search conditional errors omit red noise and selection.",
        "validation_policy": "Science flux retrieval restricted to sectors 1/2 because the lost original manifest prevents auditing old holdouts. No final season access.",
        "novelty_policy": "NASA TOI/confirmed and ExoFOP CTOI exclusions plus scoped SIMBAD types/references; additional literature checks may find known signals. ADS excluded. No exhaustive novelty clearance.",
        "failure_policy": "Retain every failure; no automatic retries or replacement stars; require explicit review of any amendment",
    }
    # Validate the entire declaration in memory before writing the pair.
    from hashlib import sha256
    serialized = json.dumps(manifest, indent=2, allow_nan=False) + "\n"
    config["manifest_sha256"] = sha256(serialized.encode()).hexdigest()
    validate_frozen(manifest, config, config["manifest_sha256"], source_hashes())
    write_json(output / "manifest.json", manifest)
    if sha256_file(output / "manifest.json") != config["manifest_sha256"]:
        raise ValueError("Manifest serialization did not preserve the reviewed hash")
    write_json(output / "config.json", config)
    print(f"Frozen {manifest['selected_count']} targets; metadata only; {output}")


if __name__ == "__main__":
    main()
