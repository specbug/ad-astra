#!/usr/bin/env python3
"""Run a reviewed, frozen QLP exploratory pilot; discovery is sectors 1/2 only.

The script intentionally has no automatic retries or resume/overwrite option.
Its only network operation is retrieving the exact frozen discovery products.
"""
import os

# Set these before numerical imports, including in spawned subprocesses.
for _thread_option in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[_thread_option] = "1"

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict
import hashlib
import importlib.metadata
import json
import math
import multiprocessing
from pathlib import Path
import platform
import re
import numpy as np

from adastra.archive import ArchiveClient, sha256_file, utc_now, write_json
from adastra.qlp import read_qlp
from adastra.screen import ScreenConfig, clean_curve, compare_sectors, screen_curve


ROOT = Path(__file__).resolve().parents[1]
SOURCE_FILES = ("scripts/qlp_survey.py", "src/adastra/archive.py", "src/adastra/qlp.py", "src/adastra/screen.py")
FROZEN_STATUS = "frozen_for_exploratory_run"
QLP_URI = re.compile(r"mast:HLSP/qlp/s(?P<sector>\d{4})/(?P<a>\d{4})/(?P<b>\d{4})/(?P<c>\d{4})/(?P<d>\d{4})/hlsp_qlp_tess_ffi_s(?P<file_sector>\d{4})-(?P<tic>\d{16})_tess_v\d+_llc\.fits")


def digest_object(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()).hexdigest()


def source_hashes(root=ROOT):
    return {name: sha256_file(root / name) for name in SOURCE_FILES}


def check_product_identity(tic, product):
    """Validate a metadata identity without opening or downloading its product."""
    sector = product.get("sector")
    if type(sector) is not int:
        raise ValueError("QLP sector must be an integer")
    if product.get("provenance_name") != "QLP":
        raise ValueError("Frozen archive metadata must explicitly identify QLP provenance")
    match = QLP_URI.fullmatch(product.get("uri", ""))
    if match is None:
        raise ValueError("Unsupported QLP product URI")
    if int(match["sector"]) != sector or int(match["file_sector"]) != sector:
        raise ValueError("QLP URI sector differs from the frozen identity")
    if int(match["tic"]) != tic or int("".join(match[k] for k in "abcd")) != tic:
        raise ValueError("QLP URI target differs from the frozen TIC identity")
    start, end = product.get("t_min"), product.get("t_max")
    if (not product.get("obsid") or type(start) not in (int, float) or type(end) not in (int, float)
            or not math.isfinite(start) or not math.isfinite(end) or not 0 < end - start <= 40):
        raise ValueError("Missing observation identity or invalid per-sector metadata span")
    return digest_object({"tic": tic, "product": product})


def check_product(tic, product, holdout_sectors):
    sector = product.get("sector")
    if type(sector) is not int or sector not in (1, 2):
        raise ValueError("Science retrieval is limited literally to discovery sectors 1 and 2")
    if sector in holdout_sectors:
        raise ValueError("Discovery product overlaps a protected holdout sector")
    return check_product_identity(tic, product)


def validate_frozen(manifest, config, manifest_sha256, current_source_hashes):
    for name, value in (("manifest", manifest), ("configuration", config)):
        if value.get("status") != FROZEN_STATUS or value.get("review_status") != "approved":
            raise ValueError(f"The {name} is provisional or has no approved scientific review")
    if config.get("manifest_sha256") != manifest_sha256:
        raise ValueError("Configuration does not bind this exact frozen manifest")
    if config.get("source_sha256") != current_source_hashes:
        raise ValueError("Current source files differ from the frozen reviewed configuration")
    science = ScreenConfig(**config["science"])
    targets = manifest.get("targets", [])
    if not targets or manifest.get("selected_count") != len(targets):
        raise ValueError("Manifest target count is empty or inconsistent")
    seen_tics, seen_uris = set(), set()
    for target in targets:
        tic = target.get("tic")
        if type(tic) is not int or tic <= 0 or tic in seen_tics:
            raise ValueError("Duplicate or invalid TIC identity")
        seen_tics.add(tic)
        holdouts = target.get("holdout", [])
        if len(holdouts) < 2:
            raise ValueError("At least two explicitly protected holdout products are required")
        held_sectors = {row["sector"] for row in holdouts}
        held_uris = {row["uri"] for row in holdouts}
        if len(held_sectors) != len(holdouts) or len(held_uris) != len(holdouts):
            raise ValueError("Protected holdout sectors and URIs must be distinct")
        for product in holdouts:
            check_product_identity(tic, product)
            if not (27 <= product["sector"] <= 73 or 80 <= product["sector"] <= 93):
                raise ValueError("Protected holdouts must be later validated legacy QLP sectors")
        discovery = target.get("discovery", [])
        if len(discovery) != 2 or {row["sector"] for row in discovery} != {1, 2}:
            raise ValueError("The discovery set must be exactly sectors 1 and 2")
        for product in discovery:
            check_product(tic, product, held_sectors)
            uri = product["uri"]
            if uri in held_uris or uri in seen_uris:
                raise ValueError("Discovery URI is repeated or overlaps a protected holdout")
            seen_uris.add(uri)
        ordered_holdouts = sorted(holdouts, key=lambda row: row["t_min"])
        if ordered_holdouts[0]["t_min"] - max(row["t_max"] for row in discovery) < 180:
            raise ValueError("Protected holdouts must start at least 180 days after discovery")
        if any(later["t_min"] - earlier["t_max"] >= 180
               for earlier, later in zip(ordered_holdouts, ordered_holdouts[1:])):
            raise ValueError("Protected holdouts must belong to one season under the 180-day gap rule")
    return science


def run_target(target, science_dict, cache, expected_sources):
    tic = target["tic"]
    result = {"tic": tic, "started_utc": utc_now(), "sectors": [], "comparison": None}
    if source_hashes() != expected_sources:
        result.update(status="failed", error="Source files changed after declaration", completed_utc=utc_now())
        return result
    client = ArchiveClient(cache)
    held_sectors = {row["sector"] for row in target["holdout"]}
    for product in sorted(target["discovery"], key=lambda row: row["sector"]):
        row = {"sector": product["sector"], "uri": product["uri"], "obsid": product["obsid"]}
        try:
            row["frozen_product_identity_sha256"] = check_product(tic, product, held_sectors)
            path = client.download(product["uri"])
            input_hash = sha256_file(path)
            lc = read_qlp(path)
            if lc.metadata["tic_id"] != tic or lc.metadata["sector"] != product["sector"]:
                raise ValueError("Downloaded FITS TIC/sector differs from the frozen product identity")
            provenance = json.loads(path.with_suffix(".json").read_text())
            if provenance.get("uri") != product["uri"] or provenance.get("sha256") != input_hash:
                raise ValueError("Downloaded provenance does not bind this URI and file hash")
            row.update(input_path=str(path), input_sha256=input_hash, ingestion=lc.metadata,
                       download_provenance=provenance, analysis=screen_curve(lc.time_btjd, lc.flux, lc.flux_error, lc.quality, science_dict))
            # Same qualifying rows, reflected about their median before the
            # additional running-median step. This tests brightening symmetry;
            # it is not an independent null sample or false-positive rate.
            try:
                t, f, e, cleaning = clean_curve(lc.time_btjd, lc.flux, lc.flux_error, lc.quality)
                inverted = 2.0 - f
                if np.any(inverted <= 0):
                    raise ValueError("Reflection is nonpositive; refusing a different control sample")
                row["inverted_analysis"] = screen_curve(t, inverted, e, np.zeros(len(t), dtype=int), science_dict)
                row["inversion"] = {"status": "complete", "reflection": "2 - median-normalized flux before additional detrending",
                                    "original_cleaning": cleaning, "same_input_rows": True,
                                    "interpretation": "Matched sign control; correlated noise and upstream corrections remain; not a null FAP"}
            except Exception as error:
                row["inversion"] = {"status": "failed", "error": f"{type(error).__name__}: {error}"}
            if sha256_file(path) != input_hash:
                raise ValueError("Input changed during analysis")
            row["status"] = "complete"
        except Exception as error:
            row.update(status="failed", error=f"{type(error).__name__}: {error}")
        result["sectors"].append(row)
    if all(row["status"] == "complete" for row in result["sectors"]):
        result["comparison"] = compare_sectors(*(row["analysis"] for row in result["sectors"]), config=science_dict)
        if all(row.get("inversion", {}).get("status") == "complete" for row in result["sectors"]):
            result["inverted_comparison"] = compare_sectors(*(row["inverted_analysis"] for row in result["sectors"]), config=science_dict)
        result["status"] = "complete"
    else:
        result["status"] = "failed"
    result["completed_utc"] = utc_now()
    return result


def run(manifest_path, config_path, output, cache, workers=4):
    if not 1 <= workers <= 4:
        raise ValueError("Use one to four process workers")
    manifest_path, config_path, output, cache = map(Path, (manifest_path, config_path, output, cache))
    manifest, config = (json.loads(path.read_text()) for path in (manifest_path, config_path))
    hashes = source_hashes()
    science = validate_frozen(manifest, config, sha256_file(manifest_path), hashes)
    # This check precedes any write or download. Existing partial runs remain
    # intact and need an explicitly reviewed recovery, not an implicit retry.
    if output.exists() and any(output.iterdir()):
        raise ValueError("Output directory is nonempty; refusing overwrite or hidden retry")
    declaration = {"created_utc": utc_now(), "status": "declared_before_discovery_retrieval",
                   "manifest_path": str(manifest_path), "manifest_sha256": sha256_file(manifest_path),
                   "configuration_path": str(config_path), "configuration_sha256": sha256_file(config_path),
                   "source_sha256": hashes, "science": asdict(science), "workers": workers,
                   "numerical_threads_per_worker": 1, "python": platform.python_version(),
                   "dependencies": {package: importlib.metadata.version(package) for package in ("numpy", "scipy", "astropy", "requests")},
                   "allowed_discovery_sectors": [1, 2], "protected_holdouts_accessed": False,
                   "matched_sign_control": "Every usable curve reflected about its median before additional detrending; identical screen; no null FAP",
                   "stage": "new reconstructed exploratory QLP pilot; historical cohort overlap is unknown",
                   "retry_policy": "One attempt per frozen product; record every error, no automatic retry"}
    write_json(output / "declaration.json", declaration)
    write_json(output / "frozen-manifest.json", manifest)
    write_json(output / "frozen-config.json", config)
    results = []
    payload = {"declaration": declaration, "status": "running", "results": results,
               "target_count": len(manifest["targets"]), "completed_targets": 0}
    write_json(output / "results.json", payload)
    with ProcessPoolExecutor(max_workers=workers, mp_context=multiprocessing.get_context("spawn")) as executor:
        futures = {executor.submit(run_target, target, asdict(science), str(cache), hashes): target for target in manifest["targets"]}
        for future in as_completed(futures):
            target = futures[future]
            try:
                row = future.result()
            except Exception as error:
                row = {"tic": target["tic"], "status": "failed", "error": f"Worker {type(error).__name__}: {error}",
                       "sectors": [], "comparison": None, "completed_utc": utc_now()}
            results.append(row)
            results.sort(key=lambda row: row["tic"])
            payload["completed_targets"] = len(results)
            write_json(output / "targets" / f"{target['tic']}.json", row)
            write_json(output / "results.json", payload)
            print(f"TIC {target['tic']}: {row['status']} ({len(results)}/{len(manifest['targets'])})", flush=True)
    unchanged = hashes == source_hashes()
    payload.update(status="complete" if unchanged else "invalid_source_changed", completed_utc=utc_now(),
                   successful_targets=sum(row["status"] == "complete" for row in results),
                   failed_targets=sum(row["status"] != "complete" for row in results),
                   input_curves_completed=sum(sector["status"] == "complete" for row in results for sector in row["sectors"]),
                   input_curves_failed=sum(sector["status"] != "complete" for row in results for sector in row["sectors"]),
                   sign_controls_completed=sum(sector.get("inversion", {}).get("status") == "complete" for row in results for sector in row["sectors"]),
                   original_sources_unchanged=unchanged,
                   all_requested_targets_accounted_for={row["tic"] for row in results} == {row["tic"] for row in manifest["targets"]})
    write_json(output / "results.json", payload)
    return payload


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--cache", type=Path, default=Path("data/cache"))
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    run(args.manifest, args.config, args.output, args.cache, args.workers)
