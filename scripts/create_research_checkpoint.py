#!/usr/bin/env python3
"""Preserve the combined completed QLP research without environment or credential files.

This command performs no network requests, science fitting, or tests. The external
inventory describes every archive member; SHA256SUMS protects the inventory too.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import io
import json
from pathlib import Path
import tarfile


ROOT_FILES = {"README.md", ".gitignore", "pyproject.toml", "requirements-validated.txt"}
ROOT_DIRS = {"src", "tests", "scripts", "docs", "data", "runs"}
SKIP_PARTS = {".git", ".venv", "__pycache__", ".pytest_cache"}
EXCLUDE_PREFIXES = ()
EXCLUDE_FILES = set()
FIGURE_SUFFIXES = {".png", ".pdf", ".jpg", ".jpeg", ".svg", ".webp"}


def now() -> str:
    return datetime.now(timezone.utc).isoformat()


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_json(path: Path):
    return json.loads(path.read_text())


def audit_single_pilot(root: Path, run_name: str) -> dict:
    run = root / "runs" / run_name
    results = read_json(run / "results.json")
    expected = {"status": "complete", "target_count": 24, "completed_targets": 24,
                "successful_targets": 24, "failed_targets": 0,
                "input_curves_completed": 48, "input_curves_failed": 0,
                "sign_controls_completed": 48, "original_sources_unchanged": True,
                "all_requested_targets_accounted_for": True}
    for key, value in expected.items():
        if results.get(key) != value:
            raise ValueError(f"Completed pilot invariant failed: {key}")
    config = read_json(run / "frozen-config.json")
    manifest = read_json(run / "frozen-manifest.json")
    if digest(run / "frozen-manifest.json") != config["manifest_sha256"]:
        raise ValueError("Frozen manifest hash changed")
    for relative, expected_sha in config["source_sha256"].items():
        for candidate in [root / relative, run / "source-snapshot" / relative]:
            if digest(candidate) != expected_sha:
                raise ValueError(f"Frozen source hash changed: {candidate}")
    if digest(root / config["control_results"]) != config["control_results_sha256"]:
        raise ValueError("Pre-screen injection results changed")
    science_paths = set()
    controls = 0
    for target in results["results"]:
        if target["status"] != "complete" or len(target["sectors"]) != 2:
            raise ValueError("Incomplete target")
        for sector in target["sectors"]:
            if sector["sector"] not in (1, 2) or sector["status"] != "complete":
                raise ValueError("Unapproved or incomplete science sector")
            if sector["analysis"]["status"] != "complete" or sector["inverted_analysis"]["status"] != "complete":
                raise ValueError("Incomplete original or sign control")
            path = root / sector["input_path"]
            if digest(path) != sector["input_sha256"]:
                raise ValueError(f"Science input changed: {path}")
            science_paths.add(sector["input_path"])
            controls += 1
    if len(science_paths) != 48 or controls != 48:
        raise ValueError("Expected 48 distinct science inputs and sign controls")
    protected = {p["uri"] for target in manifest["targets"] for p in target["holdout"]}
    cached_uris = {read_json(path)["uri"] for path in (root / "data/cache").glob("*.json")}
    if protected & cached_uris:
        raise ValueError("Protected later-sector light curve found in cache")
    return {"audited_utc": now(), "checks": expected,
            "science_input_sha256_verified": 48, "sign_controls_verified": controls,
            "current_and_snapshotted_source_sha256_verified": config["source_sha256"],
            "named_protected_uris_cached": 0,
            "limitation": "Cached QLP light curves checked; only previously recorded sector-1/2 pixel data are included. No tests rerun."}


def audit_pilot(root: Path) -> dict:
    names = ("qlp-pilot", "qlp-overlapfield-pilot")
    audits = {name: audit_single_pilot(root, name) for name in names}
    inputs = {}
    targets = set()
    for name in names:
        run = root / "runs" / name
        results = read_json(run / "results.json")
        manifest = read_json(run / "frozen-manifest.json")
        config = read_json(run / "frozen-config.json")
        # The archived declaration binds the exact frozen config and manifest.
        declaration = read_json(run / "declaration.json")
        for key, file in (("manifest_sha256", "frozen-manifest.json"),
                          ("configuration_sha256", "frozen-config.json")):
            if declaration.get(key) != digest(run / file):
                raise ValueError(f"Declaration binding changed: {name}/{key}")
        if declaration.get("source_sha256") != config["source_sha256"]:
            raise ValueError(f"Declaration source binding changed: {name}")
        for recorded_path, recorded_hash in ((declaration["manifest_path"], declaration["manifest_sha256"]),
                                             (declaration["configuration_path"], declaration["configuration_sha256"])):
            if digest(root / recorded_path) != recorded_hash:
                raise ValueError(f"Original frozen input changed: {recorded_path}")
        if manifest.get("status") != "frozen_for_exploratory_run" or config.get("review_status") != "approved":
            raise ValueError(f"Pilot was not frozen and approved: {name}")
        for target in manifest["targets"]:
            tic = str(target["tic"])
            if tic in targets:
                raise ValueError("Target appears in both pilot manifests")
            targets.add(tic)
        for target in results["results"]:
            for sector in target["sectors"]:
                path = sector["input_path"]
                if path in inputs:
                    raise ValueError("Science input reused across pilots")
                inputs[path] = sector["input_sha256"]
    if len(targets) != 48 or len(inputs) != 96 or len(set(inputs.values())) != 96:
        raise ValueError("Combined checkpoint requires 48 targets and 96 distinct science inputs")
    return {"audited_utc": now(), "pilots": audits, "distinct_targets": len(targets),
            "distinct_science_inputs": len(inputs), "sign_controls_verified": 96,
            "science_input_sha256": inputs, "named_protected_uris_cached": 0}


def selected_files(root: Path):
    for path in sorted(root.rglob("*")):
        relative = path.relative_to(root)
        text = relative.as_posix()
        if relative.parts[0] not in ROOT_DIRS and text not in ROOT_FILES:
            continue
        if any(part in SKIP_PARTS or part.endswith(".egg-info") for part in relative.parts):
            continue
        if any(text == p or text.startswith(p + "/") for p in EXCLUDE_PREFIXES) or text in EXCLUDE_FILES:
            continue
        if path.suffix in {".pyc", ".part"}:
            continue
        if path.is_symlink():
            raise ValueError(f"Unexpected symlink in checkpoint scope: {relative}")
        if not path.is_file():
            continue
        if any(part.startswith(".env") or part in {".netrc", "credentials", "secrets"} for part in relative.parts):
            raise ValueError(f"Unexpected credential-like filename in checkpoint scope: {relative}")
        yield path


def verify_archive(archive: Path, entries: list[dict]) -> int:
    expected = {e["archive_path"]: e for e in entries}
    seen = set()
    with tarfile.open(archive, "r:gz") as tar:
        for member in tar:
            if member.name not in expected or member.name in seen or not member.isfile():
                raise ValueError(f"Unexpected archive member: {member.name}")
            entry = expected[member.name]
            handle = tar.extractfile(member)
            if handle is None:
                raise ValueError("Missing archive member bytes")
            h = hashlib.sha256()
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                h.update(block)
            if member.size != entry["size_bytes"] or h.hexdigest() != entry["sha256"]:
                raise ValueError(f"Archive member does not match inventory: {member.name}")
            seen.add(member.name)
    if seen != set(expected):
        raise ValueError("Archive is missing inventoried files")
    return len(seen)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--output", type=Path, default=Path("/workspace/shared/downloads"))
    parser.add_argument("--prefix", default="ad-astra-research-20261012")
    args = parser.parse_args()
    if Path(args.prefix).name != args.prefix:
        parser.error("prefix must be a filename component")
    root = args.root.resolve()
    audit = audit_pilot(root)
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    original_checkpoint = {str(path): digest(path) for path in sorted(output.glob("ad-astra-qlp-pilot-20261012-*")) if path.is_file()}
    destinations = [output / f"{args.prefix}-{part}.tar.gz" for part in ("core", "figures")]
    inventory_path = output / f"{args.prefix}-inventory.json"
    sums_path = output / f"{args.prefix}-SHA256SUMS.txt"
    if any(p.exists() or p.with_suffix(p.suffix + ".part").exists() for p in destinations + [inventory_path, sums_path]):
        raise FileExistsError("Checkpoint output already exists; choose a new prefix")
    entries = []
    for path in selected_files(root):
        stat = path.stat()
        relative = path.relative_to(root).as_posix()
        entries.append({"archive_path": "ad-astra/" + relative,
                        "source_path": str(path), "size_bytes": stat.st_size,
                        "source_modified_utc": datetime.fromtimestamp(stat.st_mtime, timezone.utc).isoformat(),
                        "sha256": digest(path),
                        "archive_group": "figures" if path.suffix.lower() in FIGURE_SUFFIXES else "core"})
    archives = []
    for group, destination in zip(("core", "figures"), destinations):
        members = [entry for entry in entries if entry["archive_group"] == group]
        partial = destination.with_suffix(destination.suffix + ".part")
        with tarfile.open(partial, "w:gz", compresslevel=6) as tar:
            for entry in members:
                source = Path(entry["source_path"])
                content = source.read_bytes()
                if len(content) != entry["size_bytes"] or hashlib.sha256(content).hexdigest() != entry["sha256"]:
                    raise ValueError(f"File changed while checkpointing: {source}")
                info = tarfile.TarInfo(entry["archive_path"])
                info.size = len(content)
                info.mode = 0o755 if source.stat().st_mode & 0o111 else 0o644
                info.mtime = int(source.stat().st_mtime)
                tar.addfile(info, io.BytesIO(content))
        verified = verify_archive(partial, members)
        partial.rename(destination)
        archives.append({"path": str(destination), "group": group,
                         "sha256": digest(destination), "size_bytes": destination.stat().st_size,
                         "verified_file_count": verified,
                         "uncompressed_member_bytes": sum(entry["size_bytes"] for entry in members)})
    for entry in entries:
        if digest(Path(entry["source_path"])) != entry["sha256"]:
            raise ValueError(f"Source changed before completion: {entry['source_path']}")
    for preserved_path, preserved_sha in original_checkpoint.items():
        if digest(Path(preserved_path)) != preserved_sha:
            raise ValueError(f"First checkpoint changed: {preserved_path}")
    inventory = {"original_checkpoint_preserved_sha256": original_checkpoint, "created_utc": now(), "source_root": str(root), "scope": "Both completed TESS QLP pilots and all current public research artifacts",
                 "excluded": {"directories": sorted(SKIP_PARTS), "suffixes": [".egg-info", ".pyc", ".part"],
                              "prefixes": list(EXCLUDE_PREFIXES), "files": sorted(EXCLUDE_FILES),
                              "root_allowlist": sorted(ROOT_FILES | ROOT_DIRS)},
                 "pilot_integrity_audit": audit, "archives": archives,
                 "member_verification": "Every regular-file member read back from each archive and matched byte size and SHA256; no extra or duplicate members",
                 "source_verification": "Every included source rehashed after both archives were written",
                 "files": entries,
                 "limitations": ["Local archives can be lost if the environment is replaced. Retain the downloads elsewhere.",
                                 "No .git history, virtual environment, credentials, or configuration secrets included.",
                                 "The earlier lost 88-star artifacts are absent. Named protected later-sector light curves remain unopened.",
                                 "SHA256 establishes file identity, not the scientific validity or novelty of results."]}
    inventory_path.write_text(json.dumps(inventory, indent=2) + "\n")
    sums_path.write_text("".join(f"{digest(path)}  {path.name}\n" for path in destinations + [inventory_path]))
    print(json.dumps({"archives": archives, "inventory": str(inventory_path),
                      "checksums": str(sums_path), "file_count": len(entries)}, indent=2))


if __name__ == "__main__":
    main()
