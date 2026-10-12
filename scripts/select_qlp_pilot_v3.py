#!/usr/bin/env python3
"""Apply recorded prior-literature exclusions to the preserved metadata selection.

This wrapper adds no flux reads. It keeps the v2 selector and its history exclusions
unchanged, reuses verified raw metadata, and records the additional rule before
running the inherited deterministic selector.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil
import sys

import select_qlp_pilot_v2 as selector

base = selector.base


def load_exclusions(path):
    artifact = json.loads(path.read_text())
    if not isinstance(artifact.get("rules"), list) or not artifact["rules"]:
        raise ValueError("Literature artifact must declare its exclusion rules")
    if not isinstance(artifact.get("sources"), list) or not artifact["sources"]:
        raise ValueError("Literature artifact must cite preserved source files")
    for source in artifact["sources"]:
        source_path = Path(source["path"])
        if not source_path.is_file() or base.sha(source_path) != source["sha256"]:
            raise ValueError(f"Literature source hash mismatch: {source_path}")
    exclusions = artifact.get("exclusions")
    if not isinstance(exclusions, dict):
        raise ValueError("Literature exclusions must map TIC identities to reasons")
    normalized = {}
    for identity, reasons in exclusions.items():
        if not str(identity).isdigit() or int(identity) <= 0:
            raise ValueError("Literature exclusion has an invalid TIC identity")
        if not isinstance(reasons, list) or not reasons or not all(isinstance(reason, str) and reason for reason in reasons):
            raise ValueError(f"Literature exclusion reasons are invalid: {identity}")
        normalized[int(identity)] = sorted(set(reasons))
    return artifact, normalized


def copy_verified_metadata(source, destination):
    """Copy intact raw responses and sidecars without overwriting prior outputs."""
    records = []
    for sidecar in sorted(source.rglob("*.provenance.json")):
        relative = sidecar.relative_to(source)
        if relative.parts[0] not in {"simbad", "observations", "simbad-type-definitions.json.provenance.json"}:
            continue
        raw = sidecar.with_name(sidecar.name.removesuffix(".provenance.json"))
        provenance = json.loads(sidecar.read_text())
        if not raw.is_file() or base.sha(raw) != provenance.get("sha256") or not isinstance(provenance.get("url"), str):
            raise ValueError(f"Metadata cache provenance mismatch: {raw}")
        for origin in (raw, sidecar):
            target = destination / origin.relative_to(source)
            digest = base.sha(origin)
            if target.exists():
                if base.sha(target) != digest:
                    raise ValueError(f"Refusing to overwrite a changed metadata cache: {target}")
            else:
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(origin, target)
            records.append({"source": str(origin), "copy": str(target), "sha256": digest})
    return records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--previous", type=Path, default=Path("data/qlp-pilot"))
    parser.add_argument("--metadata-cache", type=Path, default=Path("data/qlp-pilot-v2"))
    parser.add_argument("--output", type=Path, default=Path("data/qlp-pilot-v3"))
    parser.add_argument("--literature-exclusions", type=Path, default=Path("data/qlp-prior-review/literature-exclusions.json"))
    parser.add_argument("--count", type=int, default=24)
    args = parser.parse_args()
    if args.count <= 0:
        raise ValueError("Target count must be positive")
    if args.output.resolve() in {args.previous.resolve(), args.metadata_cache.resolve()}:
        raise ValueError("V3 output must preserve both prior selection directories")
    completed_cache = args.metadata_cache / "provisional-manifest.json"
    previous = json.loads(completed_cache.read_text())
    if previous.get("selected_count") != previous.get("requested_count"):
        raise ValueError("Wait for the prior metadata selection to complete")
    artifact, exclusions = load_exclusions(args.literature_exclusions)
    copied = copy_verified_metadata(args.metadata_cache, args.output)
    declaration = {
        "created_utc": base.utc(),
        "status": "declared_before_literature_refinement_and_before_science_flux",
        "wrapper_script_sha256": base.sha(Path(__file__)),
        "inherited_selector_sha256": base.sha(Path(selector.__file__)),
        "base_selector_sha256": base.sha(Path(base.__file__)),
        "previous_completed_manifest_sha256": base.sha(completed_cache),
        "original_manifest_sha256": base.sha(args.previous / "provisional-manifest.json"),
        "literature_exclusions_path": str(args.literature_exclusions),
        "literature_exclusions_sha256": base.sha(args.literature_exclusions),
        "literature_rules": artifact["rules"],
        "literature_sources": artifact["sources"],
        "literature_exclusion_identity_count": len(exclusions),
        "requested_count": args.count,
        "copied_raw_metadata": copied,
        "selection_order": "Original seed-20261009 shuffled pool, unchanged",
        "history_exclusions": sorted(selector.HISTORY_EXCLUSIONS),
        "additional_exclusions": "Append artifact-backed literature reasons to each SIMBAD audit; preserve inherited history reasons separately",
        "science_flux_accessed": False,
        "historical_88_overlap": "unknown",
    }
    declaration_path = args.output / "wrapper-declaration.json"
    if declaration_path.exists():
        old = json.loads(declaration_path.read_text())
        if any(old.get(key) != value for key, value in declaration.items() if key != "created_utc"):
            raise ValueError("V3 wrapper declaration changed; choose a new output directory")
    else:
        base.write_json(declaration_path, declaration)
    inherited_audit = selector.audit_target

    def augmented_audit(identity, metadata, out, definitions, excluded_types):
        audit = inherited_audit(identity, metadata, out, definitions, excluded_types)
        literature_reasons = exclusions.get(identity, [])
        audit["literature_exclusion_reasons"] = literature_reasons
        audit["literature_exclusions_sha256"] = declaration["literature_exclusions_sha256"]
        audit["excluded_reasons"] = sorted(set(audit["excluded_reasons"]) | {"literature:" + reason for reason in literature_reasons})
        base.write_json(out / str(identity) / "audit.json", audit)
        return audit

    selector.audit_target = augmented_audit
    saved_argv = sys.argv
    try:
        sys.argv = [str(Path(selector.__file__)), "--previous", str(args.previous), "--output", str(args.output), "--count", str(args.count)]
        selector.main()
    finally:
        selector.audit_target = inherited_audit
        sys.argv = saved_argv
    manifest_path = args.output / "provisional-manifest.json"
    manifest = json.loads(manifest_path.read_text())
    manifest.update({
        "wrapper_declaration_sha256": base.sha(declaration_path),
        "wrapper_script_sha256": declaration["wrapper_script_sha256"],
        "literature_exclusions_sha256": declaration["literature_exclusions_sha256"],
        "selection_wrapper_complete": True,
        "additional_literature_exclusion_rules": artifact["rules"],
    })
    base.write_json(manifest_path, manifest)
    provenance_path = args.output / "provisional-manifest.provenance.json"
    provenance = json.loads(provenance_path.read_text())
    provenance.update({"sha256": base.sha(manifest_path), "wrapper_declaration_sha256": base.sha(declaration_path)})
    base.write_json(provenance_path, provenance)
    print(f"Literature-refined manifest: {manifest_path} (SHA-256 {base.sha(manifest_path)})", flush=True)


if __name__ == "__main__":
    main()
