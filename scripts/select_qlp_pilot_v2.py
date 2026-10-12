#!/usr/bin/env python3
"""Refine the preserved first metadata selection using SIMBAD and prior-history exclusions."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import re

import select_qlp_pilot as base

HISTORY_EXCLUSIONS = {
    149930150, 149975676, 150064525, 31512254, 167571875, 150250788, 167011368,
    293270903, 167419543, 167422005, 31528780, 150031766, 149121369,
    375036675, 375035828, 41112771, 167571315, 150389853,
}


def simbad_query(query, path):
    url = "https://simbad.cds.unistra.fr/simbad/sim-tap/sync?" + base.urllib.parse.urlencode({
        "request": "doQuery", "lang": "adql", "format": "json", "query": query,
    })
    payload = json.loads(base.retrieve(url, path, timeout=90))
    if "metadata" not in payload or "data" not in payload:
        raise ValueError("SIMBAD result lacks metadata/data")
    names = [column["name"] for column in payload["metadata"]]
    return [dict(zip(names, row)) for row in payload["data"]]


def type_exclusions(definitions):
    selected = {}
    for row in definitions:
        ancestry = set(re.split(r"\s*>\s*", row.get("path") or ""))
        description = (row.get("description") or "").lower()
        if ancestry.intersection({"V*", "EB*", "SB*", "**"}) or re.search(r"variable|binary|eruptive|eclipsing|pulsat", description):
            selected[row["otype"]] = row
    return selected


def audit_target(identity, metadata, out, definitions, excluded_types):
    ra, dec = metadata["ra"], metadata["dec"]
    query = (
        "SELECT b.main_id,b.otype,b.ra,b.dec,b.oid,ids.ids,alltypes.otypes "
        "FROM basic AS b LEFT JOIN ids ON b.oid=ids.oidref "
        "LEFT JOIN alltypes ON b.oid=alltypes.oidref "
        f"WHERE b.oid IN (SELECT oidref FROM ident WHERE id='TIC {identity}') OR "
        f"1=CONTAINS(POINT('ICRS',b.ra,b.dec),CIRCLE('ICRS',{ra},{dec},0.001388888888889))"
    )
    objects = simbad_query(query, out / str(identity) / "objects.json")
    exclusion_reasons = []
    exact = []
    for obj in objects:
        aliases = (obj.get("ids") or "").split("|")
        obj["exact_tic_alias"] = f"TIC {identity}" in aliases
        obj["separation_arcsec"] = base.angular_separation(ra, dec, obj["ra"], obj["dec"])
        if obj["exact_tic_alias"]:
            exact.append(obj)
        types = set((obj.get("otypes") or obj["otype"]).split("|")) | {obj["otype"]}
        obj["type_definitions"] = [definitions[code] for code in sorted(types) if code in definitions]
        obj["excluded_types"] = [code for code in sorted(types) if code in excluded_types]
        refs_query = (
            "SELECT r.bibcode,r.title FROM has_ref AS h JOIN ref AS r ON h.oidbibref=r.oidbib "
            f"WHERE h.oidref={int(obj['oid'])}"
        )
        try:
            obj["references"] = simbad_query(refs_query, out / str(identity) / f"references-{obj['oid']}.json")
            obj["reference_query_status"] = "complete"
        except Exception as error:
            obj["references"] = []
            obj["reference_query_status"] = "failed"
            obj["reference_error"] = f"{type(error).__name__}: {error}"
    matched = exact if len(exact) == 1 else objects
    if len(matched) > 1:
        exclusion_reasons.append("ambiguous_simbad_objects")
    for obj in matched:
        exclusion_reasons.extend("simbad_type:" + code for code in obj["excluded_types"])
    report = {
        "tic": identity, "query_status": "complete", "objects": objects,
        "identity_match": "exact_tic_alias" if len(exact) == 1 else ("position_only" if len(objects) == 1 else "ambiguous_or_absent"),
        "excluded_reasons": sorted(set(exclusion_reasons)),
        "limitations": "SIMBAD types and associated references do not provide an exhaustive prior-literature or variability audit; parent-sample membership is not a prior signal report.",
    }
    base.write_json(out / str(identity) / "audit.json", report)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--previous", type=Path, default=Path("data/qlp-pilot"))
    parser.add_argument("--output", type=Path, default=Path("data/qlp-pilot-v2"))
    parser.add_argument("--count", type=int, default=24)
    args = parser.parse_args()
    old, out = args.previous, args.output
    out.mkdir(parents=True, exist_ok=True)
    original = json.loads((old / "provisional-manifest.json").read_text())
    previous_targets = {row["tic"]: row for row in original["targets"]}
    pool = json.loads((old / "shuffled-pool.json").read_text())
    rows = {row["tic"]: row["metadata"] for row in json.loads((old / "catalog-ledger.json").read_text())}
    declaration = {
        "created_utc": base.utc(), "requested_count": args.count,
        "script_sha256": base.sha(Path(__file__)), "base_script_sha256": base.sha(Path(base.__file__)),
        "previous_manifest_sha256": base.sha(old / "provisional-manifest.json"),
        "previous_pool_sha256": base.sha(old / "shuffled-pool.json"),
        "history_exclusions": sorted(HISTORY_EXCLUSIONS),
        "history_exclusion_meaning": "Previously mentioned targets excluded conservatively; does not assert current catalog classifications or reconstruct full historical cohort.",
        "simbad_exclusion": "Exact TIC alias or 5 arcsec coordinate match with variable, eclipsing/spectroscopic/general binary, eruptive or pulsating type ancestry/description; ambiguous positional matches and failed object queries excluded.",
        "seed_order": "Reuse original seed-20261009 shuffled pool without reshuffling after further exclusions",
        "discovery_sectors": [1, 2], "holdout_policy": "Unchanged from the first provisional selection",
        "status": "provisional_before_control_review", "flux_accessed": False,
    }
    path = out / "selection-declaration.json"
    if path.exists():
        previous_declaration = json.loads(path.read_text())
        if any(previous_declaration[key] != declaration[key] for key in declaration if key != "created_utc"):
            raise ValueError("Selection declaration changed; use a new output directory")
    else:
        base.write_json(path, declaration)
    typedefs = simbad_query("SELECT * FROM otypedef", out / "simbad-type-definitions.json")
    definitions = {row["otype"]: row for row in typedefs}
    excluded_types = type_exclusions(typedefs)
    base.write_json(out / "excluded-type-definitions.json", excluded_types)
    # Audit all first-version targets, including history-excluded controls, before replacements.
    audits = {}
    for identity, row in previous_targets.items():
        try:
            audits[identity] = audit_target(identity, row["metadata"], out / "simbad", definitions, excluded_types)
        except Exception as error:
            audits[identity] = {"tic": identity, "query_status": "failed", "error": f"{type(error).__name__}: {error}", "excluded_reasons": ["simbad_query_failed"]}
        print(f"Initial SIMBAD TIC={identity} status={audits[identity]['query_status']} reasons={audits[identity]['excluded_reasons']}", flush=True)
    selected, ledger = [], []
    for rank, identity in enumerate(pool, 1):
        row = rows[identity]
        reasons = ["history_previously_mentioned"] if identity in HISTORY_EXCLUSIONS else []
        observations_directory = old / "observations" / str(identity)
        if not observations_directory.exists():
            observations_directory = out / "observations" / str(identity)
        query = {"service": "Mast.Caom.Filtered", "params": {
            "columns": "obsid,obs_id,target_name,sequence_number,t_min,t_max,provenance_name,dataURL",
            "filters": [{"paramName": "provenance_name", "values": ["QLP"]}, {"paramName": "target_name", "values": [str(identity)]}]}}
        usable = None
        entry = {"tic": identity, "rank": rank, "history_excluded": bool(reasons)}
        if not reasons:
            try:
                observations = base.mast_pages(query, observations_directory, page_size=1000)
                usable, metadata_reason = base.eligible_observations(observations)
                entry["observation_query_status"] = "complete"
                if metadata_reason:
                    reasons.append(metadata_reason)
            except Exception as error:
                reasons.append("observation_query_failed")
                entry["observation_query_status"] = "failed"
                entry["error"] = f"{type(error).__name__}: {error}"
        if usable and identity not in audits:
            try:
                audits[identity] = audit_target(identity, row, out / "simbad", definitions, excluded_types)
            except Exception as error:
                audits[identity] = {"tic": identity, "query_status": "failed", "error": f"{type(error).__name__}: {error}", "excluded_reasons": ["simbad_query_failed"]}
        if identity in audits:
            reasons.extend(audits[identity]["excluded_reasons"])
            entry["simbad_query_status"] = audits[identity]["query_status"]
        if usable and not reasons:
            selected.append({"tic": identity, "selection_rank": rank, "metadata": row, **usable,
                             "simbad_audit": str(out / "simbad" / str(identity) / "audit.json"),
                             "simbad_identity_match": audits[identity]["identity_match"]})
        entry["selected"] = usable is not None and not reasons
        entry["excluded_reasons"] = sorted(set(reasons))
        ledger.append(entry)
        base.write_json(out / "selection-ledger.json", ledger)
        print(f"rank={rank} TIC={identity} selected={len(selected)}/{args.count} reasons={reasons}", flush=True)
        if len(selected) >= args.count:
            break
    base.write_json(out / "simbad-audits.json", audits)
    manifest = {
        "created_utc": base.utc(), "status": "provisional_metadata_selection_pending_control_review",
        "previous_manifest_sha256": declaration["previous_manifest_sha256"],
        "selection_declaration_sha256": base.sha(path), "selection_script_sha256": base.sha(Path(__file__)),
        "complete_tic_rows": original["complete_tic_rows"], "catalog_rows": original["catalog_rows"],
        "post_first_exclusion_pool": len(pool), "seed_ranks_examined": len(ledger),
        "selected_count": len(selected), "requested_count": args.count,
        "historical_88_overlap": "unknown", "old_cohort_overlap_unknown": True, "flux_accessed": False,
        "targets": selected,
    }
    base.write_json(out / "provisional-manifest.json", manifest)
    base.write_json(out / "provisional-manifest.provenance.json", {"created_utc": base.utc(), "sha256": base.sha(out / "provisional-manifest.json"), "declaration_sha256": base.sha(path)})
    print(f"Selected {len(selected)} targets: {out / 'provisional-manifest.json'}", flush=True)


if __name__ == "__main__":
    main()
