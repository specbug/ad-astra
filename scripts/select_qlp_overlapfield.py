#!/usr/bin/env python3
"""Resume-safe, metadata-only selection in the fixed overlapping sectors 1/2 field.

This script does not retrieve light curves. Query failures remain exclusions on
resume; a later retry or change in policy requires a separate output directory.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
import io
import json
import math
from pathlib import Path
import random
import urllib.parse

import select_qlp_pilot as base
import select_qlp_pilot_v2 as simbad
from select_qlp_pilot_v3 import load_exclusions


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + ".part")
    temp.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    temp.replace(path)


def log(out, message):
    line = f"{base.utc()} {message}"
    print(line, flush=True)
    with (out / "progress.log").open("a") as stream:
        stream.write(line + "\n")


def cached_catalogs(directory):
    """Read the exact saved NASA/CTOI snapshots without refreshing them."""
    rows, records = {}, []
    for name, suffix in (("toi", "json"), ("confirmed", "json"), ("ctoi", "csv")):
        path = directory / f"{name}.{suffix}"
        provenance = json.loads(path.with_name(path.name + ".provenance.json").read_text())
        if provenance["sha256"] != base.sha(path):
            raise ValueError(f"Catalog snapshot hash mismatch: {path}")
        rows[name] = (list(csv.DictReader(io.StringIO(path.read_text(encoding="utf-8-sig"))))
                      if suffix == "csv" else json.loads(path.read_text()))
        records.append({"path": str(path), "sha256": base.sha(path), "provenance": provenance})
    identities, coordinates = {}, []
    for source, values in rows.items():
        for row in values:
            tic = base.tic_number(row.get({"toi": "tid", "confirmed": "tic_id", "ctoi": "TIC ID"}[source]))
            if tic:
                identities.setdefault(tic, set()).add(source)
            ra = base.degrees(row.get("ra", row.get("RA")), hours=source == "ctoi")
            dec = base.degrees(row.get("dec", row.get("Dec")))
            if ra is not None and dec is not None and math.isfinite(ra) and math.isfinite(dec):
                coordinates.append((ra, dec, source))
    return identities, coordinates, records, {name: len(values) for name, values in rows.items()}


def gaia_query(query, path):
    url = "https://gea.esac.esa.int/tap-server/tap/sync?" + urllib.parse.urlencode({
        "REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "json", "MAXREC": 100000, "QUERY": query,
    })
    result = json.loads(base.retrieve(url, path, timeout=120))
    if "metadata" not in result or "data" not in result:
        raise ValueError("Gaia TAP response lacks metadata/data")
    names = [column["name"] for column in result["metadata"]]
    return [dict(zip(names, row)) for row in result["data"]]


def audit_gaia(tic, meta, out):
    condition = ("1=CONTAINS(POINT('ICRS',ra,dec),"
                 f"CIRCLE('ICRS',{meta['ra']},{meta['dec']},0.025))")
    directory = out / "gaia" / str(tic)
    rows = gaia_query("SELECT source_id,ra,dec,phot_g_mean_mag,phot_variable_flag "
                      "FROM gaiadr3.gaia_source WHERE " + condition, directory / "cone.json")
    counts = gaia_query("SELECT COUNT(*) AS n FROM gaiadr3.gaia_source WHERE " + condition,
                        directory / "count.json")
    if len(counts) != 1 or counts[0]["n"] != len(rows):
        raise ValueError("Gaia cone independent count does not match returned rows")
    for row in rows:
        row["separation_arcsec"] = base.angular_separation(meta["ra"], meta["dec"], row["ra"], row["dec"])
    variable = [row for row in rows if row["phot_variable_flag"] == "VARIABLE"]
    audit = {
        "tic": tic, "query_status": "complete", "radius_arcsec": 90,
        "complete_cone_rows": len(rows), "objects": rows, "variable_neighbors": variable,
        "excluded_reasons": ["gaia_dr3_variable_within_90_arcsec"] if variable else [],
        "target_coverage": "position_match_within_5_arcsec" if any(row["separation_arcsec"] <= 5 for row in rows) else "unchecked_no_close_gaia_source",
        "unchecked_variability_flag_rows": sum(row["phot_variable_flag"] not in {"VARIABLE", "CONSTANT"} for row in rows),
        "limitations": "Gaia positions are epoch 2016 and are not propagated. CONSTANT/NOT_AVAILABLE or no catalog entry does not demonstrate absence of variability; this conservative neighborhood exclusion biases the sample.",
    }
    write(directory / "audit.json", audit)
    return audit


def files_for_target(out, tic):
    files = []
    for subdir in ("observations", "simbad", "gaia"):
        files.extend(path for path in (out / subdir / str(tic)).rglob("*") if path.is_file())
    return {str(path.relative_to(out)): base.sha(path) for path in sorted(files)}


def candidate(rank, meta, out, declaration_sha, definitions, excluded_types):
    tic = int(meta["ID"])
    path = out / "target-state" / f"{tic}.json"
    if path.exists():
        state = json.loads(path.read_text())
        if state["declaration_sha256"] != declaration_sha or state["rank"] != rank or state["metadata"] != meta:
            raise ValueError(f"Saved target state differs from declared selection: {tic}")
        for relative, digest in state["metadata_files_sha256"].items():
            if base.sha(out / relative) != digest:
                raise ValueError(f"Saved metadata hash differs: {relative}")
        return state
    state = {"tic": tic, "rank": rank, "metadata": meta, "declaration_sha256": declaration_sha,
             "started_utc": base.utc(), "excluded_reasons": [], "eligible": False, "target": None}
    query = {"service": "Mast.Caom.Filtered", "params": {
        "columns": "obsid,obs_id,target_name,sequence_number,t_min,t_max,provenance_name,dataURL",
        "filters": [{"paramName": "provenance_name", "values": ["QLP"]},
                    {"paramName": "target_name", "values": [str(tic)]}]}}
    stage = "observation"
    try:
        observations = base.mast_pages(query, out / "observations" / str(tic), page_size=1000)
        usable, reason = base.eligible_observations(observations)
        state["observation_query_status"] = "complete"
        state["qlp_metadata_rows"] = len(observations)
        state["available_qlp_sectors"] = sorted({int(row["sequence_number"]) for row in observations})
        if reason:
            state["excluded_reasons"].append(reason)
        if usable:
            # Preserve observed provenance explicitly for the future frozen runner.
            originals = {row["dataURL"]: row for row in observations if row.get("dataURL")}
            for product in usable["discovery"] + usable["holdout"]:
                product["provenance_name"] = originals[product["uri"]]["provenance_name"]
            stage = "simbad"
            audit = simbad.audit_target(tic, meta, out / "simbad", definitions, excluded_types)
            state["simbad_query_status"] = audit["query_status"]
            state["simbad_identity_match"] = audit["identity_match"]
            state["excluded_reasons"].extend(audit["excluded_reasons"])
            if not state["excluded_reasons"]:
                stage = "gaia"
                gaia = audit_gaia(tic, meta, out)
                state["gaia_query_status"] = gaia["query_status"]
                state["gaia_target_coverage"] = gaia["target_coverage"]
                state["excluded_reasons"].extend(gaia["excluded_reasons"])
                if not state["excluded_reasons"]:
                    state["eligible"] = True
                    state["target"] = {"tic": tic, "selection_rank": rank, "metadata": meta, **usable,
                        "simbad_audit": str(out / "simbad" / str(tic) / "audit.json"),
                        "simbad_identity_match": audit["identity_match"],
                        "gaia_audit": str(out / "gaia" / str(tic) / "audit.json"),
                        "gaia_target_coverage": gaia["target_coverage"]}
    except Exception as error:
        state[stage + "_query_status"] = "failed"
        state["error"] = f"{type(error).__name__}: {error}"
        state["excluded_reasons"].append(stage + "_query_failed")
    state["excluded_reasons"] = sorted(set(state["excluded_reasons"]))
    state["completed_utc"] = base.utc()
    state["metadata_files_sha256"] = files_for_target(out, tic)
    write(path, state)
    return state


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("data/qlp-overlapfield"))
    parser.add_argument("--workers", type=int, default=4, choices=range(1, 5))
    args = parser.parse_args()
    out = args.output
    out.mkdir(parents=True, exist_ok=True)
    identities, coordinates, catalog_sources, catalog_counts = cached_catalogs(Path("data/qlp-pilot/catalogs"))
    literature_path = Path("data/qlp-prior-review/literature-exclusions.json")
    artifact, literature = load_exclusions(literature_path)
    known_manifest = Path("runs/qlp-pilot/frozen-manifest.json")
    history = simbad.HISTORY_EXCLUSIONS | base.EXPLICIT_EXCLUSIONS
    history |= {row["tic"] for row in json.loads(known_manifest.read_text())["targets"]}
    declaration = {
        "created_utc": base.utc(), "status": "declared_before_metadata_queries_and_science_flux",
        "sources_sha256": {str(path): base.sha(path) for path in (Path(__file__), Path(base.__file__), Path(simbad.__file__), Path("scripts/select_qlp_pilot_v3.py"))},
        "center_deg": [7.1992341073908435, -59.95377502446895], "radius_deg": 3, "seed": 20261011, "requested_count": 24,
        "cuts": {"Tmag": [8, 12], "rad": [0.3, 1.5], "contratio": "0 <= value < 0.1", "abs_gallat": ">15"},
        "catalog_snapshot_sources": catalog_sources, "coordinate_catalog_exclusion_arcsec": 5,
        "literature_exclusions_path": str(literature_path), "literature_exclusions_sha256": base.sha(literature_path),
        "literature_rules": artifact["rules"], "history_exclusions": sorted(history),
        "known_science_manifest": {"path": str(known_manifest), "sha256": base.sha(known_manifest)},
        "discovery_sectors": [1, 2],
        "holdout_policy": "Latest legacy QLP season >=27 with >=2 sectors and >=180-day discovery gap; seasons split at >=180-day gaps; exclude sectors74–79 and >=94",
        "simbad_policy": "Reuse v2 exact TIC alias or 5 arcsec object and variable/binary type audit; ambiguous matches and failed queries exclude; absent catalog entry remains unchecked",
        "gaia_policy": "Exclude any Gaia DR3 phot_variable_flag VARIABLE within90 arcsec; query all cone objects and independent count; failed/incomplete queries exclude; missing/NOT_AVAILABLE coverage remains unchecked",
        "selection_order": "Sort TIC IDs then fixed-seed shuffle; first24 eligible; batches up to4 queried and ledgered even if quota reached midbatch",
        "failure_policy": "Persist completed target failures and exclude; never retry them on resume; cache must match recorded SHA256; interrupted unfinished queries can resume from verified raw responses",
        "sample_limitations": "The neighborhood-variable cut is conservative and biased; this is exploratory target selection, not a representative population or a proof of novelty.",
        "historical_88_overlap": "unknown: original manifest lost; literal sectors1/2 safeguard retained",
        "science_flux_accessed": False,
        "previous_coverage_pilot": {"path": "data/qlp-cleanfield/coverage-stop.json", "sha256": base.sha(Path("data/qlp-cleanfield/coverage-stop.json"))},
        "field_change_rationale": "Previous RA0/Dec-80/r5 field lacks sector2 in all348 completed metadata queries (including first8); new center derives from ecliptic longitude330 latitude-55, geometrically midway the S1/S2 boresight longitudes315.8/343.36. New radius3 chosen before querying this field. This is a separate metadata experiment, with no flux used for field choice.",
        "early_coverage_stop": "Stop for investigation if all first8 seed ranks lack literal sector1 or2; do not substitute different sectors or alter field after observing light curves.",
    }
    declaration_path = out / "selection-declaration.json"
    if declaration_path.exists():
        previous = json.loads(declaration_path.read_text())
        if any(previous.get(key) != value for key, value in declaration.items() if key != "created_utc"):
            raise ValueError("Declaration changed; select a separate output directory")
    else:
        write(declaration_path, declaration)
    declaration_sha = base.sha(declaration_path)
    log(out, "Declaration verified; obtaining fixed-field TIC metadata")
    tic_request = {"service": "Mast.Catalogs.Filtered.Tic.Position", "params": {
        "columns": "ID,ra,dec,Tmag,rad,contratio,gallat", "ra": 7.1992341073908435, "dec": -59.95377502446895, "radius": 3,
        "filters": [{"paramName": "Tmag", "values": [{"min": 8, "max": 12}]},
                    {"paramName": "rad", "values": [{"min": 0.3, "max": 1.5}]},
                    {"paramName": "contratio", "values": [{"min": 0, "max": 0.1}]}]}}
    rows = base.mast_pages(tic_request, out / "tic-query")
    if len({int(row["ID"]) for row in rows}) != len(rows):
        raise ValueError("Duplicate TIC rows")
    ledger, pool = [], []
    for row in rows:
        tic, reasons = int(row["ID"]), []
        for field, low, high in (("Tmag", 8, 12), ("rad", 0.3, 1.5)):
            if row.get(field) is None or not low <= row[field] <= high:
                reasons.append(field + "_cut")
        if row.get("contratio") is None or not 0 <= row["contratio"] < 0.1:
            reasons.append("contamination_cut")
        if row.get("gallat") is None or abs(row["gallat"]) <= 15:
            reasons.append("galactic_latitude_cut")
        if base.angular_separation(row["ra"], row["dec"], 7.1992341073908435, -59.95377502446895) > 10800:
            reasons.append("cone_cut")
        if tic in history:
            reasons.append("history_previously_examined_or_mentioned")
        reasons.extend("literature:" + reason for reason in literature.get(tic, []))
        reasons.extend("catalog_tic:" + source for source in identities.get(tic, []))
        reasons.extend("catalog_position:" + source for ra, dec, source in coordinates
                       if abs(row["dec"] - dec) * 3600 <= 5 and base.angular_separation(row["ra"], row["dec"], ra, dec) <= 5)
        ledger.append({"tic": tic, "metadata": row, "excluded_reasons": sorted(set(reasons))})
        if not reasons:
            pool.append(row)
    pool.sort(key=lambda row: int(row["ID"]))
    random.Random(20261011).shuffle(pool)
    write(out / "catalog-ledger.json", ledger)
    write(out / "shuffled-pool.json", [int(row["ID"]) for row in pool])
    log(out, f"Complete TIC field rows={len(rows)} post-exclusion pool={len(pool)}; S1/2 availability will be recorded for every examined target")
    definitions_list = simbad.simbad_query("SELECT * FROM otypedef", out / "simbad-type-definitions.json")
    definitions = {row["otype"]: row for row in definitions_list}
    excluded_types = simbad.type_exclusions(definitions_list)
    write(out / "excluded-type-definitions.json", excluded_types)
    states, selected = [], []
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        for start in range(0, len(pool), args.workers):
            batch = list(enumerate(pool[start:start + args.workers], start + 1))
            futures = [executor.submit(candidate, rank, row, out, declaration_sha, definitions, excluded_types) for rank, row in batch]
            for future in futures:
                state = future.result()
                state["selected"] = state["eligible"] and len(selected) < 24
                if state["selected"]:
                    selected.append(state["target"])
                states.append(state)
                write(out / "selection-ledger.json", states)
                log(out, f"rank={state['rank']} TIC={state['tic']} selected={len(selected)}/24 available={state.get('available_qlp_sectors', [])} exclusions={state['excluded_reasons']}")
            manifest = {
                "created_utc": base.utc(), "status": "provisional_metadata_selection_pending_control_review",
                "selection_complete": len(selected) == 24 or start + args.workers >= len(pool),
                "selection_declaration_sha256": declaration_sha, "selection_script_sha256": base.sha(Path(__file__)),
                "complete_tic_rows": len(rows), "catalog_rows": catalog_counts, "post_exclusion_pool": len(pool),
                "seed_ranks_examined": len(states), "selected_count": len(selected), "requested_count": 24,
                "missing_sector_1_or_2": sum("missing_sector_1_or_2" in row["excluded_reasons"] for row in states),
                "query_failures": sum(any(reason.endswith("query_failed") for reason in row["excluded_reasons"]) for row in states),
                "historical_88_overlap": "unknown", "old_cohort_overlap_unknown": True, "flux_accessed": False,
                "targets": selected,
            }
            write(out / "provisional-manifest.json", manifest)
            write(out / "provisional-manifest.provenance.json", {"created_utc": base.utc(), "sha256": base.sha(out / "provisional-manifest.json"), "declaration_sha256": declaration_sha})
            if len(states) >= 8 and all("missing_sector_1_or_2" in row["excluded_reasons"] for row in states[:8]):
                write(out / "coverage-stop.json", {"created_utc": base.utc(), "reason": "All first8 seed ranks lack literal S1/2; investigate geometry before further metadata or flux", "seed_ranks_examined": len(states), "flux_accessed": False})
                break
            if len(selected) >= 24:
                break
    log(out, f"Complete: selected={len(selected)} examined={len(states)} fixed-field pool={len(pool)}; flux_accessed=False")


if __name__ == "__main__":
    main()
