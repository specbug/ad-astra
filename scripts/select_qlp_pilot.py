#!/usr/bin/env python3
"""Freeze a metadata-only QLP pilot, preserving raw catalog and query provenance."""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import io
import json
import math
from pathlib import Path
import random
import re
import urllib.parse
import urllib.request

MAST = "https://mast.stsci.edu/api/v0/invoke"
EXPLICIT_EXCLUSIONS = {150250788, 167011368, 370612437, 252721481, 350883383}


def utc():
    return datetime.now(timezone.utc).isoformat()


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def retrieve(url, path, timeout=180):
    """Reuse an intact, recorded snapshot; never silently replace its provenance."""
    sidecar = path.with_name(path.name + ".provenance.json")
    if path.exists() and sidecar.exists():
        previous = json.loads(sidecar.read_text())
        if previous["url"] != url or previous["sha256"] != sha(path):
            raise ValueError(f"Cache provenance mismatch: {path}")
        return path.read_bytes()
    request = urllib.request.Request(url, headers={"User-Agent": "ad-astra-research/0.1"})
    with urllib.request.urlopen(request, timeout=timeout) as response:
        payload = response.read()
        provenance = {
            "retrieved_utc": utc(), "url": url, "http_status": response.status,
            "sha256": hashlib.sha256(payload).hexdigest(), "bytes": len(payload),
            "content_type": response.headers.get("Content-Type"),
        }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(payload)
    write_json(sidecar, provenance)
    return payload


def mast_pages(query, directory, page_size=10000, max_rows=100000):
    """Save every raw page and fail closed on partial/inconsistent pagination."""
    rows = []
    seen_pages = set()
    page = 1
    total = None
    while True:
        q = dict(query, format="json", pagesize=page_size, page=page)
        write_json(directory / f"request-{page:03d}.json", q)
        url = MAST + "?" + urllib.parse.urlencode({"request": json.dumps(q)})
        response = json.loads(retrieve(url, directory / f"response-{page:03d}.json", timeout=240))
        if response.get("status") != "COMPLETE":
            raise RuntimeError(f"Incomplete MAST response: {response.get('status')} {response.get('msg')}")
        values = response.get("data", [])
        digest = hashlib.sha256(json.dumps(values, sort_keys=True).encode()).hexdigest()
        if values and digest in seen_pages:
            raise ValueError("MAST repeated a page; completeness cannot be established")
        seen_pages.add(digest)
        rows.extend(values)
        if len(rows) > max_rows:
            raise ValueError("Metadata cap exceeded; refusing a truncated population")
        paging = response.get("paging", {})
        reported = paging.get("rowsFiltered")
        if reported is not None:
            if total is not None and reported != total:
                raise ValueError("MAST filtered row count changed during pagination")
            total = reported
        pages = paging.get("pagesFiltered")
        if pages is not None and page >= pages:
            break
        if total is not None and len(rows) >= total:
            break
        if len(values) < page_size:
            break
        page += 1
    if total is not None and len(rows) != total:
        raise ValueError(f"MAST expected {total} filtered rows, retrieved {len(rows)}")
    write_json(directory / "complete.json", {"completed_utc": utc(), "rows": len(rows), "pages": page})
    return rows


def tic_number(value):
    if value is None or str(value).strip() == "":
        return None
    match = re.fullmatch(r"(?:TIC\s*)?(\d+)(?:\.0)?", str(value).strip())
    return int(match.group(1)) if match else None


def degrees(value, hours=False):
    try:
        if value is None or str(value).strip() == "":
            return None
        text = str(value).strip()
        if ":" not in text and " " not in text:
            return float(text)
        parts = [float(x) for x in re.split(r"[:\s]+", text)]
        sign = -1 if text.startswith("-") else 1
        return sign * (abs(parts[0]) + parts[1] / 60 + parts[2] / 3600) * (15 if hours else 1)
    except (ValueError, IndexError):
        return None


def angular_separation(ra, dec, other_ra, other_dec):
    radians = math.pi / 180
    sd = math.sin((dec - other_dec) * radians / 2)
    sr = math.sin((ra - other_ra) * radians / 2)
    h = sd * sd + math.cos(dec * radians) * math.cos(other_dec * radians) * sr * sr
    return 2 * math.asin(math.sqrt(min(1, max(0, h)))) / radians * 3600


def catalog_exclusions(directory):
    queries = {
        "toi": "select tid,toi,tfopwg_disp,ra,dec from toi",
        "confirmed": "select pl_name,tic_id,ra,dec from pscomppars",
    }
    rows = {}
    for name, query in queries.items():
        url = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync?" + urllib.parse.urlencode({"query": query, "format": "json"})
        rows[name] = json.loads(retrieve(url, directory / f"{name}.json"))
    payload = retrieve("https://exofop.ipac.caltech.edu/tess/download_ctoi.php?sort=ctoi&output=csv", directory / "ctoi.csv")
    rows["ctoi"] = list(csv.DictReader(io.StringIO(payload.decode("utf-8-sig"))))
    ids = {}
    coordinates = []
    for source, values in rows.items():
        for row in values:
            identity = tic_number(row.get({"toi": "tid", "confirmed": "tic_id", "ctoi": "TIC ID"}[source]))
            if identity:
                ids.setdefault(identity, set()).add(source)
            ra = degrees(row.get("ra", row.get("RA")), hours=source == "ctoi")
            dec = degrees(row.get("dec", row.get("Dec")))
            if ra is not None and dec is not None and math.isfinite(ra) and math.isfinite(dec):
                coordinates.append((ra, dec, source, identity))
    counts = {name: len(values) for name, values in rows.items()}
    write_json(directory / "summary.json", {"created_utc": utc(), "rows": counts, "unique_tics": len(ids), "coordinate_rows": len(coordinates)})
    return ids, coordinates, counts


def eligible_observations(observations):
    """Use only legacy QLP schema sectors; reserve the last separated season."""
    by_sector = {}
    for row in observations:
        sector = int(row["sequence_number"])
        uri = row.get("dataURL") or ""
        if not (1 <= sector <= 73 or 80 <= sector <= 93):
            continue
        if row.get("provenance_name") != "QLP" or not uri.endswith("_llc.fits"):
            continue
        if row.get("t_min") is None or row.get("t_max") is None:
            continue
        duration = row["t_max"] - row["t_min"]
        if not 0 < duration <= 40:
            continue
        normalized = {"sector": sector, "uri": uri, "obsid": row["obsid"], "obs_id": row["obs_id"], "t_min": row["t_min"], "t_max": row["t_max"]}
        # Deterministic URI ordering selects the latest named processing version.
        if sector not in by_sector or uri > by_sector[sector]["uri"]:
            by_sector[sector] = normalized
    if not {1, 2}.issubset(by_sector):
        return None, "missing_sector_1_or_2"
    later = sorted((row for sector, row in by_sector.items() if sector >= 27), key=lambda row: row["t_min"])
    if len(later) < 2:
        return None, "fewer_than_two_later_sectors"
    seasons = [[]]
    for row in later:
        if seasons[-1] and row["t_min"] - seasons[-1][-1]["t_max"] >= 180:
            seasons.append([])
        seasons[-1].append(row)
    holdout = seasons[-1]
    if len(holdout) < 2:
        return None, "latest_eligible_season_has_fewer_than_two_sectors"
    discovery = [by_sector[1], by_sector[2]]
    if holdout[0]["t_min"] - discovery[-1]["t_max"] < 180:
        return None, "insufficient_discovery_holdout_gap"
    return {"discovery": discovery, "holdout": holdout, "available_safe_sectors": sorted(by_sector)}, None


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("data/qlp-pilot"))
    parser.add_argument("--count", type=int, default=24)
    parser.add_argument("--seed", type=int, default=20261009)
    args = parser.parse_args()
    out = args.output
    out.mkdir(parents=True, exist_ok=True)
    declaration = {
        "created_utc": utc(), "script_sha256": sha(Path(__file__)), "seed": args.seed,
        "requested_stars": args.count, "center_deg": [90, -66.5], "radius_deg": 5,
        "cuts": {"Tmag": [8, 12], "rad": [0.3, 1.5], "contratio": "0 <= value < 0.1", "abs_gallat": ">15"},
        "discovery_sectors": [1, 2], "holdout_policy": "latest available QLP season >=27 with >=2 sectors; season gap >=180 days; all sectors 74–79 and >=94 excluded",
        "coordinate_catalog_exclusion_arcsec": 5,
        "explicit_exclusions": sorted(EXPLICIT_EXCLUSIONS),
        "historical_88_overlap": "unknown: historical manifest unavailable after runtime reset",
        "selection_basis": "metadata only; no flux downloaded by this script",
        "prior_pipeline_audit": "No SPOC DV or literature audit in this selection step; archive catalog exclusions do not establish novelty.",
    }
    declaration_path = out / "selection-declaration.json"
    if declaration_path.exists():
        old = json.loads(declaration_path.read_text())
        if any(old[key] != declaration[key] for key in declaration if key != "created_utc"):
            raise ValueError("Selection declaration changed: choose a new output directory")
    else:
        write_json(declaration_path, declaration)
    ids, coordinates, catalog_counts = catalog_exclusions(out / "catalogs")
    print("Catalogs:", catalog_counts, flush=True)
    tic_request = {
        "service": "Mast.Catalogs.Filtered.Tic.Position",
        "params": {"columns": "ID,ra,dec,Tmag,rad,contratio,gallat", "ra": 90, "dec": -66.5, "radius": 5,
                   "filters": [{"paramName": "Tmag", "values": [{"min": 8, "max": 12}]},
                               {"paramName": "rad", "values": [{"min": 0.3, "max": 1.5}]},
                               {"paramName": "contratio", "values": [{"min": 0, "max": 0.1}]}]},
    }
    tic_rows = mast_pages(tic_request, out / "tic-query")
    if len({str(row["ID"]) for row in tic_rows}) != len(tic_rows):
        raise ValueError("Duplicate TIC identities in cone query")
    ledger = []
    candidates = []
    for row in tic_rows:
        identity = int(row["ID"])
        reasons = []
        if row.get("Tmag") is None or not 8 <= float(row["Tmag"]) <= 12:
            reasons.append("magnitude_cut")
        if row.get("rad") is None or not 0.3 <= float(row["rad"]) <= 1.5:
            reasons.append("stellar_radius_cut")
        if angular_separation(float(row["ra"]), float(row["dec"]), 90, -66.5) > 5 * 3600:
            reasons.append("cone_cut")
        if row.get("gallat") is None or abs(float(row["gallat"])) <= 15:
            reasons.append("galactic_latitude_cut")
        if row.get("contratio") is None or not 0 <= float(row["contratio"]) < 0.1:
            reasons.append("contamination_cut")
        if identity in EXPLICIT_EXCLUSIONS:
            reasons.append("explicit_previous_lead_or_reference")
        reasons.extend("catalog_tic:" + source for source in sorted(ids.get(identity, [])))
        # Include positional exclusions for catalog entries without a usable TIC ID.
        for ra, dec, source, catalog_tic in coordinates:
            if abs(float(row["dec"]) - dec) * 3600 > 5:
                continue
            if angular_separation(float(row["ra"]), float(row["dec"]), ra, dec) <= 5:
                label = "catalog_position:" + source
                if label not in reasons:
                    reasons.append(label)
        ledger.append({"tic": identity, "catalog_exclusions": reasons, "metadata": row})
        if not reasons:
            candidates.append(row)
    candidates.sort(key=lambda row: int(row["ID"]))
    random.Random(args.seed).shuffle(candidates)
    write_json(out / "catalog-ledger.json", ledger)
    write_json(out / "shuffled-pool.json", [int(row["ID"]) for row in candidates])
    print(f"TIC rows {len(tic_rows)}, remaining {len(candidates)}", flush=True)
    selected = []
    observation_ledger = []
    for rank, row in enumerate(candidates, 1):
        identity = int(row["ID"])
        query = {"service": "Mast.Caom.Filtered", "params": {
            "columns": "obsid,obs_id,target_name,sequence_number,t_min,t_max,provenance_name,dataURL",
            "filters": [{"paramName": "provenance_name", "values": ["QLP"]}, {"paramName": "target_name", "values": [str(identity)]}]}}
        try:
            observations = mast_pages(query, out / "observations" / str(identity), page_size=1000)
            usable, reason = eligible_observations(observations)
            entry = {"tic": identity, "rank": rank, "query_status": "complete", "observations": len(observations), "excluded_reason": reason}
            if usable:
                selected.append({"tic": identity, "selection_rank": rank, "metadata": row, **usable})
                entry["selected"] = True
            else:
                entry["selected"] = False
        except Exception as error:
            entry = {"tic": identity, "rank": rank, "query_status": "failed", "error": f"{type(error).__name__}: {error}", "selected": False}
        observation_ledger.append(entry)
        write_json(out / "observation-ledger.json", observation_ledger)
        print(f"rank={rank} TIC={identity} selected={len(selected)}/{args.count} status={entry.get('excluded_reason',entry['query_status'])}", flush=True)
        if len(selected) >= args.count:
            break
    manifest = {
        "created_utc": utc(), "status": "provisional_metadata_selection_pending_control_review", "selection_declaration_sha256": sha(declaration_path),
        "selection_script_sha256": sha(Path(__file__)), "complete_tic_rows": len(tic_rows),
        "catalog_rows": catalog_counts, "post_exclusion_pool": len(candidates),
        "observation_queries": len(observation_ledger), "query_failures": sum(row["query_status"] == "failed" for row in observation_ledger),
        "selected_count": len(selected), "requested_count": args.count,
        "historical_88_overlap": "unknown", "old_cohort_overlap_unknown": True, "flux_accessed": False,
        "targets": selected,
    }
    write_json(out / "provisional-manifest.json", manifest)
    write_json(out / "provisional-manifest.provenance.json", {"created_utc": utc(), "sha256": sha(out / "provisional-manifest.json"), "declaration_sha256": sha(declaration_path)})
    print(f"Selected {len(selected)} targets provisionally: {out / 'provisional-manifest.json'}", flush=True)


if __name__ == "__main__":
    main()
