#!/usr/bin/env python3
"""Read-only prior-catalog queries; downloads no light curves or science images."""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import hashlib
import io
import json
from pathlib import Path
import urllib.parse
import urllib.request

from astropy.io.votable import parse


def utc():
    return datetime.now(timezone.utc).isoformat()


def retrieve(url, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    sidecar = path.with_suffix(path.suffix + ".provenance.json")
    if path.exists() and sidecar.exists():
        body = path.read_bytes()
        saved = json.loads(sidecar.read_text())
        if saved["url"] != url or saved["sha256"] != hashlib.sha256(body).hexdigest():
            raise ValueError(f"Cache provenance mismatch: {path}")
        return body
    with urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "ad-astra-research/0.2 metadata-audit"}), timeout=90) as response:
        body = response.read()
        status = response.status
        content_type = response.headers.get("Content-Type")
    path.write_bytes(body)
    sidecar.write_text(json.dumps({"retrieved_utc": utc(), "url": url, "status": status,
                                  "content_type": content_type, "sha256": hashlib.sha256(body).hexdigest(),
                                  "bytes": len(body)}, indent=2) + "\n")
    return body


def primitive(value):
    if getattr(value, "mask", False):
        return None
    return value.item() if hasattr(value, "item") else value


def audit_variables(target, output):
    tic, metadata = target["tic"], target["metadata"]
    requested_catalogs = {"II/366/catv2021", "B/vsx/vsx", "B/gcvs/gcvs_cat"}
    query = {"-source": ",".join(sorted(requested_catalogs)),
             "-c": f"{metadata['ra']} {metadata['dec']}", "-c.rs": "10", "-out.max": "100", "-out.all": "1"}
    url = "https://vizier.cds.unistra.fr/viz-bin/votable?" + urllib.parse.urlencode(query)
    path = output / "variable-cones" / f"{tic}.xml"
    result = {"tic": tic, "ra": metadata["ra"], "dec": metadata["dec"],
              "query_radius_arcsec": 10, "target_match_radius_arcsec": 5, "raw_path": str(path)}
    try:
        body = retrieve(url, path)
        document = parse(io.BytesIO(body), verify="exception")
        infos = list(document.infos) + [info for resource in document.resources for info in resource.infos]
        errors = [str(info.value) for info in infos
                  if (info.name or "").lower() == "error" or (info.value or "").upper() in {"ERROR", "OVERFLOW"}]
        if errors:
            raise ValueError(f"VizieR errors/overflow: {errors}")
        tables = []
        for table in document.iter_tables():
            rows = [{name: primitive(row[name]) for name in table.array.dtype.names} for row in table.array]
            if len(rows) >= 100:
                raise ValueError("VizieR table reached row cap")
            tables.append({"catalog": table.name, "description": table.description, "rows": rows})
        returned = {table["catalog"] for table in tables}
        missing = sorted(requested_catalogs - returned)
        result.update(status="complete" if not missing else "partial_catalog_coverage", tables=tables,
                      missing_catalog_tables=missing,
                      coverage_note="An absent requested table is unchecked, not a zero-match result")
    except Exception as error:
        result.update(status="failed", error=f"{type(error).__name__}: {error}")
    print(f"Variable catalogs TIC={tic} {result['status']} rows={sum(len(t['rows']) for t in result.get('tables', []))}", flush=True)
    return result


def audit_gaia(targets, output):
    clauses = [f"1=CONTAINS(POINT('ICRS',g.ra,g.dec),CIRCLE('ICRS',{row['metadata']['ra']},{row['metadata']['dec']},0.001388888888889))" for row in targets]
    query = ("SELECT g.source_id,g.ra,g.dec,g.phot_variable_flag,c.best_class_name,c.best_class_score,"
             "r.best_rotation_period,r.best_rotation_period_error FROM gaiadr3.gaia_source AS g "
             "LEFT JOIN gaiadr3.vari_classifier_result AS c ON g.source_id=c.source_id "
             "LEFT JOIN gaiadr3.vari_rotation_modulation AS r ON g.source_id=r.source_id WHERE " + " OR ".join(clauses))
    url = "https://gea.esac.esa.int/tap-server/tap/sync?" + urllib.parse.urlencode({"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "json", "QUERY": query})
    key = "-".join(str(row["tic"]) for row in targets)
    path = output / "gaia-variable-cones" / f"{key}.json"
    result = {"requested_tics": [row["tic"] for row in targets], "raw_path": str(path), "query_radius_arcsec": 5}
    try:
        payload = json.loads(retrieve(url, path))
        names = [column["name"] for column in payload["metadata"]]
        result.update(status="complete", rows=[dict(zip(names, row)) for row in payload["data"]])
    except Exception as error:
        result.update(status="failed", error=f"{type(error).__name__}: {error}")
    print(f"Gaia catalog batch {result['requested_tics']} {result['status']}", flush=True)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=Path("data/qlp-prior-review"))
    parser.add_argument("--summary", default="catalog-review.json")
    args = parser.parse_args()
    targets = json.loads(args.manifest.read_text())["targets"]
    declaration = {"created_utc": utc(), "manifest": str(args.manifest),
                   "manifest_sha256": hashlib.sha256(args.manifest.read_bytes()).hexdigest(),
                   "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                   "target_tics": [row["tic"] for row in targets],
                   "scope": "ASAS-SN VizieR catv2021 snapshot; VSX; GCVS, 10 arcsec; Gaia DR3 variability class and rotation catalogs, 5 arcsec. Coordinate matches are candidates for identity review, not automatic identities. No epoch propagation."}
    args.output.mkdir(parents=True, exist_ok=True)
    with ThreadPoolExecutor(max_workers=4) as pool:
        variables = list(pool.map(lambda target: audit_variables(target, args.output), targets))
        batches = [targets[index:index + 6] for index in range(0, len(targets), 6)]
        gaia = list(pool.map(lambda batch: audit_gaia(batch, args.output), batches))
    result = {"declaration": declaration, "completed_utc": utc(), "variable_catalogs": variables, "gaia_catalogs": gaia,
              "limitations": "Catalog non-matches do not establish novelty. ASAS-SN snapshot is dated in catalog metadata; SIMBAD indexing and Gaia classifications are incomplete, and blending/astrometry can misidentify a cone match. No ADS. No scientific light curves or images downloaded."}
    path = args.output / args.summary
    if path.exists():
        raise ValueError(f"Refusing to overwrite prior summary {path}")
    path.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
