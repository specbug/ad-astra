#!/usr/bin/env python3
"""Public metadata audit of two exploratory QLP leads; never retrieves flux."""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import json
import math
from pathlib import Path
import re

import requests

from adastra.archive import sha256_file, utc_now, write_json

ROOT = Path(__file__).resolve().parents[1]
GAIA = "https://gea.esac.esa.int/tap-server/tap/sync"
SIMBAD = "https://simbad.cds.unistra.fr/simbad/sim-tap/sync"
LEADS = {149392671: {"period_days": 0.477, "semiamplitudes_ppm": [283, 336]},
         150104580: {"period_days": 1.20, "semiamplitudes_ppm": [107, 101]}}


def tap(endpoint, query, path):
    provenance_path = path.with_suffix(".provenance.json")
    if path.exists():
        provenance = json.loads(provenance_path.read_text())
        if (provenance["query"] != query or provenance["endpoint"] != endpoint
                or provenance["sha256"] != sha256_file(path)):
            raise ValueError(f"Cache provenance mismatch: {path}")
        payload = json.loads(path.read_text())
    else:
        response = requests.get(endpoint, params={"REQUEST": "doQuery", "LANG": "ADQL",
                               "FORMAT": "json", "MAXREC": 100000, "QUERY": query}, timeout=120)
        response.raise_for_status()
        payload = response.json()
        if "metadata" not in payload or "data" not in payload:
            raise ValueError("TAP JSON lacks metadata/data")
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(response.content)
        write_json(provenance_path, {"retrieved_utc": utc_now(), "endpoint": endpoint,
                   "resolved_url": response.url, "query": query,
                   "sha256": sha256_file(path), "bytes": path.stat().st_size,
                   "tls_verification": "requests defaults; verified", "rows": len(payload["data"])})
    names = [column["name"] for column in payload["metadata"]]
    return [dict(zip(names, row)) for row in payload["data"]]


def separation(ra1, dec1, ra2, dec2):
    a, b = math.radians(dec1), math.radians(dec2)
    value = math.sin((b-a)/2)**2 + math.cos(a)*math.cos(b)*math.sin(math.radians(ra2-ra1)/2)**2
    return 3600 * math.degrees(2*math.asin(min(1, math.sqrt(value))))


def audit(row, output):
    tic, meta = row["tic"], row["metadata"]
    directory = output / str(tic)
    directory.mkdir(parents=True, exist_ok=True)
    prior_audit_path = ROOT / row["simbad_audit"]
    prior_audit = json.loads(prior_audit_path.read_text())
    exact = [obj for obj in prior_audit["objects"] if f"TIC {tic}" in obj["ids"].split("|")]
    if len(exact) != 1:
        raise ValueError("A unique exact TIC alias is required")
    aliases = exact[0]["ids"].split("|")
    source_ids = [int(alias.split()[-1]) for alias in aliases if re.fullmatch(r"Gaia DR3 \d+", alias)]
    if len(source_ids) != 1:
        raise ValueError("A unique exact Gaia DR3 alias is required")
    target_id = source_ids[0]
    source = tap(GAIA, f"SELECT * FROM gaiadr3.gaia_source WHERE source_id={target_id}", directory / "gaia-target.json")
    if len(source) != 1 or source[0]["source_id"] != target_id:
        raise ValueError("Exact Gaia source missing or ambiguous")
    target = source[0]
    condition = ("1=CONTAINS(POINT('ICRS',ra,dec),"
                 f"CIRCLE('ICRS',{target['ra']},{target['dec']},0.025))")
    columns = ("source_id,ra,dec,ref_epoch,parallax,parallax_error,pmra,pmdec,phot_g_mean_mag,"
               "phot_bp_mean_mag,phot_rp_mean_mag,bp_rp,ruwe,astrometric_excess_noise,"
               "astrometric_excess_noise_sig,duplicated_source,phot_variable_flag,non_single_star,"
               "radial_velocity,radial_velocity_error,rv_nb_transits,rv_chisq_pvalue,"
               "rv_renormalised_gof,vbroad,vbroad_error,phot_g_n_obs,phot_g_mean_flux_over_error")
    neighbors = tap(GAIA, f"SELECT {columns} FROM gaiadr3.gaia_source WHERE {condition}", directory / "gaia-cone.json")
    count = tap(GAIA, f"SELECT COUNT(*) AS n FROM gaiadr3.gaia_source WHERE {condition}", directory / "gaia-cone-count.json")[0]["n"]
    if count != len(neighbors):
        raise ValueError("Gaia cone row count differs from independent count")
    ids = ",".join(str(obj["source_id"]) for obj in neighbors)
    auxiliary = {}
    failures = {}
    for name in ("vari_summary", "vari_classifier_result", "vari_rotation_modulation", "vari_short_timescale", "vari_eclipsing_binary"):
        try:
            auxiliary[name] = tap(GAIA, f"SELECT * FROM gaiadr3.{name} WHERE source_id IN ({ids})", directory / f"gaia-{name}.json")
        except Exception as error:
            failures[name] = f"{type(error).__name__}: {error}"
    simbad = tap(SIMBAD, "SELECT b.main_id,b.otype,b.ra,b.dec,b.oid,ids.ids,alltypes.otypes "
                 "FROM basic AS b LEFT JOIN ids ON b.oid=ids.oidref "
                 "LEFT JOIN alltypes ON b.oid=alltypes.oidref WHERE "
                 f"1=CONTAINS(POINT('ICRS',b.ra,b.dec),CIRCLE('ICRS',{target['ra']},{target['dec']},0.025))",
                 directory / "simbad-cone.json")
    for obj in simbad:
        obj["separation_arcsec_from_gaia_target_2016"] = separation(target["ra"], target["dec"], obj["ra"], obj["dec"])
    for obj in neighbors:
        obj["is_exact_target"] = obj["source_id"] == target_id
        obj["separation_arcsec"] = separation(target["ra"], target["dec"], obj["ra"], obj["dec"])
        obj["separation_tess_pixels_21arcsec"] = obj["separation_arcsec"] / 21
        if obj["phot_g_mean_mag"] is not None:
            ratio = 10 ** (-0.4 * (obj["phot_g_mean_mag"] - target["phot_g_mean_mag"]))
            obj["gaia_g_flux_ratio_to_target"] = ratio
            if not obj["is_exact_target"]:
                obj["conditional_required_intrinsic_semiamplitude_fraction_by_sector"] = [
                    amp * 1e-6 * (1+ratio)/ratio for amp in LEADS[tic]["semiamplitudes_ppm"]]
        obj["gaia_variability_tables"] = {name: [record for record in records if record["source_id"] == obj["source_id"]]
                                          for name, records in auxiliary.items()}
    neighbors.sort(key=lambda obj: obj["separation_arcsec"])
    report = {"created_utc": utc_now(), "tic": tic, "exact_gaia_dr3_source_id": target_id,
              "identity_evidence": {"simbad_main_id": exact[0]["main_id"], "aliases": aliases,
                                    "audit_path": str(prior_audit_path.relative_to(ROOT)),
                                    "audit_sha256": sha256_file(prior_audit_path)},
              "discovery_signal_rounded_for_context": LEADS[tic], "target": target,
              "cone_radius_arcsec": 90, "cone_reference_epoch": target["ref_epoch"],
              "cone_count_verified": count, "sources": neighbors, "simbad_sources": simbad,
              "auxiliary_query_failures": failures,
              "conditional_rotation_speed_kms": 2*math.pi*meta["rad"]*695700/(LEADS[tic]["period_days"]*86400),
              "limitations": [
                  "No flux products are retrieved or read; no protected season is accessed.",
                  "Gaia G and TESS have different passbands, and aperture transmission is unknown.",
                  "Conditional blend semiamplitudes assume only target+one neighbor, equal throughput and Gaia-G flux ratio in the TESS band; other flux raises the required amplitude under those assumptions.",
                  "Throughput or color differences can change the bound in either direction; these values do not localize or exclude a contaminant.",
                  "A missing Gaia variability row or NOT_AVAILABLE flag is not evidence of constancy.",
                  "Catalog RV constancy, RUWE, and vbroad are diagnostics, not a definitive variability class or rotation measurement.",
                  "The rotation speed assumes TIC radius and that the exploratory photometric period is the equatorial rotation period.",
                  "SIMBAD cone coordinates can have different reference epochs; exact aliases identify target, while separation alone does not establish crossmatches.",
              ]}
    write_json(directory / "report.json", report)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=ROOT / "data/qlp-frozen/manifest.json")
    parser.add_argument("--output", type=Path, default=ROOT / "data/qlp-stellar-metadata")
    args = parser.parse_args()
    manifest = json.loads(args.manifest.read_text())
    rows = [row for row in manifest["targets"] if row["tic"] in LEADS]
    if {row["tic"] for row in rows} != set(LEADS):
        raise ValueError("Manifest lacks either lead")
    args.output.mkdir(parents=True, exist_ok=True)
    declaration = {"created_utc": utc_now(), "source_sha256": sha256_file(Path(__file__)),
                   "manifest_sha256": sha256_file(args.manifest), "purpose": __doc__, "leads": LEADS}
    write_json(args.output / "declaration.json", declaration)
    with ThreadPoolExecutor(max_workers=2) as executor:
        reports = list(executor.map(lambda row: audit(row, args.output), rows))
    write_json(args.output / "summary.json", {"completed_utc": utc_now(), "reports": reports})
    for report in reports:
        target = report["target"]
        print(json.dumps({"tic": report["tic"], "sources": report["cone_count_verified"],
                          "ruwe": target["ruwe"], "rv_chisq_pvalue": target["rv_chisq_pvalue"],
                          "vbroad": target["vbroad"], "failures": report["auxiliary_query_failures"]}), flush=True)


if __name__ == "__main__":
    main()
