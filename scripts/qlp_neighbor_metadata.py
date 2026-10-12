#!/usr/bin/env python3
"""Audit published periods of cataloged neighbors and target spectroscopy; no flux access."""
from __future__ import annotations

import io
import json
from pathlib import Path
import urllib.parse
import warnings

from astropy.io.votable import parse
from astropy.io.votable.exceptions import W03, W06, W50

from adastra.archive import sha256_file, utc_now, write_json
import qlp_stellar_metadata as base
import review_qlp_catalogs as prior

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data/qlp-stellar-metadata/completion"


def parse_spectroscopy(body):
    """Record these reviewed catalog-metadata defects without changing any values."""
    with warnings.catch_warnings(record=True) as recorded:
        warnings.simplefilter("always")
        document = parse(io.BytesIO(body), verify="warn")
    accepted = []
    for warning in recorded:
        message = str(warning.message)
        allowed = (
            warning.category is W03 and "Implicitly generating an ID from a name" in message
            or warning.category is W06 and message.endswith(
                "W06: Invalid UCD 'obs.field': Secondary word 'obs.field' is not valid as a primary word")
            or warning.category is W50 and any(message.endswith(f"W50: Invalid unit string '{unit}'")
                                              for unit in ("%", "log(cm.s**-2)"))
        )
        if not allowed:
            raise ValueError(f"Unreviewed spectroscopy metadata warning: {message}")
        accepted.append({"code": warning.category.__name__, "message": message})
    infos = list(document.infos) + [i for r in document.resources for i in r.infos]
    if any((i.name or "").lower() == "error" or (i.value or "").upper() in {"ERROR", "OVERFLOW"}
           for i in infos):
        raise ValueError("Spectroscopy catalog error or overflow")
    tables = []
    for table in document.iter_tables():
        if len(table.array) >= 100:
            raise ValueError("Spectroscopy query reached row cap")
        tables.append({"catalog": table.name, "description": table.description,
                       "fields": [{"name": f.name, "unit": str(f.unit), "description": f.description}
                                  for f in table.fields],
                       "rows": [{name: prior.primitive(row[name]) for name in table.array.dtype.names}
                                for row in table.array]})
    return {"recorded_metadata_warnings": accepted, "tables": tables}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    metadata_path = ROOT / "data/qlp-stellar-metadata/summary.json"
    metadata = json.loads(metadata_path.read_text())
    write_json(OUT / "declaration.json", {
        "created_utc": utc_now(), "source_sha256": sha256_file(Path(__file__)),
        "metadata_sha256": sha256_file(metadata_path),
        "purpose": __doc__, "chronology": "Exploratory follow-up after the QLP pilot and wider Gaia cone were inspected.",
        "helper_sources_sha256": {str(Path(module.__file__).relative_to(ROOT)): sha256_file(Path(module.__file__))
                                  for module in (base, prior)},
    })
    rr_ids = sorted(obj["source_id"] for report in metadata["reports"] for obj in report["sources"]
                    if obj["phot_variable_flag"] == "VARIABLE"
                    and any(row.get("best_class_name") == "RR"
                            for row in obj["gaia_variability_tables"]["vari_classifier_result"]))
    query = "SELECT * FROM gaiadr3.vari_rrlyrae WHERE source_id IN (" + ",".join(map(str, rr_ids)) + ")"
    rr = base.tap(base.GAIA, query, OUT / "gaia-neighbor-rrlyrae.json")
    report_rows = []
    for report in metadata["reports"]:
        tic, target = report["tic"], report["target"]
        query = {"-source": "J/MNRAS/506/150/stars,J/other/PASA/42.51/allstar",
                 "-c": f"{target['ra']} {target['dec']}", "-c.rs": 5,
                 "-out.all": "1", "-out.max": 100}
        url = "https://vizier.cds.unistra.fr/viz-bin/votable?" + urllib.parse.urlencode(query)
        path = OUT / f"galah-{tic}.xml"
        spectra = parse_spectroscopy(prior.retrieve(url, path))
        for table in spectra["tables"]:
            for row in table["rows"]:
                gaia_id = row.get("GaiaDR3", row.get("GaiaEDR3"))
                row["exact_target_gaia_identity"] = gaia_id == report["exact_gaia_dr3_source_id"]
        requested = {"J/MNRAS/506/150/stars", "J/other/PASA/42.51/allstar"}
        spectra["missing_requested_tables"] = sorted(requested - {t["catalog"] for t in spectra["tables"]})
        write_json(OUT / f"galah-{tic}-parsed.json", spectra)
        neighbors = []
        for source in report["sources"]:
            published = next((r for r in rr if r["source_id"] == source["source_id"]), None)
            if published is None:
                continue
            selected = {key: published[key] for key in ("source_id", "pf", "pf_error", "p1_o",
                        "best_classification", "peak_to_peak_g", "peak_to_peak_rp", "epoch_g")}
            selected.update(separation_arcsec=source["separation_arcsec"],
                            gaia_g_mean_mag=source["phot_g_mean_mag"],
                            conditional_required_intrinsic_semiamplitude_fraction_by_sector=
                            source["conditional_required_intrinsic_semiamplitude_fraction_by_sector"])
            if published["pf"] and 1 / published["pf"] > 1:
                selected["one_cycle_per_day_difference_alias_days"] = 1 / (1 / published["pf"] - 1)
            neighbors.append(selected)
        report_rows.append({"tic": tic, "target_gaia_dr3": report["exact_gaia_dr3_source_id"],
                            "published_rrlyrae_neighbors": neighbors,
                            "conditional_equatorial_rotation_speed_kms": report["conditional_rotation_speed_kms"],
                            "spectroscopy_path": str((OUT / f"galah-{tic}-parsed.json").relative_to(ROOT))})
    write_json(OUT / "summary.json", {
        "completed_utc": utc_now(), "reports": report_rows,
        "limitations": [
            "This is post-selection exploratory metadata follow-up; it retrieves no light curves or images.",
            "Known neighbor periods identify contamination hypotheses; they do not independently localize the TESS source.",
            "One-cycle-per-day frequency arithmetic does not prove a sampling alias, particularly for nearly continuous TESS observations.",
            "Gaia epoch_g time origin is not inferred here; no cross-mission phase comparison uses it.",
            "GALAH broadening includes rotational and other line broadening; conditional equatorial speeds use TIC radii and assume the photometric period is rotation.",
            "Spectroscopy flags and per-row exact Gaia source identities must be considered with all numeric values.",
            "Catalog non-matches or missing classifications never establish constancy or novelty.",
        ],
    })
    print(json.dumps(report_rows, indent=2))


if __name__ == "__main__":
    main()
