#!/usr/bin/env python3
"""Describe an existing S1/S2 pilot without downloading or opening holdouts.

All rankings are exploratory.  Fixed S1 ephemerides, the strongest peaks, and
every period-matched peak pair are retained even when their diagnostics fail.
No period/epoch optimization, clipping, or additional signal search occurs.
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import tempfile

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "ad-astra-matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(Path(tempfile.gettempdir()) / "ad-astra-cache"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

from adastra.archive import sha256_file, utc_now, write_json
from adastra.qlp import read_qlp
from adastra.screen import ScreenConfig, clean_curve, compare_sectors, event_windows, flatten_curve


ROOT = Path(__file__).resolve().parents[1]
LIMITATIONS = [
    "Exploratory inspection priorities, not detections, classifications, or novelty claims.",
    "Both discovery sectors already participated in screening; these diagnostics are not untouched validation.",
    "Positive event depths and odd/even differences have no calibrated statistical significance.",
    "An S1 ephemeris is held fixed here, but no timing uncertainty is propagated into S2.",
    "SAP has no assigned KSPSAP uncertainty; SAP depths are unweighted local median differences.",
    "QLP upstream processing and this screen's running median can change astrophysical shapes.",
    "A period match can arise from blends, spacecraft systematics, harmonics, or chance.",
    "Sign inversion is an asymmetric-noise check, not a calibrated false-positive rate.",
    "Historical cohort overlap is unknown. Protected holdout flux is never opened.",
]
RANKING = [
    "First: direct period-matched, event-triaged pairs with median S2 center separation <= half the sum of durations and >= 3 positive fixed-ephemeris S2 event depths.",
    "Second: other direct period matches; third: LS period agreement; fourth: harmonic-only matches; last: remaining targets.",
    "Within a tier, sort by the number of positive fixed-ephemeris S2 depths for the displayed S1 peak, then TIC ID.",
    "These thresholds are descriptive inspection rules; they do not imply validation or calibrated likelihood.",
]


def phase_days(time, period, epoch):
    return (np.asarray(time) - epoch + period / 2) % period - period / 2


def depth_robustness(events):
    """Summarize local event medians without pretending they are independent."""
    valid = [event for event in events["events"] if event["local_median_depth_ppm"] is not None]
    depths = np.asarray([event["local_median_depth_ppm"] for event in valid], dtype=float)
    parity = {}
    for remainder, name in ((0, "even"), (1, "odd")):
        group = [event["local_median_depth_ppm"] for event in valid if event["index"] % 2 == remainder]
        parity[name] = {"measurable_events": len(group), "median_event_depth_ppm": float(np.median(group)) if group else None}
    positive = depths[depths > 0]
    loo = [float(np.mean(np.delete(depths, index))) for index in range(len(depths))] if len(depths) > 1 else []
    return {
        "sampled_windows": events["sampled_windows"], "measurable_windows": len(valid),
        "positive_windows": int(np.count_nonzero(depths > 0)),
        "positive_fraction": float(np.mean(depths > 0)) if len(depths) else None,
        "median_event_depth_ppm": float(np.median(depths)) if len(depths) else None,
        "minimum_event_depth_ppm": float(np.min(depths)) if len(depths) else None,
        "maximum_event_depth_ppm": float(np.max(depths)) if len(depths) else None,
        "odd_even": parity,
        "odd_minus_even_median_depth_ppm": (parity["odd"]["median_event_depth_ppm"] - parity["even"]["median_event_depth_ppm"]
                                                if all(parity[name]["median_event_depth_ppm"] is not None for name in ("odd", "even")) else None),
        "largest_positive_depth_fraction_of_positive_depth_sum": float(np.max(positive) / np.sum(positive)) if len(positive) else None,
        "leave_one_event_out_mean_depth_range_ppm": [min(loo), max(loo)] if loo else None,
        "interpretation": "Unweighted summaries of local event sample medians; depth-sum fraction is not a share of integrated signal or significance.",
    }


def fixed_ephemeris(curve, peak):
    result = {"period_days": peak["period_days"], "epoch_btjd": peak["epoch_btjd"], "duration_days": peak["duration_days"]}
    for name, time, flux in (("kspsap_flattened", curve["t"], curve["flat"]),
                             ("kspsap_no_additional_flattening", curve["t"], curve["f"]),
                             ("sap", curve["sap_t"], curve["sap_f"])):
        if len(time):
            events = event_windows(time, flux, peak["period_days"], peak["epoch_btjd"], peak["duration_days"])
            result[name] = {"events": events, "depth_robustness": depth_robustness(events)}
        else:
            result[name] = {"status": "no_usable_samples"}
    return result


def pair_phase(first, second, time):
    """Compare independently selected centers across S2; no phase fitting."""
    p1, p2 = first["period_days"], second["period_days"]
    start = int(np.ceil((time[0] - second["epoch_btjd"]) / p2))
    end = int(np.floor((time[-1] - second["epoch_btjd"]) / p2))
    centers = second["epoch_btjd"] + np.arange(start, end + 1) * p2
    offsets = phase_days(centers, p1, first["epoch_btjd"])
    threshold = (first["duration_days"] + second["duration_days"]) / 2
    at_midpoint = second["epoch_btjd"] + np.round((np.median(time) - second["epoch_btjd"]) / p2) * p2
    return {
        "comparison_centers_btjd": centers.tolist(), "signed_offset_from_nearest_s1_prediction_days": offsets.tolist(),
        "median_absolute_center_offset_days": float(np.median(np.abs(offsets))) if len(offsets) else None,
        "maximum_absolute_center_offset_days": float(np.max(np.abs(offsets))) if len(offsets) else None,
        "offset_at_s2_center_days": float(phase_days(at_midpoint, p1, first["epoch_btjd"])),
        "window_overlap_tolerance_days": threshold,
        "fraction_centers_with_window_overlap": float(np.mean(np.abs(offsets) <= threshold)) if len(offsets) else None,
        "median_centers_overlap": bool(len(offsets) and np.median(np.abs(offsets)) <= threshold),
        "interpretation": "Geometric center differences of independent box fits, not a statistical coherence test; no ephemeris uncertainty is propagated.",
    }


def load_cached_sector(row, tic, repo_root):
    if type(row.get("sector")) is not int or row["sector"] not in (1, 2):
        raise ValueError("Diagnostics may open literal sectors 1 and 2 only")
    if row.get("status") != "complete":
        raise ValueError("Cannot open an incomplete sector result")
    uri = row["uri"]
    match = re.fullmatch(r"mast:HLSP/qlp/s(\d{4})/(\d{4})/(\d{4})/(\d{4})/(\d{4})/hlsp_qlp_tess_ffi_s(\d{4})-(\d{16})_tess_v\d+_llc\.fits", uri)
    if (match is None or int(match[1]) != row["sector"] or int(match[6]) != row["sector"]
            or int("".join(match.group(index) for index in (2, 3, 4, 5))) != tic or int(match[7]) != tic):
        raise ValueError("Recorded URI does not identify this exact S1/S2 product")
    raw_path = Path(row["input_path"])
    path = raw_path if raw_path.is_absolute() else repo_root / raw_path
    if path.name != hashlib.sha256(uri.encode()).hexdigest() + ".fits":
        raise ValueError("Cached filename does not bind this exact archive URI")
    expected = row["input_sha256"]
    if sha256_file(path) != expected:
        raise ValueError("Cached input differs from the pilot input hash")
    provenance = json.loads(path.with_suffix(".json").read_text())
    if provenance.get("sha256") != expected or provenance.get("uri") != row["uri"]:
        raise ValueError("Cached provenance differs from the pilot input")
    lc, sap = read_qlp(path), read_qlp(path, "SAP_FLUX")
    for product in (lc, sap):
        if product.metadata["tic_id"] != tic or product.metadata["sector"] != row["sector"]:
            raise ValueError("Cached FITS identity differs from the pilot target/sector")
    config = ScreenConfig(**row["analysis"]["config"])
    t, f, error, _ = clean_curve(lc.time_btjd, lc.flux, lc.flux_error, lc.quality)
    flat, _, trend, _ = flatten_curve(t, f, error, config)
    good = sap.good_mask()
    sap_t, sap_f = sap.time_btjd[good], sap.flux[good]
    order = np.argsort(sap_t, kind="stable")
    sap_t, sap_f = sap_t[order], sap_f[order]
    keep = np.r_[True, np.diff(sap_t) > 0] if len(sap_t) else np.empty(0, dtype=bool)
    sap_t, sap_f = sap_t[keep], sap_f[keep]
    if len(sap_f):
        sap_f = sap_f / np.median(sap_f)
    if sha256_file(path) != expected:
        raise ValueError("Cached input changed while reading")
    return {"t": t, "f": f, "error": error, "flat": flat, "trend": trend,
            "sap_t": sap_t, "sap_f": sap_f, "path": str(path), "sha256": expected, "sector": row["sector"]}


def summarize_target(row, curves):
    result = {"tic": row["tic"], "status": row["status"], "tags": [], "limitations": LIMITATIONS,
              "sector_summaries": [], "fixed_s1_ephemerides": [], "matched_pair_phases": []}
    if row["status"] != "complete" or set(curves) != {1, 2}:
        result.update(status="incomplete", errors=[sector.get("error") for sector in row.get("sectors", []) if sector.get("status") != "complete"], inspection_tier=9)
        return result
    sectors = {sector["sector"]: sector for sector in row["sectors"]}
    comparison = row["comparison"]
    for number, sector in sorted(sectors.items()):
        result["sector_summaries"].append({"sector": number, "input_path": curves[number]["path"],
            "input_sha256": curves[number]["sha256"], "bls_peaks": sector["analysis"]["bls_peaks"],
            "lomb_scargle": sector["analysis"]["lomb_scargle"], "inverted_analysis": sector.get("inverted_analysis")})
    first_peaks = {peak["rank"]: peak for peak in sectors[1]["analysis"]["bls_peaks"]}
    second_peaks = {peak["rank"]: peak for peak in sectors[2]["analysis"]["bls_peaks"]}
    for rank, peak in first_peaks.items():
        result["fixed_s1_ephemerides"].append({"s1_rank": rank, "sector1": fixed_ephemeris(curves[1], peak),
                                              "sector2": fixed_ephemeris(curves[2], peak)})
    fixed = {item["s1_rank"]: item for item in result["fixed_s1_ephemerides"]}
    coherent = []
    tolerance = sectors[1]["analysis"]["config"]["period_agreement_fraction"]
    for pair in comparison["pairs"]:
        period_matched = pair["fractional_direct_difference"] <= tolerance
        harmonic_matched = pair["fractional_harmonic_difference"] <= tolerance
        if period_matched or harmonic_matched:
            copy = dict(pair)
            copy["direct_period_match_before_event_triage"] = period_matched
            copy["harmonic_period_match_before_event_triage"] = harmonic_matched
            if period_matched:
                copy["phase"] = pair_phase(first_peaks[pair["first_rank"]], second_peaks[pair["second_rank"]], curves[2]["t"])
                positives = fixed[pair["first_rank"]]["sector2"]["kspsap_flattened"]["depth_robustness"]["positive_windows"]
                copy["exploratory_phase_and_positive_window_priority"] = bool(pair["direct_inspection_priority"] and copy["phase"]["median_centers_overlap"] and positives >= 3)
                if copy["exploratory_phase_and_positive_window_priority"]:
                    coherent.append(copy)
            result["matched_pair_phases"].append(copy)
    direct = [pair for pair in comparison["pairs"] if pair["direct_inspection_priority"]]
    harmonic = [pair for pair in comparison["pairs"] if pair["harmonic_inspection_priority"]]
    candidates = coherent or direct or harmonic
    if candidates:
        representative = min(candidates, key=lambda pair: (pair["fractional_direct_difference"] if pair["direct_inspection_priority"] else pair["fractional_harmonic_difference"], pair["first_rank"], pair["second_rank"]))
        displayed_rank = representative["first_rank"]
        result["displayed_pair"] = {key: representative[key] for key in ("first_rank", "second_rank", "first_period_days", "second_period_days")}
    else:
        displayed_rank = min(first_peaks) if first_peaks else None
        result["displayed_pair"] = None
    result["displayed_s1_rank"] = displayed_rank
    result["displayed_s1_peak"] = first_peaks.get(displayed_rank)
    result["comparison"] = comparison
    result["coherent_direct_pairs"] = len(coherent)
    if coherent:
        result["tags"].append("direct_period_match_with_geometric_phase_and_positive_window_support")
    elif direct:
        result["tags"].append("direct_period_match_fails_phase_or_fixed_positive_window_rule")
    if harmonic:
        result["tags"].append("harmonic_period_match")
    if comparison["ls_period_agreement"]:
        result["tags"].append("ls_period_match")
    if not result["tags"]:
        result["tags"].append("no_cross_sector_period_priority")
    result["inspection_tier"] = 0 if coherent else 1 if direct else 2 if comparison["ls_period_agreement"] else 3 if harmonic else 4
    result["displayed_fixed_s2_positive_windows"] = fixed[displayed_rank]["sector2"]["kspsap_flattened"]["depth_robustness"]["positive_windows"] if displayed_rank else 0
    inverted = [sectors[number].get("inverted_analysis") for number in (1, 2)]
    if all(item and item.get("status") == "complete" for item in inverted):
        result["inverted_comparison"] = compare_sectors(*inverted, config=sectors[1]["analysis"]["config"])
        if result["inverted_comparison"]["any_direct_inspection_priority"]:
            result["tags"].append("inverted_flux_also_has_direct_period_priority")
    else:
        result["inverted_comparison"] = None
    return result


def draw_series(axis, time, flux, title):
    axis.scatter(time, (flux - 1) * 1e6, s=1.5, alpha=0.6, rasterized=True)
    axis.set_title(title, fontsize=8)
    axis.set_xlabel("BTJD (TDB)", fontsize=7)
    axis.set_ylabel("Median-normalized ppm", fontsize=7)
    axis.tick_params(labelsize=7)


def draw_fold(axis, time, flux, period, epoch, title, duration=None):
    phase = phase_days(time, period, epoch) / period
    axis.scatter(phase, (flux - 1) * 1e6, s=2, alpha=0.25, color="#4c6a87", rasterized=True)
    edges = np.linspace(-0.5, 0.5, 61)
    assignments = np.searchsorted(edges, phase, side="right") - 1
    centers, medians = [], []
    for index in range(60):
        samples = flux[assignments == index]
        if len(samples) >= 3:
            centers.append((edges[index] + edges[index + 1]) / 2)
            medians.append((np.median(samples) - 1) * 1e6)
    axis.plot(centers, medians, ".", color="#af4e2f", markersize=3, label="60-bin medians, >=3 samples")
    if duration:
        axis.axvspan(-duration / period / 2, duration / period / 2, alpha=0.14, color="#e3b24a")
    axis.set_xlim(-0.5, 0.5)
    axis.set_title(title, fontsize=8)
    axis.set_xlabel("Phase from stated epoch (cycles)", fontsize=7)
    axis.set_ylabel("ppm", fontsize=7)
    axis.tick_params(labelsize=7)


def plot_target(summary, curves, output):
    if summary["status"] != "complete":
        return None
    figure, axes = plt.subplots(6, 2, figsize=(13, 17), constrained_layout=True)
    peak = summary["displayed_s1_peak"]
    fixed = next((item for item in summary["fixed_s1_ephemerides"] if item["s1_rank"] == summary["displayed_s1_rank"]), None)
    for column, number in enumerate((1, 2)):
        curve = curves[number]
        draw_series(axes[0, column], curve["sap_t"], curve["sap_f"], f"Sector {number}: SAP, no additional flattening")
        draw_series(axes[1, column], curve["t"], curve["f"], f"Sector {number}: KSPSAP, upstream corrections remain")
        draw_series(axes[2, column], curve["t"], curve["flat"], f"Sector {number}: KSPSAP / fixed-config running median")
        if peak:
            draw_fold(axes[3, column], curve["t"], curve["flat"], peak["period_days"], peak["epoch_btjd"],
                      f"S{number}, fixed S1 rank {peak['rank']} ephemeris: P={peak['period_days']:.7f} d", peak["duration_days"])
            event_list = fixed[f"sector{number}"]["kspsap_flattened"]["events"]["events"]
            measured = [item for item in event_list if item["local_median_depth_ppm"] is not None]
            axes[5, column].axhline(0, color="grey", linewidth=0.5)
            for parity, marker, label in ((0, "o", "even index"), (1, "x", "odd index")):
                subset = [item for item in measured if item["index"] % 2 == parity]
                axes[5, column].scatter([item["center_btjd"] for item in subset], [item["local_median_depth_ppm"] for item in subset], marker=marker, label=label, s=18)
            axes[5, column].legend(fontsize=7)
            axes[5, column].set_title(f"S{number}: local median depths at fixed S1 predictions", fontsize=8)
            axes[5, column].set_xlabel("Predicted center BTJD", fontsize=7)
            axes[5, column].set_ylabel("Depth ppm (no significance assigned)", fontsize=7)
            axes[5, column].tick_params(labelsize=7)
        ls = next(item["lomb_scargle"] for item in summary["sector_summaries"] if item["sector"] == number)
        if ls.get("period_days"):
            draw_fold(axes[4, column], curve["t"], curve["f"], ls["period_days"], curve["t"][0],
                      f"S{number}: independent strongest LS P={ls['period_days']:.7f} d (arbitrary phase zero)")
        else:
            axes[4, column].text(0.5, 0.5, "No LS peak", ha="center")
    figure.suptitle(f"TIC {summary['tic']} | inspection rank {summary.get('inspection_rank', '?')} | S1 rank {summary['displayed_s1_rank']} fixed in both sectors\n"
                     + "; ".join(summary["tags"]) + "\nExploratory. Full sample ranges shown; no clipping; no holdout data.", fontsize=10)
    path = output / "targets" / f"tic-{summary['tic']}.png"
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=145)
    plt.close(figure)
    return str(path)


def contact_sheets(summaries, all_curves, output):
    successful = [item for item in summaries if item["status"] == "complete"]
    paths = []
    with PdfPages(output / "all-targets-contact-sheets.pdf") as pdf:
        for offset in range(0, len(successful), 6):
            group = successful[offset:offset + 6]
            figure, axes = plt.subplots(len(group), 4, figsize=(18, 2.35 * len(group)), squeeze=False, constrained_layout=True)
            for row_number, summary in enumerate(group):
                for column, (sector, kind) in enumerate(((1, "SAP"), (1, "KSPSAP"), (2, "SAP"), (2, "KSPSAP"))):
                    curve = all_curves[summary["tic"]][sector]
                    time, flux = (curve["sap_t"], curve["sap_f"]) if kind == "SAP" else (curve["t"], curve["f"])
                    draw_series(axes[row_number, column], time, flux, f"#{summary['inspection_rank']} TIC {summary['tic']} S{sector} {kind}")
            figure.suptitle("Original archive SAP and KSPSAP, each median normalized; independent full-range y axes; no clipping", fontsize=11)
            path = output / f"contact-time-{offset // 6 + 1}.png"
            figure.savefig(path, dpi=145)
            pdf.savefig(figure)
            paths.append(str(path))
            plt.close(figure)
            figure, axes = plt.subplots(len(group), 2, figsize=(12, 2.35 * len(group)), squeeze=False, constrained_layout=True)
            for row_number, summary in enumerate(group):
                peak = summary["displayed_s1_peak"]
                for column, sector in enumerate((1, 2)):
                    curve = all_curves[summary["tic"]][sector]
                    if peak:
                        draw_fold(axes[row_number, column], curve["t"], curve["flat"], peak["period_days"], peak["epoch_btjd"],
                                  f"#{summary['inspection_rank']} TIC {summary['tic']} S{sector}: S1 rank {peak['rank']}, P={peak['period_days']:.6f}d", peak["duration_days"])
            figure.suptitle("Fixed S1 ephemeris applied to both sectors; full phase; 60-bin medians; exploratory only", fontsize=11)
            path = output / f"contact-fold-{offset // 6 + 1}.png"
            figure.savefig(path, dpi=145)
            pdf.savefig(figure)
            paths.append(str(path))
            plt.close(figure)
    return paths


def compact_row(item):
    sectors = {row["sector"]: row for row in item["sector_summaries"]}
    peak = item.get("displayed_s1_peak") or {}
    selected = next((row for row in item["fixed_s1_ephemerides"] if row["s1_rank"] == item.get("displayed_s1_rank")), None)
    fixed_depth = selected["sector2"]["kspsap_flattened"]["depth_robustness"] if selected else {}
    sap_depth = selected["sector2"]["sap"].get("depth_robustness", {}) if selected else {}
    return {"inspection_rank": item.get("inspection_rank"), "tic": item["tic"], "status": item["status"],
            "inspection_tier": item["inspection_tier"], "tags": "; ".join(item["tags"]),
            "displayed_s1_rank": item.get("displayed_s1_rank"), "s1_period_days": peak.get("period_days"),
            "s1_epoch_btjd": peak.get("epoch_btjd"), "s1_depth_ppm": peak.get("depth_ppm"),
            "fixed_s2_measurable_windows": fixed_depth.get("measurable_windows"),
            "fixed_s2_positive_windows": fixed_depth.get("positive_windows"),
            "fixed_s2_median_local_depth_ppm": fixed_depth.get("median_event_depth_ppm"),
            "fixed_s2_odd_minus_even_ppm": fixed_depth.get("odd_minus_even_median_depth_ppm"),
            "fixed_s2_largest_positive_depth_fraction": fixed_depth.get("largest_positive_depth_fraction_of_positive_depth_sum"),
            "fixed_s2_sap_median_local_depth_ppm": sap_depth.get("median_event_depth_ppm"),
            "coherent_direct_pairs": item.get("coherent_direct_pairs", 0),
            "s1_ls_period_days": sectors.get(1, {}).get("lomb_scargle", {}).get("period_days"),
            "s2_ls_period_days": sectors.get(2, {}).get("lomb_scargle", {}).get("period_days"),
            "s1_ls_amplitude_ppm": sectors.get(1, {}).get("lomb_scargle", {}).get("amplitude_ppm"),
            "s2_ls_amplitude_ppm": sectors.get(2, {}).get("lomb_scargle", {}).get("amplitude_ppm"),
            "plot_path": item.get("plot_path")}


def run(results_path, output, repo_root=ROOT):
    results_path, output, repo_root = map(Path, (results_path, output, repo_root))
    payload = json.loads(results_path.read_text())
    declaration = payload["declaration"]
    if payload.get("status") != "complete" or not payload.get("original_sources_unchanged") or not payload.get("all_requested_targets_accounted_for"):
        raise ValueError("Require a completed, unchanged, fully accounted-for pilot")
    if declaration.get("allowed_discovery_sectors") != [1, 2] or declaration.get("protected_holdouts_accessed") is not False:
        raise ValueError("Pilot declaration does not protect holdouts and literal sectors 1/2")
    for name in ("src/adastra/archive.py", "src/adastra/qlp.py", "src/adastra/screen.py"):
        if declaration.get("source_sha256", {}).get(name) != sha256_file(repo_root / name):
            raise ValueError("Current ingestion or flattening source differs from the completed pilot")
    targets = payload["results"]
    if len(targets) != payload["target_count"] or len({row["tic"] for row in targets}) != len(targets):
        raise ValueError("Duplicate or missing targets in completed pilot")
    # Validate all recorded sectors before opening any cached FITS.
    for row in targets:
        sector_numbers = [sector.get("sector") for sector in row.get("sectors", [])]
        if any(type(sector) is not int or sector not in (1, 2) for sector in sector_numbers) or len(set(sector_numbers)) != len(sector_numbers):
            raise ValueError("Recorded sector is outside literal S1/S2, or repeated")
    if output.exists() and any(output.iterdir()):
        raise ValueError("Diagnostics output directory is nonempty; refusing overwrite")
    declaration = {"created_utc": utc_now(), "stage": "exploratory diagnostics of completed screening, no holdout validation",
                   "input_results": str(results_path), "input_results_sha256": sha256_file(results_path),
                   "source_sha256": {"scripts/qlp_diagnostics.py": sha256_file(Path(__file__)),
                                     "src/adastra/archive.py": sha256_file(repo_root / "src/adastra/archive.py"),
                                     "src/adastra/qlp.py": sha256_file(repo_root / "src/adastra/qlp.py"),
                                     "src/adastra/screen.py": sha256_file(repo_root / "src/adastra/screen.py")},
                   "allowed_sectors": [1, 2], "protected_holdouts_accessed": False,
                   "ranking_rules": RANKING, "limitations": LIMITATIONS,
                   "phase_rule": "Direct-matched independent box centers overlap if separation <= half sum of durations; no retuning or uncertainty model.",
                   "plot_rules": "All quality-zero finite positive samples retained, no clipping; 60 equal phase bins with >=3 samples; axes retain full flux range."}
    write_json(output / "declaration.json", declaration)
    summaries, all_curves = [], {}
    for row in targets:
        curves = {}
        for sector in row.get("sectors", []):
            if sector.get("status") == "complete":
                curves[sector["sector"]] = load_cached_sector(sector, row["tic"], repo_root)
        all_curves[row["tic"]] = curves
        summaries.append(summarize_target(row, curves))
    summaries.sort(key=lambda item: (item["inspection_tier"], -item.get("displayed_fixed_s2_positive_windows", 0), item["tic"]))
    for rank, summary in enumerate(summaries, 1):
        summary["inspection_rank"] = rank
        summary["plot_path"] = plot_target(summary, all_curves[summary["tic"]], output)
        write_json(output / "targets" / f"tic-{summary['tic']}.json", summary)
        print(f"Diagnosed TIC {summary['tic']}: tier {summary['inspection_tier']} ({rank}/{len(summaries)})", flush=True)
    sheets = contact_sheets(summaries, all_curves, output)
    with (output / "summary.csv").open("w", newline="") as handle:
        rows = [compact_row(item) for item in summaries]
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    unchanged = sha256_file(results_path) == declaration["input_results_sha256"]
    for name, digest in declaration["source_sha256"].items():
        unchanged &= sha256_file(repo_root / name) == digest
    report = {"declaration": declaration, "completed_utc": utc_now(),
              "status": "complete" if unchanged else "invalid_changed_sources_or_results",
              "targets": summaries, "target_count": len(summaries),
              "tier_counts": {str(tier): sum(item["inspection_tier"] == tier for item in summaries) for tier in (0, 1, 2, 3, 4, 9)},
              "contact_sheets": sheets, "protected_holdouts_accessed": False}
    write_json(output / "summary.json", report)
    return report


def self_test():
    """Check phase geometry and local-event bookkeeping on synthetic arrays."""
    peak = {"period_days": 2.0, "epoch_btjd": 0.25, "duration_days": 0.2}
    time = np.arange(10.0, 20.0, 0.02)
    same = pair_phase(peak, {**peak, "epoch_btjd": 4.25}, time)
    shifted = pair_phase(peak, {**peak, "epoch_btjd": 4.75}, time)
    assert same["median_centers_overlap"] and same["maximum_absolute_center_offset_days"] < 1e-10
    assert not shifted["median_centers_overlap"] and abs(shifted["median_absolute_center_offset_days"] - 0.5) < 1e-10
    flux = np.ones(len(time))
    flux[np.abs(phase_days(time, 2.0, 0.25)) <= 0.1] -= 0.003
    events = event_windows(time, flux, 2.0, 0.25, 0.2)
    metrics = depth_robustness(events)
    assert metrics["positive_windows"] == 5 and abs(metrics["median_event_depth_ppm"] - 3000) < 1e-6
    assert abs(metrics["odd_minus_even_median_depth_ppm"]) < 1e-6
    assert depth_robustness({"sampled_windows": 0, "events": []})["median_event_depth_ppm"] is None
    try:
        load_cached_sector({"sector": 27, "status": "complete"}, 1, ROOT)
    except ValueError as error:
        assert "sectors 1 and 2" in str(error)
    else:
        raise AssertionError("Holdout sector must be rejected before any file access")
    peak = {**peak, "rank": 1, "depth_ppm": 3000.0, "passes_event_window_triage": True}
    config = vars(ScreenConfig())
    analyses = [{"status": "complete", "config": config, "bls_peaks": [dict(peak)],
                 "lomb_scargle": {"period_days": 2.0, "amplitude_ppm": 1000.0}} for _ in (1, 2)]
    row = {"tic": 999, "status": "complete", "sectors": [{"sector": sector, "analysis": analyses[sector - 1]} for sector in (1, 2)],
           "comparison": compare_sectors(*analyses)}
    curves = {sector: {"t": time - (10 if sector == 1 else 0), "f": flux, "flat": flux,
                       "sap_t": time - (10 if sector == 1 else 0), "sap_f": flux,
                       "path": "synthetic-only", "sha256": "synthetic-only"} for sector in (1, 2)}
    summary = summarize_target(row, curves)
    assert summary["inspection_tier"] == 0 and summary["displayed_fixed_s2_positive_windows"] == 5
    summary["inspection_rank"] = 1
    assert compact_row(summary)["fixed_s2_positive_windows"] == 5
    with tempfile.TemporaryDirectory(prefix="qlp-diagnostics-synthetic-") as temporary:
        output = Path(temporary)
        assert Path(plot_target(summary, curves, output)).exists()
        assert len(contact_sheets([summary], {999: curves}, output)) == 2
    analyses[1]["bls_peaks"][0]["epoch_btjd"] += 0.5
    mismatched = summarize_target(row, curves)
    assert mismatched["inspection_tier"] == 1 and mismatched["coherent_direct_pairs"] == 0
    print("Synthetic phase, event-depth, odd/even, empty-depth, holdout guards, ranking, and plot integration passed.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--repo-root", type=Path, default=ROOT)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
    else:
        if args.results is None or args.output is None:
            parser.error("--results and --output are required unless --self-test")
        run(args.results, args.output, args.repo_root)
