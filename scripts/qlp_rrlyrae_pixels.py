#!/usr/bin/env python3
"""Exploratory, fixed-period pixel check of a catalogued RR Lyrae neighbour.

Only sectors 1 and 2 are permitted. This is a post-selection contamination
diagnostic, not a new independent detection or a source-localization claim.
"""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import argparse
import json
import os
import zipfile
from urllib.parse import urlparse

import requests
import numpy as np
from astropy.time import Time
from astropy.wcs import WCS

from adastra.archive import sha256_file, utc_now, validated_fits, write_json

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data/qlp-rrlyrae-pixels"
OUT = ROOT / "runs/qlp-rrlyrae-pixels"
API = "https://mast.stsci.edu/tesscut/api/v0.1/astrocut"
PERIOD = 0.47711359612012993
RA, DEC = 84.265499658, -62.44919474
MAX_BYTES = 128 * 1024 * 1024


def fetch(sector):
    if type(sector) is not int or sector not in (1, 2):
        raise ValueError("Literal sector 1/2 guard protects the old holdouts")
    directory = DATA / f"s{sector:04d}"
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / "tesscut.zip"
    sidecar = directory / "provenance.json"
    parameters = dict(ra=RA, dec=DEC, x=7, y=7, units="px", sector=sector)
    if path.exists():
        saved = json.loads(sidecar.read_text())
        if saved["parameters"] != parameters or saved["zip_sha256"] != sha256_file(path):
            raise ValueError("Cached archive mismatch")
    else:
        part = path.with_suffix(".zip.part")
        if part.exists():
            raise ValueError("Preserve existing partial for explicit recovery")
        initial = requests.get(API, params=parameters, stream=True, timeout=(30, 150), allow_redirects=False)
        initial_url = initial.url
        delivery_url = initial.headers.get("Location")
        regional_delivery = False
        if initial.is_redirect:
            location = urlparse(delivery_url)
            if (location.scheme != "https" or location.hostname != "stpubdata.s3.amazonaws.com"
                    or location.query or not location.path.startswith("/tess/public/staged_cutouts/")):
                initial.close()
                raise ValueError("Unexpected or signed TESSCut delivery redirect")
            # The same public bucket's official regional endpoint is already in
            # this environment's NASA-data allowlist. Preserve both exact URLs.
            delivery_url = "https://stpubdata.s3.us-east-1.amazonaws.com" + location.path
            regional_delivery = True
            initial.close()
            response = requests.get(delivery_url, stream=True, timeout=(30, 120))
        else:
            response = initial
        with response:
            record = dict(requested_utc=utc_now(), parameters=parameters,
                          initial_url=initial_url, delivery_url=delivery_url,
                          regional_delivery=regional_delivery,
                          resolved_url=response.url, http_status=response.status_code)
            if not response.ok:
                record["error_body"] = response.text[:4000]
                write_json(directory / "error.json", record)
                response.raise_for_status()
            size = 0
            with part.open("xb") as stream:
                for block in response.iter_content(1024 * 1024):
                    size += len(block)
                    if size > MAX_BYTES:
                        raise ValueError("TESSCut archive exceeds explicit 128 MiB cap")
                    stream.write(block)
            if not zipfile.is_zipfile(part):
                write_json(directory / "error.json", {**record, "error": "Non-ZIP body",
                           "body": part.read_bytes()[:4000].decode(errors="replace")})
                raise ValueError("TESSCut did not return a ZIP archive")
            record.update(completed_utc=utc_now(), zip_bytes=size, zip_sha256=sha256_file(part))
            part.replace(path)
            write_json(sidecar, record)
    products = []
    with zipfile.ZipFile(path) as archive:
        members = [x for x in archive.infolist() if x.filename.lower().endswith(".fits")]
        if not members or sum(x.file_size for x in members) > MAX_BYTES:
            raise ValueError("Unexpected or oversized FITS archive")
        for member in members:
            # Flat basename only: never extract user/server supplied path components.
            target = directory / Path(member.filename).name
            payload = archive.read(member)
            if target.exists() and target.read_bytes() != payload:
                raise ValueError("Existing extracted artifact differs from preserved ZIP")
            if not target.exists():
                target.write_bytes(payload)
            with validated_fits(target) as hdus:
                sectors = {int(h.header["SECTOR"]) for h in hdus if "SECTOR" in h.header}
                if sectors != {sector}:
                    raise ValueError(f"Downloaded sector mismatch: {sectors}")
                integrity = [{"hdu": i, "checksum_supplied": "CHECKSUM" in h.header,
                              "datasum_supplied": "DATASUM" in h.header}
                             for i, h in enumerate(hdus)]
            products.append(dict(path=str(target.relative_to(ROOT)), sha256=sha256_file(target),
                                 original_zip_member=member.filename, integrity=integrity))
    write_json(directory / "products.json", products)
    return products


def fit_pixels(time, cube, error, background, period=PERIOD):
    """Independent pixel regressions, with common fixed harmonic basis.

    Formal standard errors are residual-scaled white-noise diagnostics only.
    No optimized frequency, pixel selection, clipping, or centroid is used.
    """
    labels = np.cumsum(np.r_[0, np.diff(time) > .5])
    baseline = []
    for label in np.unique(labels):
        member = labels == label
        centered_time = (time - np.mean(time[member])) / max(np.ptp(time[member]), 1.)
        baseline.extend([member.astype(float), member * centered_time])
    phase = 2 * np.pi * (time - time[0]) / period
    harmonics = [f(k * phase) for k in (1, 2, 3) for f in (np.cos, np.sin)]
    if background is not None:
        bg = background - np.median(background)
        if np.std(bg) > 0:
            baseline.append(bg / np.std(bg))
    design = np.array(baseline + harmonics).T
    output = []
    for y in range(cube.shape[1]):
        for x in range(cube.shape[2]):
            value, sigma = cube[:, y, x], error[:, y, x]
            good = np.isfinite(value) & np.isfinite(sigma) & (sigma > 0)
            if good.sum() < 100:
                output.append(dict(x=x, y=y, status="insufficient_positive_uncertainties"))
                continue
            a = design[good] / sigma[good, None]
            b = value[good] / sigma[good]
            coef, _, rank, singular = np.linalg.lstsq(a, b, rcond=None)
            if rank != design.shape[1]:
                output.append(dict(x=x, y=y, status="rank_deficient"))
                continue
            residual = b - a @ coef
            reduced_chi2 = float(np.sum(residual**2) / (len(b) - len(coef)))
            covariance = np.linalg.inv(a.T @ a) * reduced_chi2
            harmonic = coef[-6:].reshape(3, 2)
            fundamental = harmonic[0]
            cov2 = covariance[-6:-4, -6:-4]
            amp = float(np.linalg.norm(fundamental))
            formal = float(np.sqrt(fundamental @ cov2 @ fundamental) / amp) if amp else None
            output.append(dict(x=x, y=y, status="ok", samples=int(good.sum()),
                harmonic_coefficients_electrons_per_s=harmonic.tolist(),
                fundamental_semiamplitude_electrons_per_s=amp,
                white_noise_semiamplitude_error_electrons_per_s=formal,
                residual_scaled_reduced_chi2=reduced_chi2,
                condition_number=float(singular[0] / singular[-1])))
    return dict(segment_count=int(len(np.unique(labels))),
                phase_reference_btjd=float(time[0]), pixels=output)


def analyse_product(product, sector):
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/ad-astra-rrlyrae-matplotlib")
    os.environ.setdefault("XDG_CACHE_HOME", "/tmp/ad-astra-rrlyrae-cache")
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib import pyplot as plt

    path = ROOT / product["path"]
    if sha256_file(path) != product["sha256"]:
        raise ValueError("Pixel FITS hash changed")
    with validated_fits(path) as hdus:
        table = hdus[1]
        if str(table.header.get("TIMESYS", "")).upper() != "TDB":
            raise ValueError("Explicit TDB time system required")
        if table.header.get("BJDREFI") != 2457000:
            raise ValueError("Unexpected TESS time reference")
        if str(table.header.get("TIMEUNIT", "")).lower() not in ("d", "day", "days"):
            raise ValueError("Explicit day unit required")
        wcs = WCS(hdus[2].header)
        if not wcs.has_celestial:
            raise ValueError("No celestial WCS")
        time = np.asarray(table.data["TIME"], dtype=float)
        quality = np.asarray(table.data["QUALITY"])
        cube = np.array(table.data["FLUX"], dtype=float)
        error = np.array(table.data["FLUX_ERR"], dtype=float)
        good = np.isfinite(time) & (quality == 0)
        time, cube, error = time[good], cube[good], error[good]
        sort = np.argsort(time)
        time, cube, error = time[sort], cube[sort], error[sort]
        if len(time) < 100 or cube.shape[1:] != (7, 7) or np.any(np.diff(time) <= 0):
            raise ValueError("Unexpected pixel/time dimensions")
        instrument = {k: hdus[0].header.get(k, table.header.get(k))
                      for k in ("SECTOR", "CAMERA", "CCD")}
    report_path = ROOT / "data/qlp-stellar-metadata/149392671/report.json"
    catalogue = json.loads(report_path.read_text())
    year = float(Time(2457000 + np.median(time), format="jd", scale="tdb").jyear)
    positions = []
    identities = {4757558780144090368: "TIC target", 4757558814501499392: "near companion",
                  4757558810204649856: "known RRab"}
    for source in catalogue["sources"]:
        if source["source_id"] not in identities:
            continue
        years = year - source["ref_epoch"]
        ra = source["ra"] + (source.get("pmra") or 0.) * years / 3.6e6 / np.cos(np.deg2rad(source["dec"]))
        dec = source["dec"] + (source.get("pmdec") or 0.) * years / 3.6e6
        x, y = wcs.world_to_pixel_values(ra, dec)
        positions.append(dict(source_id=source["source_id"], label=identities[source["source_id"]],
                              x=float(x), y=float(y), epoch_jyear=year,
                              propagated_ra_deg=ra, propagated_dec_deg=dec))
    border_mask = np.zeros((7, 7), dtype=bool)
    border_mask[[0, -1], :] = True
    border_mask[:, [0, -1]] = True
    background = np.nanmedian(cube[:, border_mask], axis=1)
    modes = {"baseline_only": fit_pixels(time, cube, error, None),
             "border_background_nuisance": fit_pixels(time, cube, error, background)}
    rr = next(x for x in positions if x["label"] == "known RRab")
    for name, mode in modes.items():
        pixels = [p for p in mode["pixels"] if p["status"] == "ok"]
        if not pixels:
            continue
        reference = [p for p in pixels if abs(p["x"] - round(rr["x"])) <= 1
                     and abs(p["y"] - round(rr["y"])) <= 1]
        ref_coef = np.sum([p["harmonic_coefficients_electrons_per_s"][0] for p in reference], axis=0)
        ref_coef /= np.linalg.norm(ref_coef)
        for pixel in pixels:
            pixel["rr_aperture_phase_projection_electrons_per_s"] = float(
                np.dot(pixel["harmonic_coefficients_electrons_per_s"][0], ref_coef))
        peak = max(pixels, key=lambda p: p["fundamental_semiamplitude_electrons_per_s"])
        mode["maximum_amplitude_pixel"] = {k: peak[k] for k in
            ("x", "y", "fundamental_semiamplitude_electrons_per_s", "white_noise_semiamplitude_error_electrons_per_s")}
        mode["maximum_pixel_distance_to_gaia_sources_pixels"] = {
            p["label"]: float(np.hypot(peak["x"] - p["x"], peak["y"] - p["y"])) for p in positions}
        mode["nearest_pixel_to_gaia_sources"] = {
            p["label"]: min(pixels, key=lambda pixel:
                (pixel["x"] - p["x"])**2 + (pixel["y"] - p["y"])**2)
            for p in positions}
        amplitude = np.full((7, 7), np.nan)
        projection = amplitude.copy()
        for p in pixels:
            amplitude[p["y"], p["x"]] = p["fundamental_semiamplitude_electrons_per_s"]
            projection[p["y"], p["x"]] = p["rr_aperture_phase_projection_electrons_per_s"]
        fig, axes = plt.subplots(1, 3, figsize=(12, 4), constrained_layout=True)
        for ax, values, title in zip(axes, [np.nanmedian(cube, axis=0), amplitude, projection],
                ["Median FFI pixel flux", "Fixed-period fundamental amplitude", "Signed RR-aperture phase projection"]):
            artist = ax.imshow(values, origin="lower", cmap="viridis")
            for p, marker, color in zip(positions, ["+", "x", "o"], ["white", "cyan", "red"]):
                ax.scatter(p["x"], p["y"], marker=marker, edgecolors=color if marker == "o" else None,
                           facecolors="none" if marker == "o" else color, s=90, label=p["label"])
            ax.set(title=title, xlabel="Cutout column (0-based)", ylabel="Cutout row (0-based)")
            fig.colorbar(artist, ax=ax, label="electrons / second", shrink=.72)
        axes[0].legend(fontsize=7)
        fig.suptitle(f"TIC 149392671 • sector {sector} • P = {PERIOD:.9f} d • {name}")
        fig.savefig(OUT / f"s{sector:04d}-{name}.png", dpi=160)
        plt.close(fig)
    return dict(sector=sector, status="analysed", product=product, instrument=instrument,
        samples_quality_zero=int(len(time)), time_range_btjd=[float(time.min()), float(time.max())],
        period_days=PERIOD, gaia_positions=positions, models=modes,
        limitations=["Post-selection exploratory check at a known neighbour's fixed published period.",
            "Pixel maxima are descriptive and have no source-localization confidence interval.",
            "White-noise formal errors neglect temporal and interpixel correlations.",
            "Outer border can contain stellar wings; compare baseline-only and nuisance-background models.",
            "No PRF deblending or independent holdout validation; nearby sources can remain unresolved."],
        catalogue_sha256=sha256_file(report_path))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fetch-only", action="store_true")
    args = parser.parse_args()
    DATA.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    declaration = OUT / "declaration.json"
    if not declaration.exists():
        write_json(declaration, dict(declared_utc=utc_now(), sectors=[1, 2], tic=149392671,
            period_days=PERIOD, period_source="Gaia DR3 vari_rrlyrae source_id 4757558810204649856",
            source_record="data/qlp-stellar-metadata/neighbor-rrlyrae.json",
            source_record_sha256=sha256_file(ROOT / "data/qlp-stellar-metadata/neighbor-rrlyrae.json"),
            script_sha256_at_first_fetch=sha256_file(__file__),
            classification="post-selection exploratory contamination check",
            plan="Quality zero; fixed published period with three harmonics; per-orbit linear baselines and outer-border background nuisance; pixel modulation maps with proper-motion-adjusted Gaia positions; compare raw/border-background models. No period optimization or holdout flux."))
    results = []
    with ThreadPoolExecutor(max_workers=2) as pool:
        futures = {s: pool.submit(fetch, s) for s in (1, 2)}
        for sector, future in futures.items():
            try:
                products = future.result()
                row = dict(sector=sector, status="downloaded", products=products)
            except Exception as error:
                row = dict(sector=sector, status="failed", error_type=type(error).__name__, error=str(error))
            results.append(row)
            print(json.dumps(row), flush=True)
    write_json(OUT / "download-summary.json", dict(completed_utc=utc_now(), results=results))
    if not args.fetch_only:
        analysed = []
        for row in results:
            if row["status"] != "downloaded":
                continue
            for product in row["products"]:
                try:
                    result = analyse_product(product, row["sector"])
                except Exception as error:
                    result = dict(sector=row["sector"], status="analysis_failed",
                                  error_type=type(error).__name__, error=str(error))
                analysed.append(result)
        write_json(OUT / "results.json", dict(completed_utc=utc_now(),
                   source_sha256=sha256_file(__file__), results=analysed))


if __name__ == "__main__":
    main()
