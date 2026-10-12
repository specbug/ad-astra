"""Conservative reader for the legacy (pre-sector-94) QLP light-curve schema.

This module intentionally does not accept every QLP release.  It preserves the
measured samples; filtering for a particular analysis is an explicit operation.
"""
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .archive import validated_fits


@dataclass
class QLPLightCurve:
    time_bjd: np.ndarray
    flux: np.ndarray
    flux_error: np.ndarray | None
    quality: np.ndarray
    metadata: dict

    @property
    def time_btjd(self):
        return self.time_bjd - 2457000.0

    def good_mask(self):
        """Quality-zero, finite, positive flux; no clipping or transit removal."""
        mask = (self.quality == 0) & np.isfinite(self.time_bjd)
        mask &= np.isfinite(self.flux) & (self.flux > 0)
        if self.flux_error is not None:
            mask &= np.isfinite(self.flux_error) & (self.flux_error > 0)
        return mask


def _value(headers, key, default=None):
    values = [h[key] for h in headers if key in h]
    if not values:
        return default
    if any(v != values[0] for v in values[1:]):
        raise ValueError(f"Conflicting {key} values between FITS headers")
    return values[0]


def read_qlp(path, flux_column="KSPSAP_FLUX"):
    """Read only sectors 1--73 or 80--93, checking time conventions explicitly.

    SAP_FLUX has no formal error column in the inspected legacy files.  Its
    errors remain None; KSPSAP errors are never assigned to a SAP series.
    """
    if flux_column not in {"SAP_FLUX", "KSPSAP_FLUX"}:
        raise ValueError("Only SAP_FLUX and KSPSAP_FLUX are supported")
    with validated_fits(Path(path)) as hdus:
        matches = [h for h in hdus if getattr(h, "columns", None) is not None
                   and "TIME" in h.columns.names and flux_column in h.columns.names]
        if len(matches) != 1:
            raise ValueError("Expected exactly one QLP light-curve table")
        table = matches[0]
        headers = [hdus[0].header, table.header]
        sector = int(_value(headers, "SECTOR", -1))
        if not (1 <= sector <= 73 or 80 <= sector <= 93):
            raise ValueError(f"QLP sector {sector} is outside the validated legacy scope")
        if str(_value(headers, "TIMESYS", "")).strip().upper() != "TDB":
            raise ValueError("Expected TIMESYS=TDB")
        if str(_value(headers, "TIMEUNIT", "")).strip().lower() not in {"d", "day", "days"}:
            raise ValueError("Expected an explicit day-valued TIMEUNIT")
        integer = _value(headers, "BJDREFI")
        frac_r = _value(headers, "BJDREFR")
        frac_f = _value(headers, "BJDREFF")
        if frac_r is not None and frac_f is not None and float(frac_r) != float(frac_f):
            raise ValueError("Conflicting BJDREFR and BJDREFF fractional references")
        frac = frac_r if frac_r is not None else frac_f
        if integer is None or frac is None:
            raise ValueError("Missing explicit BJD reference components")
        reference = float(integer) + float(frac)
        if not np.isfinite(reference) or not 2_000_000 < reference < 3_000_000:
            raise ValueError("Implausible BJD reference")
        required = {"TIME", flux_column, "QUALITY"}
        if flux_column == "KSPSAP_FLUX":
            required.add("KSPSAP_FLUX_ERR")
        missing = required - set(table.columns.names)
        if missing:
            raise ValueError(f"Missing QLP columns: {sorted(missing)}")
        error_column = f"{flux_column}_ERR"
        error = (np.array(table.data[error_column], dtype=float, copy=True)
                 if error_column in table.columns.names else None)
        return QLPLightCurve(
            time_bjd=np.array(table.data["TIME"], dtype=float, copy=True) + reference,
            flux=np.array(table.data[flux_column], dtype=float, copy=True),
            flux_error=error,
            quality=np.array(table.data["QUALITY"], dtype=np.int64, copy=True),
            metadata={
                "path": str(Path(path)), "sector": sector,
                "tic_id": int(_value(headers, "TICID", -1)),
                "flux_column": flux_column, "error_column": error_column if error is not None else None,
                "bjd_reference": reference, "time_scale": "TDB", "samples": len(table.data),
                "checksum_cards": ["CHECKSUM" in h.header for h in hdus],
                "datasum_cards": ["DATASUM" in h.header for h in hdus],
                "scope": "legacy QLP sectors 1-73 or 80-93",
            },
        )


def read_qlp_lightcurve(path, flux_column="KSPSAP_FLUX"):
    """Mapping facade; ``time`` is explicitly BTJD (BJD TDB minus 2457000)."""
    lc = read_qlp(path, flux_column)
    return {
        "time": lc.time_btjd, "time_bjd": lc.time_bjd,
        "flux": lc.flux, "flux_err": lc.flux_error, "quality": lc.quality,
        "tic": lc.metadata["tic_id"], "sector": lc.metadata["sector"],
        "bjdref": 2457000.0, "source_bjdref": lc.metadata["bjd_reference"],
        "time_system": "TDB", "time_format": "BTJD", "metadata": lc.metadata,
    }
