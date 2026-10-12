"""Public archive access with explicit provenance and original-byte integrity.

Reconstructed after the cloud runtime lost the earlier uncommitted workspace.
This is new code, not a recovered historical implementation.
"""
from contextlib import contextmanager
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

from astropy.io import fits
import numpy as np
import requests

MAST_INVOKE = "https://mast.stsci.edu/api/v0/invoke"
MAST_DOWNLOAD = "https://mast.stsci.edu/api/v0.1/Download/file"


def utc_now():
    return datetime.now(timezone.utc).isoformat()


def sha256_file(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path, record):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".part")
    temporary.write_text(json.dumps(record, indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


def ones_complement_sum(stream, offset, length):
    """FITS-standard 32-bit one's-complement sum over original padded bytes."""
    if length % 4:
        raise ValueError("FITS byte span is not divisible by four")
    stream.seek(offset)
    total = 0
    remaining = length
    while remaining:
        block = stream.read(min(1024 * 1024, remaining))
        if not block or len(block) % 4:
            raise ValueError("Truncated FITS byte span")
        remaining -= len(block)
        total += int(np.frombuffer(block, dtype=">u4").sum(dtype=np.uint64))
        while total >> 32:
            total = (total & 0xFFFFFFFF) + (total >> 32)
    return total


@contextmanager
def validated_fits(path):
    """Verify structure and every supplied CHECKSUM/DATASUM before yielding.

    CHECKSUM-only HDUs are summed directly: no assumed DATASUM=0 and no header
    rewrite. Missing checksum cards are explicitly not an integrity guarantee.
    """
    path = Path(path)
    if path.stat().st_size % 2880:
        raise ValueError("FITS file is not a complete padded record stream")
    with fits.open(path, mode="readonly", memmap=False, checksum=False) as hdus:
        hdus.verify("exception")
        with path.open("rb") as raw:
            for hdu in hdus:
                info = hdu.fileinfo()
                if info is None:
                    raise ValueError("Missing original HDU byte locations")
                start, data, span = info["hdrLoc"], info["datLoc"], info["datSpan"]
                if data + span > path.stat().st_size:
                    raise ValueError("FITS HDU extends beyond the downloaded file")
                if "CHECKSUM" in hdu.header:
                    if ones_complement_sum(raw, start, data + span - start) != 0xFFFFFFFF:
                        raise ValueError("FITS CHECKSUM mismatch")
                if "DATASUM" in hdu.header:
                    if ones_complement_sum(raw, data, span) != int(hdu.header["DATASUM"]):
                        raise ValueError("FITS DATASUM mismatch")
        yield hdus


def query_mast(request, timeout=60):
    response = requests.get(MAST_INVOKE, params={"request": json.dumps(request)}, timeout=timeout)
    response.raise_for_status()
    payload = response.json()
    if payload.get("status") not in (None, "COMPLETE"):
        raise RuntimeError(f"MAST query status: {payload.get('status')}; result unchecked")
    return payload


class ArchiveClient:
    def __init__(self, cache):
        self.cache = Path(cache)
        self.cache.mkdir(parents=True, exist_ok=True)

    def download(self, uri):
        if not isinstance(uri, str) or not uri.startswith("mast:"):
            raise ValueError("Only explicit public MAST product URIs are accepted")
        key = hashlib.sha256(uri.encode()).hexdigest()
        path = self.cache / f"{key}.fits"
        sidecar = path.with_suffix(".json")
        if path.exists():
            if not sidecar.exists():
                raise ValueError("Existing file lacks provenance; refusing silent reuse")
            saved = json.loads(sidecar.read_text())
            if saved.get("uri") != uri or saved.get("sha256") != sha256_file(path):
                raise ValueError("Existing cache identity/hash mismatch")
            with validated_fits(path):
                pass
            return path
        temporary = path.with_suffix(".fits.part")
        if temporary.exists():
            raise ValueError("Partial artifact already exists; preserve it for explicit recovery")
        try:
            with requests.get(MAST_DOWNLOAD, params={"uri": uri}, timeout=(30, 120), stream=True) as response:
                response.raise_for_status()
                size = 0
                with temporary.open("xb") as stream:
                    for block in response.iter_content(1024 * 1024):
                        size += len(block)
                        if size > 512 * 1024 * 1024:
                            raise ValueError("Artifact exceeds the explicit 512 MiB download limit")
                        stream.write(block)
                with validated_fits(temporary) as hdus:
                    integrity = [{"hdu": i, "checksum_supplied": "CHECKSUM" in h.header,
                                  "datasum_supplied": "DATASUM" in h.header}
                                 for i, h in enumerate(hdus)]
                record = {"uri": uri, "retrieved_at_utc": utc_now(), "http_status": response.status_code,
                          "resolved_url": response.url, "bytes": size,
                          "sha256": sha256_file(temporary), "integrity": integrity,
                          "integrity_note": "Structural and supplied FITS checksum verification; SHA256 establishes subsequent cache identity, not an external publisher checksum."}
            temporary.replace(path)
            write_json(sidecar, record)
        except Exception:
            if temporary.exists():
                temporary.unlink()
            raise
        return path

