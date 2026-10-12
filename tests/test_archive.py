import io
import json
import numpy as np
import pytest
from astropy.io import fits

from adastra.archive import ArchiveClient, ones_complement_sum, sha256_file, validated_fits


def test_checksum_only_original_bytes_and_corruption(tmp_path):
    path = tmp_path / "test.fits"
    hdu = fits.PrimaryHDU(np.arange(20, dtype=np.int32))
    hdu.add_checksum(override_datasum=True)
    hdu.writeto(path)
    assert "CHECKSUM" in fits.getheader(path)
    assert "DATASUM" not in fits.getheader(path)
    with validated_fits(path) as hdus:
        np.testing.assert_array_equal(hdus[0].data, np.arange(20))
    data = bytearray(path.read_bytes())
    data[2883] ^= 1
    path.write_bytes(data)
    with pytest.raises(ValueError, match="CHECKSUM"):
        with validated_fits(path):
            pass


def test_supplied_datasum_and_header_corruption(tmp_path):
    path = tmp_path / "test.fits"
    fits.PrimaryHDU(np.arange(100, dtype=np.float64)).writeto(path, checksum=True)
    with validated_fits(path):
        pass
    with path.open("r+b") as stream:
        raw = bytearray(stream.read())
        index = raw.index(b"Written by") if b"Written by" in raw else raw.index(b"CHECKSUM") + 75
        raw[index] = ord("X") if raw[index] != ord("X") else ord("Y")
        stream.seek(0)
        stream.write(raw)
    with pytest.raises(ValueError):
        with validated_fits(path):
            pass


def test_cache_does_not_accept_altered_artifact(tmp_path):
    import hashlib
    uri = "mast:TESS/product/test.fits"
    path = tmp_path / (hashlib.sha256(uri.encode()).hexdigest() + ".fits")
    fits.PrimaryHDU(np.arange(10)).writeto(path)
    saved_hash = sha256_file(path)
    path.with_suffix(".json").write_text(json.dumps({"uri": uri, "sha256": saved_hash}))
    client = ArchiveClient(tmp_path)
    assert client.download(uri) == path
    with path.open("ab") as stream:
        stream.write(b"changed")
    with pytest.raises(ValueError, match="hash mismatch"):
        client.download(uri)


def test_checksum_carry_and_truncation():
    assert ones_complement_sum(io.BytesIO(bytes.fromhex("ffffffff00000001")), 0, 8) == 1
    with pytest.raises(ValueError, match="Truncated"):
        ones_complement_sum(io.BytesIO(b""), 0, 4)

