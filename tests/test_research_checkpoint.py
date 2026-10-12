"""Checkpoint read-back must reject missing, altered and duplicate members."""
import hashlib
import importlib.util
import io
from pathlib import Path
import tarfile

import pytest

spec = importlib.util.spec_from_file_location(
    "research_checkpoint", Path(__file__).parents[1] / "scripts/create_research_checkpoint.py")
checkpoint = importlib.util.module_from_spec(spec)
spec.loader.exec_module(checkpoint)


def archive(path, contents):
    with tarfile.open(path, "w:gz") as tar:
        for name, payload in contents:
            info = tarfile.TarInfo(name)
            info.size = len(payload)
            tar.addfile(info, io.BytesIO(payload))


@pytest.mark.parametrize("contents", [[], [("ad-astra/file", b"bad")],
                                     [("ad-astra/file", b"ok"), ("ad-astra/file", b"ok")]])
def test_readback_rejects_missing_altered_or_duplicate(tmp_path, contents):
    path = tmp_path / "archive.tar.gz"
    archive(path, contents)
    inventory = [{"archive_path": "ad-astra/file", "size_bytes": 2,
                  "sha256": hashlib.sha256(b"ok").hexdigest()}]
    with pytest.raises(ValueError):
        checkpoint.verify_archive(path, inventory)


def test_readback_accepts_exact_bytes(tmp_path):
    path = tmp_path / "archive.tar.gz"
    archive(path, [("ad-astra/file", b"ok")])
    assert checkpoint.verify_archive(path, [{"archive_path": "ad-astra/file", "size_bytes": 2,
                                           "sha256": hashlib.sha256(b"ok").hexdigest()}]) == 1
