import hashlib
import importlib.util
from pathlib import Path
import sys


def test_multipart_archive_preserves_every_byte_and_refuses_overwrite(tmp_path):
    scripts = Path(__file__).resolve().parents[1] / 'scripts'
    sys.path.insert(0, str(scripts))
    try:
        spec = importlib.util.spec_from_file_location('incremental_checkpoint', scripts / 'create_incremental_checkpoint.py')
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        archive = tmp_path / 'science.tar.gz'
        content = bytes(range(256)) * 3 + b'end'
        archive.write_bytes(content)
        chunks = module.split_archive(archive, chunk_bytes=113)
        restored = b''.join((tmp_path / c['filename']).read_bytes() for c in chunks)
        assert restored == content
        assert [c['order'] for c in chunks] == list(range(1, len(chunks) + 1))
        assert all(c['size_bytes'] <= 113 for c in chunks)
        assert all(hashlib.sha256((tmp_path / c['filename']).read_bytes()).hexdigest() == c['sha256'] for c in chunks)
        import pytest
        with pytest.raises(FileExistsError):
            module.split_archive(archive, chunk_bytes=113)
        assert b''.join((tmp_path / c['filename']).read_bytes() for c in chunks) == content
    finally:
        sys.path.remove(str(scripts))
