"""Fast column-projected reads of large parquet files on S3 through an HTTPS proxy.

pyarrow's native S3 client is slow through the session proxy, and fsspec's block cache re-downloads
data on parquet's scattered access pattern. Here the footer is fetched once, the byte ranges of just
the requested column chunks are computed, merged and fetched in parallel with plain HTTP range
requests, and pyarrow then reads from an in-memory sparse file.
"""
import io
import concurrent.futures as cf

import pyarrow.parquet as pq

import s3io


class SparseFile(io.RawIOBase):
    """Seekable read-only file backed by a dict of fetched byte ranges (fetches misses on demand)."""

    def __init__(self, url, size, ranges=None):
        self.url, self.size, self.pos = url, size, 0
        self.chunks = []  # (start, bytes)
        if ranges:
            for a, b in ranges.items():
                self.chunks.append((a, b))

    def readable(self):
        return True

    def seekable(self):
        return True

    def tell(self):
        return self.pos

    def seek(self, off, whence=0):
        self.pos = off if whence == 0 else (self.pos + off if whence == 1 else self.size + off)
        return self.pos

    def _get(self, a, n):
        for s, data in self.chunks:
            if s <= a and a + n <= s + len(data):
                return data[a - s:a - s + n]
        data = _range(self.url, a, a + n - 1)
        self.chunks.append((a, data))
        return data

    def read(self, n=-1):
        if n is None or n < 0:
            n = self.size - self.pos
        n = max(0, min(n, self.size - self.pos))
        out = self._get(self.pos, n)
        self.pos += n
        return out

    def readinto(self, b):
        data = self.read(len(b))
        b[:len(data)] = data
        return len(data)


def _range(url, a, b):
    for attempt in range(5):
        try:
            r = s3io.session().get(url, headers={"Range": f"bytes={a}-{b}"}, timeout=120)
            r.raise_for_status()
            return r.content
        except Exception:
            if attempt == 4:
                raise
    return b""


def file_size(url):
    r = s3io.session().head(url, timeout=60)
    r.raise_for_status()
    return int(r.headers["Content-Length"])


def read_columns(key, columns, threads=16, gap=1 << 20):
    url = key if key.startswith("http") else f"{s3io.BUCKET}/{key}"
    size = file_size(url)
    tail = min(size, 1 << 20)
    sf = SparseFile(url, size, {size - tail: _range(url, size - tail, size - 1)})
    md = pq.ParquetFile(sf).metadata
    names = [md.schema.column(i).name for i in range(md.num_columns)]
    want = [i for i, n in enumerate(names) if n in set(columns)]
    spans = []
    for rg in range(md.num_row_groups):
        g = md.row_group(rg)
        for i in want:
            c = g.column(i)
            start = c.dictionary_page_offset if c.has_dictionary_page and c.dictionary_page_offset else c.data_page_offset
            spans.append((start, start + c.total_compressed_size))
    spans.sort()
    merged = []
    for a, b in spans:
        if merged and a - merged[-1][1] <= gap:
            merged[-1][1] = max(merged[-1][1], b)
        else:
            merged.append([a, b])
    # split very large spans so they download in parallel, then re-join each span into one blob
    pieces = []
    for j, (a, b) in enumerate(merged):
        while b - a > 8 << 20:
            pieces.append((j, a, a + (8 << 20)))
            a += 8 << 20
        pieces.append((j, a, b))
    with cf.ThreadPoolExecutor(threads) as ex:
        datas = list(ex.map(lambda p: _range(url, p[1], p[2] - 1), pieces))
    blobs = {}
    for (j, a, _), d in zip(pieces, datas):
        blobs.setdefault(j, []).append(d)
    for j, (a, b) in enumerate(merged):
        sf.chunks.append((a, b"".join(blobs[j])))
    sf.seek(0)
    return pq.ParquetFile(sf).read(columns=list(columns))
