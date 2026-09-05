"""The ``.satz`` container: one compressed stream per column, not per line.

Why columnar
------------
The panel's big tables are row-major TSV run through gzip. A row of
``monomers.tsv`` interleaves seventeen different value types — an array id, a
handful of small integers, three fixed-point numbers, a DNA string — and gzip's
32 KB window sees all of them mixed together. Split the same table by column
and each stream becomes homogeneous: ``array_id`` repeats roughly thirteen
times per array and collapses from 45.6 MB to 0.9 MB, which no window size can
achieve on the interleaved form (measured: ``--long`` windows buy 1%, the split
buys 2×).

Measured against whole-file compression at the same codec level, on six real
tables: **1.98× byte-weighted**, and the gain grows with table size, because a
bigger table gives the codec longer homogeneous runs per column while the
row-major form only interleaves more.

Exactness
---------
The container reproduces its input **byte for byte**, which is what makes
``expand`` a real inverse rather than an approximation:

* values are length-prefixed, so a field may contain anything at all;
* rows whose field count differs from the header are stored verbatim as
  *oddities* rather than being forced into the column layout;
* the comment preamble, the header line and the presence or absence of a final
  newline are all recorded.

The footer carries the md5 and byte length of the content the container
reproduces, so a reader can verify itself without consulting anything else.

Layouts
-------
``tsv``    columnar table, optionally filtered on write.
``fasta``  two columns — the header line (with its terminator) and the record
           body verbatim — so concatenation is exact and the highly repetitive
           header lines compress apart from the sequence.
``blob``   whole file, one stream. For HTML reports and logs, where there are no
           columns to split.
"""

import base64
import hashlib
import io
import json
import os
import struct

MAGIC = b"SATZ1\n"
CONTAINER_VERSION = 1

#: Rows per chunk. Large enough that a column stream is worth compressing on its
#: own, small enough that a 3.2 GB table never lands in RAM whole.
DEFAULT_CHUNK_ROWS = 100_000

#: Bytes per blob chunk in ``blob`` layout.
BLOB_CHUNK_BYTES = 8 << 20

CHUNK_MARKER = 1
END_MARKER = 0

DEFAULT_LEVEL = 12


class ContainerError(Exception):
    """The container is unusable: wrong magic, truncated, or a codec mismatch."""


def _zstd():
    """Import python-zstandard, or say exactly what to install.

    Compaction has no gzip fallback on purpose. Silently writing a weaker codec
    would produce a corpus that is half one format and half another, with the
    difference invisible until someone measured the volume again.
    """
    try:
        import zstandard
    except ImportError as e:  # pragma: no cover - exercised by the message test
        raise ContainerError(
            "satellome compact needs the 'zstandard' package for the columnar "
            "codec. Install it with: pip install zstandard"
        ) from e
    return zstandard


def _compress(data, level):
    return _zstd().ZstdCompressor(level=level).compress(data)


def _decompress(data):
    return _zstd().ZstdDecompressor().decompress(data, max_output_size=0)


def _b64(raw):
    return base64.b64encode(raw).decode("ascii")


def _unb64(text):
    return base64.b64decode(text.encode("ascii"))


def _write_blob(fh, payload, level):
    packed = _compress(payload, level)
    fh.write(struct.pack("<Q", len(packed)))
    fh.write(packed)
    return len(packed)


def _read_blob(fh):
    raw = fh.read(8)
    if len(raw) != 8:
        raise ContainerError("truncated container: expected a block length")
    (length,) = struct.unpack("<Q", raw)
    payload = fh.read(length)
    if len(payload) != length:
        raise ContainerError(
            f"truncated container: block claims {length} bytes, {len(payload)} present"
        )
    return _decompress(payload)


def _lengths_blob(values):
    return b"".join(b"%d\n" % len(v) for v in values)


def _parse_lengths(blob):
    if not blob:
        return []
    return [int(x) for x in blob.split(b"\n") if x]


def _split_values(data, lengths):
    values = []
    pos = 0
    for length in lengths:
        values.append(data[pos:pos + length])
        pos += length
    if pos != len(data):
        raise ContainerError(
            f"column data does not match its lengths ({pos} consumed of {len(data)})"
        )
    return values


class ColumnarWriter:
    """Write a ``.satz`` container, one chunk at a time.

    Not used directly by callers; :func:`write_table`, :func:`write_fasta` and
    :func:`write_blob` are the entry points.
    """

    def __init__(self, fh, header, level=DEFAULT_LEVEL, chunk_rows=DEFAULT_CHUNK_ROWS):
        self.fh = fh
        self.level = level
        self.chunk_rows = chunk_rows
        self.rows = 0
        self._digest = hashlib.md5()
        self._decoded_bytes = 0
        header = dict(header)
        header["container"] = CONTAINER_VERSION
        header["codec"] = "zstd"
        header["level"] = level
        fh.write(MAGIC)
        blob = json.dumps(header, sort_keys=True).encode("utf-8")
        fh.write(struct.pack("<I", len(blob)))
        fh.write(blob)

    def account(self, raw):
        """Record bytes the container must reproduce, for the footer digest."""
        self._digest.update(raw)
        self._decoded_bytes += len(raw)

    def write_chunk(self, columns, oddities, n_rows):
        """Write one chunk.

        Args:
            columns: list of per-column value lists (bytes), all the same length.
            oddities: list of ``(local_index, raw_line_bytes)`` for rows that do
                not fit the column layout.
            n_rows: total rows in this chunk, oddities included.
        """
        self.fh.write(struct.pack("<BII", CHUNK_MARKER, n_rows, len(oddities)))
        payload = json.dumps([[i, _b64(line)] for i, line in oddities]).encode("utf-8")
        _write_blob(self.fh, payload, self.level)
        for values in columns:
            _write_blob(self.fh, _lengths_blob(values), self.level)
            _write_blob(self.fh, b"".join(values), self.level)
        self.rows += n_rows

    def close(self, extra=None):
        self.fh.write(struct.pack("<B", END_MARKER))
        footer = {
            "rows": self.rows,
            "bytes": self._decoded_bytes,
            "md5": self._digest.hexdigest(),
        }
        if extra:
            footer.update(extra)
        blob = json.dumps(footer, sort_keys=True).encode("utf-8")
        self.fh.write(blob)
        self.fh.write(struct.pack("<I", len(blob)))
        return footer


def read_header(path):
    """Return the header dict of a container without decoding any data."""
    with open(path, "rb") as fh:
        magic = fh.read(len(MAGIC))
        if magic != MAGIC:
            raise ContainerError(
                f"{path} is not a satellome columnar container "
                f"(magic {magic!r}, expected {MAGIC!r})"
            )
        raw = fh.read(4)
        if len(raw) != 4:
            raise ContainerError(f"{path}: truncated header length")
        (length,) = struct.unpack("<I", raw)
        blob = fh.read(length)
        if len(blob) != length:
            raise ContainerError(f"{path}: truncated header")
        return json.loads(blob.decode("utf-8"))


def header_line(path):
    """The original header row of a ``tsv`` container, terminator included."""
    return _unb64(read_header(path).get("header_line") or "")


def read_footer(path):
    """Return the footer dict — rows, decoded byte length and md5."""
    size = os.path.getsize(path)
    with open(path, "rb") as fh:
        if size < 4:
            raise ContainerError(f"{path}: too short to be a container")
        fh.seek(-4, os.SEEK_END)
        (length,) = struct.unpack("<I", fh.read(4))
        if length <= 0 or length + 4 > size:
            raise ContainerError(f"{path}: footer length {length} is not plausible")
        fh.seek(-(4 + length), os.SEEK_END)
        return json.loads(fh.read(length).decode("utf-8"))


def _iter_chunks(fh, n_columns):
    """Yield ``(n_rows, oddities, columns)`` for each chunk in an open container."""
    while True:
        marker = fh.read(1)
        if not marker:
            raise ContainerError("truncated container: no end marker")
        (kind,) = struct.unpack("<B", marker)
        if kind == END_MARKER:
            return
        if kind != CHUNK_MARKER:
            raise ContainerError(f"unknown chunk marker {kind}")
        n_rows, n_odd = struct.unpack("<II", fh.read(8))
        oddities = {
            int(i): _unb64(line) for i, line in json.loads(_read_blob(fh).decode("utf-8"))
        }
        if len(oddities) != n_odd:
            raise ContainerError(
                f"chunk declares {n_odd} verbatim rows but carries {len(oddities)}"
            )
        columns = []
        for _ in range(n_columns):
            lengths = _parse_lengths(_read_blob(fh))
            columns.append(_split_values(_read_blob(fh), lengths))
        yield n_rows, oddities, columns


# --------------------------------------------------------------------------
# TSV tables
# --------------------------------------------------------------------------

def write_table(
    rows,
    out_fh,
    columns,
    preamble=b"",
    header_line=b"",
    source=None,
    filter_note=None,
    level=DEFAULT_LEVEL,
    chunk_rows=DEFAULT_CHUNK_ROWS,
    line_terminator=b"\n",
):
    """Write TSV *rows* into a columnar container.

    A row is stored column-wise only when it fits the layout exactly: the right
    number of fields *and* the file's own line terminator. Anything else — a
    ragged row, a lone CRLF in an LF file, a last line with no terminator at all
    — is kept verbatim, terminator included, so restoring is a concatenation and
    no flag has to describe the end of the file.

    Args:
        rows: iterable of ``(line_without_terminator, terminator_bytes)``.
        out_fh: binary file handle, opened for writing.
        columns: column names, from the header row.
        preamble: raw bytes preceding the header row (``#`` comment lines).
        header_line: the header row's raw bytes, terminator included.
        source: provenance dict recorded in the header.
        filter_note: description of the row filter applied, or None.

    Returns:
        The footer dict (rows, decoded bytes, md5 of the reproduced content).
    """
    header = {
        "layout": "tsv",
        "columns": list(columns),
        "preamble": _b64(preamble),
        "header_line": _b64(header_line),
        "line_terminator": _b64(line_terminator),
        "source": source or {},
        "filter": filter_note,
    }
    writer = ColumnarWriter(out_fh, header, level=level, chunk_rows=chunk_rows)
    writer.account(preamble)
    writer.account(header_line)

    n_cols = len(columns)
    buffers = [[] for _ in range(n_cols)]
    oddities = []
    in_chunk = 0

    def flush():
        nonlocal in_chunk
        if in_chunk == 0:
            return
        writer.write_chunk(buffers, oddities, in_chunk)
        for buf in buffers:
            buf.clear()
        oddities.clear()
        in_chunk = 0

    for line, terminator in rows:
        writer.account(line + terminator)
        fields = line.split(b"\t")
        if len(fields) == n_cols and terminator == line_terminator:
            for i, value in enumerate(fields):
                buffers[i].append(value)
        else:
            oddities.append((in_chunk, line + terminator))
        in_chunk += 1
        if in_chunk >= chunk_rows:
            flush()
    flush()
    return writer.close()


def read_table(path):
    """Yield raw line bytes (no terminator) for each data row of a ``tsv`` container.

    A row that was stored verbatim keeps its own terminator, which is what makes
    it verbatim; :func:`restore_table` re-emits those as they are.
    """
    header = read_header(path)
    if header.get("layout") != "tsv":
        raise ContainerError(f"{path}: layout is {header.get('layout')!r}, expected 'tsv'")
    for line, _verbatim in _iter_rows(path):
        yield line


def _iter_rows(path):
    """Yield ``(bytes, verbatim)`` per data row; verbatim rows carry a terminator."""
    header = read_header(path)
    n_cols = len(header["columns"])
    with open(path, "rb") as fh:
        _skip_header(fh)
        for n_rows, oddities, columns in _iter_chunks(fh, n_cols):
            cursor = 0
            for local in range(n_rows):
                if local in oddities:
                    yield oddities[local], True
                    continue
                yield b"\t".join(col[cursor] for col in columns), False
                cursor += 1


def _skip_header(fh):
    magic = fh.read(len(MAGIC))
    if magic != MAGIC:
        raise ContainerError(f"not a satellome columnar container (magic {magic!r})")
    (length,) = struct.unpack("<I", fh.read(4))
    fh.read(length)


def restore_table(path, out_fh):
    """Write the exact original TSV bytes of a ``tsv`` container to *out_fh*."""
    header = read_header(path)
    footer = read_footer(path)
    terminator = _unb64(header.get("line_terminator") or _b64(b"\n"))
    digest = hashlib.md5()
    written = 0

    def emit(raw):
        nonlocal written
        digest.update(raw)
        written += len(raw)
        out_fh.write(raw)

    emit(_unb64(header.get("preamble") or ""))
    emit(_unb64(header.get("header_line") or ""))

    for line, verbatim in _iter_rows(path):
        emit(line if verbatim else line + terminator)
    del footer
    return {"bytes": written, "md5": digest.hexdigest()}


# --------------------------------------------------------------------------
# FASTA
# --------------------------------------------------------------------------

def write_fasta(records, out_fh, source=None, filter_note=None, level=DEFAULT_LEVEL,
                chunk_rows=DEFAULT_CHUNK_ROWS):
    """Write FASTA *records* as two columns: header line and body, both verbatim.

    Args:
        records: iterable of ``(header_bytes, body_bytes)`` where the header
            includes its own terminator when the source had one, and the body is
            everything up to the next ``>``. Concatenating the pair reproduces
            the record exactly.
    """
    header = {
        "layout": "fasta",
        "columns": ["header", "body"],
        "source": source or {},
        "filter": filter_note,
    }
    writer = ColumnarWriter(out_fh, header, level=level, chunk_rows=chunk_rows)
    heads, bodies = [], []
    in_chunk = 0

    def flush():
        nonlocal in_chunk
        if in_chunk == 0:
            return
        writer.write_chunk([heads, bodies], [], in_chunk)
        heads.clear()
        bodies.clear()
        in_chunk = 0

    for head, body in records:
        writer.account(head)
        writer.account(body)
        heads.append(head)
        bodies.append(body)
        in_chunk += 1
        if in_chunk >= chunk_rows:
            flush()
    flush()
    return writer.close()


def read_fasta(path):
    """Yield ``(header_bytes, body_bytes)`` for each record of a ``fasta`` container."""
    header = read_header(path)
    if header.get("layout") != "fasta":
        raise ContainerError(
            f"{path}: layout is {header.get('layout')!r}, expected 'fasta'"
        )
    with open(path, "rb") as fh:
        _skip_header(fh)
        for n_rows, _oddities, columns in _iter_chunks(fh, 2):
            heads, bodies = columns
            for i in range(n_rows):
                yield heads[i], bodies[i]


def restore_fasta(path, out_fh):
    """Write the exact original FASTA bytes of a ``fasta`` container to *out_fh*."""
    digest = hashlib.md5()
    written = 0
    for head, body in read_fasta(path):
        for raw in (head, body):
            digest.update(raw)
            written += len(raw)
            out_fh.write(raw)
    return {"bytes": written, "md5": digest.hexdigest()}


# --------------------------------------------------------------------------
# Opaque blobs
# --------------------------------------------------------------------------

def write_blob_stream(in_fh, out_fh, source=None, level=DEFAULT_LEVEL):
    """Write an arbitrary file into a ``blob`` container, streaming."""
    header = {"layout": "blob", "columns": [], "source": source or {}, "filter": None}
    writer = ColumnarWriter(out_fh, header, level=level)
    while True:
        buf = in_fh.read(BLOB_CHUNK_BYTES)
        if not buf:
            break
        writer.account(buf)
        writer.write_chunk([[buf]], [], 1)
    return writer.close()


def restore_blob(path, out_fh):
    """Write the exact original bytes of a ``blob`` container to *out_fh*."""
    header = read_header(path)
    if header.get("layout") != "blob":
        raise ContainerError(
            f"{path}: layout is {header.get('layout')!r}, expected 'blob'"
        )
    digest = hashlib.md5()
    written = 0
    with open(path, "rb") as fh:
        _skip_header(fh)
        for _n_rows, _oddities, columns in _iter_chunks(fh, 1):
            raw = columns[0][0]
            digest.update(raw)
            written += len(raw)
            out_fh.write(raw)
    return {"bytes": written, "md5": digest.hexdigest()}


def restore(path, out_fh):
    """Restore any container to *out_fh*, dispatching on its recorded layout."""
    layout = read_header(path).get("layout")
    if layout == "tsv":
        return restore_table(path, out_fh)
    if layout == "fasta":
        return restore_fasta(path, out_fh)
    if layout == "blob":
        return restore_blob(path, out_fh)
    raise ContainerError(f"{path}: unknown layout {layout!r}")


def verify(path):
    """Decode a container and check it against its own footer.

    Returns ``(ok, detail)``. Used before the original file is unlinked, so a
    container that does not read back never replaces anything.
    """
    footer = read_footer(path)
    sink = _Counter()
    seen = restore(path, sink)
    if seen["md5"] != footer.get("md5"):
        return False, (
            f"decoded md5 {seen['md5']} does not match the recorded "
            f"{footer.get('md5')}"
        )
    if seen["bytes"] != footer.get("bytes"):
        return False, (
            f"decoded {seen['bytes']} bytes, footer records {footer.get('bytes')}"
        )
    return True, f"{seen['bytes']} bytes, md5 {seen['md5']}"


class _Counter(io.RawIOBase):
    """A sink that accepts writes and keeps nothing."""

    def writable(self):
        return True

    def write(self, data):
        return len(data)
