"""Reading and writing the concrete file shapes a satellome run leaves behind.

This is the layer between :mod:`satellome.compact.columnar`, which knows about
columns and nothing about biology, and :mod:`satellome.compact.engine`, which
knows about policy and nothing about bytes.

Two shapes of trap live here, both found on live data rather than reasoned
about:

* ``.sat`` files begin with ``#`` comment lines before the header row, and the
  number of them is not fixed. They are preserved verbatim, not reconstructed.
* ``array_id`` is ``<chrom>_<start>_<end>_<length>_<period>`` and chromosome
  names contain underscores and dots — ``NC_073249.2``, ``ABRL03000001.1``.
  The array length is the **second-to-last** field, so it is parsed from the
  right. Splitting from the left is wrong on every RefSeq-style name.
"""

import gzip
import hashlib
import io

from satellome.compact import columnar

#: Read/write buffer for whole-file digests and blob streaming.
IO_CHUNK = 4 << 20


class FormatError(Exception):
    """A file does not have the shape its name claims."""


def open_maybe_gzip(path, compression):
    """Open *path* for reading in binary, transparently un-gzipping."""
    if compression == "gzip":
        return gzip.open(path, "rb")
    return open(path, "rb")


def md5_of_file(path):
    """md5 of the bytes on disk (the compressed bytes, for a ``.gz``)."""
    digest = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(IO_CHUNK), b""):
            digest.update(block)
    return digest.hexdigest()


#: Level used when a file is written back as gzip.
GZIP_LEVEL = 6


def gzip_framing(path):
    """The parts of a gzip header that carry information, not encoding choices.

    A gzip stream cannot be reproduced byte for byte from its content: the
    deflate output depends on which compressor and which zlib build made it, and
    the corpus was gzipped by at least one tool that is not this Python. So the
    guarantee compaction makes is about **content**, and what is preserved from
    the container is what a reader can actually observe — the embedded original
    file name and modification time.

    Returns None when *path* is not a gzip member.
    """
    try:
        with open(path, "rb") as fh:
            head = fh.read(10)
            if len(head) < 10 or head[:2] != b"\x1f\x8b":
                return None
            flags = head[3]
            mtime = int.from_bytes(head[4:8], "little")
            name = None
            if flags & 0x04:  # FEXTRA
                extra_len = int.from_bytes(fh.read(2), "little")
                fh.read(extra_len)
            if flags & 0x08:  # FNAME
                raw = bytearray()
                while True:
                    byte = fh.read(1)
                    if not byte or byte == b"\x00":
                        break
                    raw += byte
                name = raw.decode("latin-1")
            return {"mtime": mtime, "fname": name}
    except OSError:
        return None


def open_gzip_out(path, framing=None, level=GZIP_LEVEL):
    """Open *path* for writing as gzip, replaying the recorded framing."""
    framing = framing or {}
    handle = open(path, "wb")
    return gzip.GzipFile(
        filename=framing.get("fname") or "",
        mode="wb",
        compresslevel=level,
        fileobj=handle,
        mtime=framing.get("mtime"),
    )


def md5_of_content(path, compression):
    """md5 of the *decompressed* bytes, so a re-gzip does not look like a change."""
    digest = hashlib.md5()
    with open_maybe_gzip(path, compression) as fh:
        for block in iter(lambda: fh.read(IO_CHUNK), b""):
            digest.update(block)
    return digest.hexdigest()


# --------------------------------------------------------------------------
# array ids
# --------------------------------------------------------------------------

def array_length_from_id(array_id):
    """Length of the array named by ``<chrom>_<start>_<end>_<length>_<period>``.

    Parsed from the right, because chromosome names contain underscores.

    Returns:
        int, or None when the id does not have that shape — in which case the
        caller must keep the row rather than guess, since a row we cannot
        measure is not a row we may delete.
    """
    if isinstance(array_id, bytes):
        try:
            array_id = array_id.decode("utf-8", "strict")
        except UnicodeDecodeError:
            return None
    parts = array_id.rsplit("_", 2)
    if len(parts) != 3:
        return None
    try:
        return int(parts[1])
    except ValueError:
        return None


def array_length_from_fasta_header(header):
    """Length of the array named in a ``>array_id ...`` FASTA header line."""
    text = header[1:] if header.startswith(b">") else header
    token = text.split(None, 1)[0] if text.split(None, 1) else b""
    return array_length_from_id(token)


# --------------------------------------------------------------------------
# TSV tables
# --------------------------------------------------------------------------

def _split_preamble(fh, comment_prefix=b"#"):
    """Consume leading comment lines; return ``(preamble_bytes, first_data_line)``.

    ``first_data_line`` keeps its terminator so the caller can put it back
    verbatim; it is None at end of file.
    """
    preamble = bytearray()
    while True:
        line = fh.readline()
        if not line:
            return bytes(preamble), None
        if line.startswith(comment_prefix):
            preamble += line
            continue
        return bytes(preamble), line


def _iter_lines(fh):
    """Yield ``(line_without_terminator, terminator_bytes)`` for the rest of *fh*."""
    for raw in fh:
        if raw.endswith(b"\r\n"):
            yield raw[:-2], b"\r\n"
        elif raw.endswith(b"\n"):
            yield raw[:-1], b"\n"
        else:
            yield raw, b""


def _terminator_of(raw):
    """The line terminator of a raw line, or b"" when it has none."""
    if raw.endswith(b"\r\n"):
        return b"\r\n"
    if raw.endswith(b"\n"):
        return b"\n"
    return b''


def encode_tsv(
    src_path,
    compression,
    out_path,
    level,
    has_header=True,
    comment_prefix=b"#",
    row_filter=None,
    filter_note=None,
    source=None,
):
    """Encode a TSV file into a ``tsv`` container.

    Args:
        row_filter: ``callable(fields) -> bool`` applied to data rows, or None
            to keep everything. Rows whose field count differs from the header
            bypass the filter and are stored verbatim: an unparsable row is
            never silently dropped.

    Returns:
        The container footer, plus ``rows_in`` / ``rows_out`` counts.
    """
    kept = 0
    seen = 0

    with open_maybe_gzip(src_path, compression) as fh:
        preamble, first = _split_preamble(fh, comment_prefix)
        terminator = b"\n"
        if first is None:
            header_line = b""
            columns = []
            rows = iter(())
        elif has_header:
            header_line = first
            terminator = _terminator_of(first) or b"\n"
            stripped = first[: len(first) - len(_terminator_of(first))]
            columns = [c.decode("utf-8", "replace") for c in stripped.split(b"\t")]
            rows = _iter_lines(fh)
        else:
            # A headerless table (fastan/*.bed): the first line is data, and the
            # column names are positional so the container stays self-describing.
            header_line = b""
            terminator = _terminator_of(first) or b"\n"
            stripped = first[: len(first) - len(_terminator_of(first))]
            columns = [f"col{i + 1}" for i in range(len(stripped.split(b"\t")))]

            def _rows(first_line=first, handle=fh):
                mark = _terminator_of(first_line)
                yield first_line[: len(first_line) - len(mark)], mark
                for item in _iter_lines(handle):
                    yield item

            rows = _rows()

        n_cols = len(columns)

        def filtered():
            nonlocal kept, seen
            for line, mark in rows:
                seen += 1
                if row_filter is not None:
                    fields = line.split(b"\t")
                    # A row we cannot parse is never dropped by a filter: it is
                    # kept and stored verbatim, because a row whose shape we do
                    # not understand is not a row we may decide about.
                    if len(fields) == n_cols and not row_filter(fields):
                        continue
                kept += 1
                yield line, mark

        with open(out_path, "wb") as out:
            footer = columnar.write_table(
                filtered(),
                out,
                columns,
                preamble=preamble,
                header_line=header_line,
                source=source or {},
                filter_note=filter_note,
                level=level,
                line_terminator=terminator,
            )

    footer["rows_in"] = seen
    footer["rows_out"] = kept
    return footer


def decode_to_file(container_path, out_path, compression="none", framing=None,
                   level=GZIP_LEVEL):
    """Restore a container to *out_path*, re-gzipping when asked.

    Returns the digest dict of the *content* (uncompressed), so it can be
    compared with what was recorded before compaction regardless of whether the
    file was gzipped then or now.
    """
    if compression == "gzip":
        with open_gzip_out(out_path, framing, level=GZIP_LEVEL) as out:
            return columnar.restore(container_path, out)
    with open(out_path, "wb") as out:
        return columnar.restore(container_path, out)


# --------------------------------------------------------------------------
# FASTA
# --------------------------------------------------------------------------

def iter_fasta_records(fh):
    """Yield ``(header_line_with_terminator, body_bytes)`` pairs, verbatim.

    Concatenating a pair reproduces the record exactly, including the absence of
    a final newline, so no flag is needed to describe the end of the file.
    """
    header = None
    body = bytearray()
    for raw in fh:
        if raw.startswith(b">"):
            if header is not None:
                yield header, bytes(body)
            header = raw
            body = bytearray()
        else:
            if header is None:
                # Bytes before the first '>' are not a FASTA record. Emit them
                # under an empty header so the round trip stays exact rather
                # than dropping content nobody expected to be there.
                header = b""
            body += raw
    if header is not None:
        yield header, bytes(body)


def encode_fasta(src_path, compression, out_path, level, record_filter=None,
                 filter_note=None, source=None):
    """Encode a FASTA file into a ``fasta`` container."""
    seen = 0
    kept = 0

    with open_maybe_gzip(src_path, compression) as fh:
        def records():
            nonlocal seen, kept
            for header, body in iter_fasta_records(fh):
                seen += 1
                if record_filter is not None and not record_filter(header):
                    continue
                kept += 1
                yield header, body

        with open(out_path, "wb") as out:
            footer = columnar.write_fasta(
                records(), out, source=source or {},
                filter_note=filter_note, level=level,
            )

    footer["rows_in"] = seen
    footer["rows_out"] = kept
    return footer


# --------------------------------------------------------------------------
# Opaque files
# --------------------------------------------------------------------------

def encode_blob(src_path, compression, out_path, level, source=None):
    """Encode any file into a ``blob`` container."""
    with open_maybe_gzip(src_path, compression) as fh:
        with open(out_path, "wb") as out:
            return columnar.write_blob_stream(fh, out, source=source or {}, level=level)


# --------------------------------------------------------------------------
# lengths <- decomposed
# --------------------------------------------------------------------------

def lengths_from_decomposed_records(records, out_fh):
    """Write the ``.lengths`` projection of decomposed FASTA *records*.

    ``lengths`` is ``decomposed`` with every monomer sequence replaced by its
    length; the header lines are identical. Verified byte-identical on real
    files, which is the only reason ``.lengths`` may be dropped at all.
    """
    digest = hashlib.md5()
    written = 0
    for header, body in records:
        for chunk in (header, _lengths_body(body)):
            digest.update(chunk)
            written += len(chunk)
            out_fh.write(chunk)
    return {"bytes": written, "md5": digest.hexdigest()}


def _lengths_body(body):
    """Replace every whitespace-separated token of *body* with its length."""
    out = bytearray()
    for raw in body.splitlines(keepends=True):
        if raw.endswith(b"\r\n"):
            line, terminator = raw[:-2], b"\r\n"
        elif raw.endswith(b"\n"):
            line, terminator = raw[:-1], b"\n"
        else:
            line, terminator = raw, b""
        tokens = line.split(b" ")
        out += b" ".join(b"%d" % len(token) for token in tokens) + terminator
    return bytes(out)


class DigestWriter(io.RawIOBase):
    """Wrap a writable handle, keeping an md5 and a byte count as it goes."""

    def __init__(self, handle):
        self._handle = handle
        self.digest = hashlib.md5()
        self.bytes = 0

    def writable(self):
        return True

    def write(self, data):
        raw = bytes(data)
        self.digest.update(raw)
        self.bytes += len(raw)
        if self._handle is not None:
            self._handle.write(raw)
        return len(raw)

    def result(self):
        return {"bytes": self.bytes, "md5": self.digest.hexdigest()}
