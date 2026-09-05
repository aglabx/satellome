"""The columnar container must reproduce its input byte for byte.

Everything else in compaction rests on that: if the container is only
approximately lossless, ``expand`` restores an approximation and nobody finds
out, because the approximation is a valid file.
"""

import gzip
import io
import os

import pytest

from satellome.compact import columnar, formats


def encode_lines(tmp_path, raw, **kwargs):
    src = tmp_path / "table.tsv"
    src.write_bytes(raw)
    out = tmp_path / "table.satz"
    footer = formats.encode_tsv(str(src), "none", str(out), level=3, **kwargs)
    return out, footer


def restored(path):
    buf = io.BytesIO()
    columnar.restore(str(path), buf)
    return buf.getvalue()


PLAIN = (
    b"# a comment\n"
    b"# another\n"
    b"array_id\ttype\tlength\n"
    b"NC_073249.2_10_2010_2000_171\tbase\t171\n"
    b"NC_073249.2_10_2010_2000_171\tbase\t170\n"
    b"ABRL03000001.1_0_500_500_9\tbase\t9\n"
)


def test_round_trip_reproduces_every_byte(tmp_path):
    out, footer = encode_lines(tmp_path, PLAIN)
    assert restored(out) == PLAIN
    assert footer["rows"] == 3
    assert footer["bytes"] == len(PLAIN)


def test_preamble_is_preserved_not_reconstructed(tmp_path):
    out, _ = encode_lines(tmp_path, PLAIN)
    header = columnar.read_header(str(out))
    assert header["columns"] == ["array_id", "type", "length"]
    assert restored(out).startswith(b"# a comment\n# another\n")


def test_file_without_a_final_newline_round_trips(tmp_path):
    raw = PLAIN[:-1]
    out, _ = encode_lines(tmp_path, raw)
    assert restored(out) == raw


def test_crlf_terminators_survive(tmp_path):
    raw = PLAIN.replace(b"\n", b"\r\n")
    out, _ = encode_lines(tmp_path, raw)
    # The container stores CRLF rows verbatim; the bytes come back unchanged.
    assert restored(out) == raw


def test_a_row_with_the_wrong_field_count_is_kept_verbatim(tmp_path):
    raw = PLAIN + b"a\tb\n" + b"x\ty\tz\n"
    out, footer = encode_lines(tmp_path, raw)
    assert restored(out) == raw
    assert footer["rows"] == 5


def test_empty_table_round_trips(tmp_path):
    raw = b"array_id\ttype\tlength\n"
    out, footer = encode_lines(tmp_path, raw)
    assert restored(out) == raw
    assert footer["rows"] == 0


def test_values_may_contain_anything_but_a_tab(tmp_path):
    raw = b"a\tb\n" + "значение\tsecond".encode("utf-8") + b"\n" + b"\t\n"
    out, _ = encode_lines(tmp_path, raw)
    assert restored(out) == raw


def test_row_filter_keeps_only_matching_rows_and_reports_counts(tmp_path):
    out, footer = encode_lines(
        tmp_path, PLAIN,
        row_filter=lambda fields: formats.array_length_from_id(fields[0]) >= 1000,
    )
    got = restored(out)
    assert b"ABRL03000001.1_0_500_500_9" not in got
    assert got.count(b"NC_073249.2_10_2010_2000_171") == 2
    assert (footer["rows_in"], footer["rows_out"]) == (3, 2)


def test_gzipped_input_decodes_to_the_same_content(tmp_path):
    src = tmp_path / "t.tsv.gz"
    with gzip.open(src, "wb") as fh:
        fh.write(PLAIN)
    out = tmp_path / "t.satz"
    formats.encode_tsv(str(src), "gzip", str(out), level=3)
    assert restored(out) == PLAIN


FASTA = (
    b">NC_073249.2_10_2010_2000_171 cut=AC n_monomers=2\n"
    b"ACGT ACGA\n"
    b">ABRL03000001.1_0_500_500_9 cut=A n_monomers=1\n"
    b"AAAA\n"
)


def test_fasta_round_trip_is_exact(tmp_path):
    src = tmp_path / "d.fasta"
    src.write_bytes(FASTA)
    out = tmp_path / "d.satz"
    formats.encode_fasta(str(src), "none", str(out), level=3)
    assert restored(out) == FASTA


def test_fasta_without_a_final_newline_round_trips(tmp_path):
    raw = FASTA[:-1]
    src = tmp_path / "d.fasta"
    src.write_bytes(raw)
    out = tmp_path / "d.satz"
    formats.encode_fasta(str(src), "none", str(out), level=3)
    assert restored(out) == raw


def test_fasta_record_filter_drops_whole_records(tmp_path):
    src = tmp_path / "d.fasta"
    src.write_bytes(FASTA)
    out = tmp_path / "d.satz"
    footer = formats.encode_fasta(
        str(src), "none", str(out), level=3,
        record_filter=lambda h: formats.array_length_from_fasta_header(h) >= 1000,
    )
    got = restored(out)
    assert b"ABRL03000001.1" not in got
    assert b"NC_073249.2" in got
    assert (footer["rows_in"], footer["rows_out"]) == (2, 1)


def test_blob_round_trip_is_exact(tmp_path):
    raw = os.urandom(200_000)
    src = tmp_path / "x.bin"
    src.write_bytes(raw)
    out = tmp_path / "x.satz"
    formats.encode_blob(str(src), "none", str(out), level=3)
    assert restored(out) == raw


def test_verify_accepts_a_good_container(tmp_path):
    out, _ = encode_lines(tmp_path, PLAIN)
    ok, detail = columnar.verify(str(out))
    assert ok, detail
    assert "md5" in detail


def test_verify_rejects_a_truncated_container(tmp_path):
    out, _ = encode_lines(tmp_path, PLAIN)
    data = out.read_bytes()
    out.write_bytes(data[: len(data) // 2])
    with pytest.raises(columnar.ContainerError):
        columnar.verify(str(out))


def test_a_non_container_is_rejected_by_name(tmp_path):
    bogus = tmp_path / "not.satz"
    bogus.write_bytes(b"just some text, definitely not a container\n")
    with pytest.raises(columnar.ContainerError) as excinfo:
        columnar.read_header(str(bogus))
    assert "not a satellome columnar container" in str(excinfo.value)


def test_footer_records_the_content_digest(tmp_path):
    import hashlib

    out, footer = encode_lines(tmp_path, PLAIN)
    assert footer["md5"] == hashlib.md5(PLAIN).hexdigest()
    assert columnar.read_footer(str(out))["md5"] == footer["md5"]


def test_chunking_does_not_change_the_result(tmp_path):
    raw = b"a\tb\n" + b"".join(b"row%d\tvalue%d\n" % (i, i) for i in range(2500))
    src = tmp_path / "big.tsv"
    src.write_bytes(raw)
    outputs = []
    for chunk in (7, 100, 10_000):
        out = tmp_path / f"big{chunk}.satz"
        with open(out, "wb") as fh:
            columnar.write_table(
                ((line, b"\n") for line in raw.split(b"\n")[1:] if line),
                fh, ["a", "b"], header_line=b"a\tb\n", level=3, chunk_rows=chunk,
            )
        outputs.append(restored(out))
    assert outputs[0] == outputs[1] == outputs[2]
    assert outputs[0] == raw
