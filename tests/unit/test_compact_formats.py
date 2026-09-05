"""Parsing array ids, and the one projection that is pure arithmetic."""

import io

import pytest

from satellome.compact import formats


@pytest.mark.parametrize(
    "array_id,expected",
    [
        # Chromosome names contain underscores and dots; parse from the right.
        (b"NC_073249.2_10_2010_2000_171", 2000),
        (b"ABRL03000001.1_0_500_500_9", 500),
        (b"CM002639.2_27_70_43_9", 43),
        (b"MU281329.1_14998_15079_81_15", 81),
        ("NC_073249.2_10_2010_2000_171", 2000),
    ],
)
def test_array_length_is_the_second_field_from_the_right(array_id, expected):
    assert formats.array_length_from_id(array_id) == expected


def test_splitting_from_the_left_would_be_wrong():
    # Guard against the obvious wrong implementation: the left-hand split of
    # NC_073249.2_10_2010_2000_171 gives "NC" and then 073249 as a "length".
    assert formats.array_length_from_id(b"NC_073249.2_10_2010_2000_171") != 73249


@pytest.mark.parametrize("bad", [b"", b"noseparators", b"a_b", b"chr1_x_y_z_w"])
def test_an_unparsable_id_returns_none_rather_than_a_guess(bad):
    # None means "keep this row": a row we cannot measure is not a row we may
    # delete.
    assert formats.array_length_from_id(bad) is None


def test_length_is_read_from_a_fasta_header_too():
    header = b">NC_073249.2_10_2010_2000_171 cut=AC orientation=rev n_monomers=4\n"
    assert formats.array_length_from_fasta_header(header) == 2000


DECOMPOSED = (
    b">A_0_43_43_9 cut=ATG orientation=rev n_monomers=4 period=10\n"
    b"CCAAGAACCC AAGGACCC ATGGGCCCCC ATGGGCCCC CATGGG\n"
    b">B_0_81_81_15 cut=AAG orientation=rev n_monomers=6 period=15\n"
    b"TCTGATGACC AAGGCTCTAATGACC\n"
)

LENGTHS = (
    b">A_0_43_43_9 cut=ATG orientation=rev n_monomers=4 period=10\n"
    b"10 8 10 9 6\n"
    b">B_0_81_81_15 cut=AAG orientation=rev n_monomers=6 period=15\n"
    b"10 15\n"
)


def test_lengths_is_decomposed_with_every_sequence_replaced_by_its_length():
    out = io.BytesIO()
    formats.lengths_from_decomposed_records(
        formats.iter_fasta_records(io.BytesIO(DECOMPOSED)), out
    )
    assert out.getvalue() == LENGTHS


def test_lengths_projection_keeps_the_header_lines_untouched():
    out = io.BytesIO()
    formats.lengths_from_decomposed_records(
        formats.iter_fasta_records(io.BytesIO(DECOMPOSED)), out
    )
    produced = [l for l in out.getvalue().split(b"\n") if l.startswith(b">")]
    original = [l for l in DECOMPOSED.split(b"\n") if l.startswith(b">")]
    assert produced == original


def test_fasta_records_concatenate_back_to_the_original():
    joined = b"".join(
        head + body for head, body in formats.iter_fasta_records(io.BytesIO(DECOMPOSED))
    )
    assert joined == DECOMPOSED


def test_bytes_before_the_first_record_are_not_dropped():
    raw = b"stray\n>A_0_10_10_1 x\nACGT\n"
    joined = b"".join(
        head + body for head, body in formats.iter_fasta_records(io.BytesIO(raw))
    )
    assert joined == raw


def test_gzip_framing_reports_the_embedded_name_and_time(tmp_path):
    import gzip

    path = tmp_path / "x.sat.gz"
    with gzip.GzipFile(filename="x.sat", mode="wb", fileobj=open(path, "wb"),
                       mtime=1787879765) as fh:
        fh.write(b"content\n")
    framing = formats.gzip_framing(str(path))
    assert framing == {"mtime": 1787879765, "fname": "x.sat"}


def test_gzip_framing_is_none_for_a_plain_file(tmp_path):
    path = tmp_path / "x.sat"
    path.write_bytes(b"content\n")
    assert formats.gzip_framing(str(path)) is None


def test_replayed_framing_survives_a_write(tmp_path):
    import gzip

    out = tmp_path / "y.sat.gz"
    with formats.open_gzip_out(str(out), {"mtime": 42, "fname": "y.sat"}) as fh:
        fh.write(b"hello\n")
    assert formats.gzip_framing(str(out)) == {"mtime": 42, "fname": "y.sat"}
    assert gzip.open(out, "rb").read() == b"hello\n"
