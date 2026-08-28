"""Tests for the .2bit reader.

The fixtures are encoded here by hand, straight from the UCSC format
description, rather than by calling the reader's own helpers — otherwise a
misreading of the spec would be tested against itself and pass.
"""

import struct

import pytest

from satellome.core_functions.io.twobit_file import (
    SIGNATURE,
    TwoBitError,
    TwoBitFile,
    is_twobit,
    sc_iter_twobit,
    twobit_to_fasta,
)

# The packing order is part of the format: T=0, C=1, A=2, G=3, four bases per
# byte, most significant bits first.
CODE = {"T": 0, "C": 1, "A": 2, "G": 3}


def pack_dna(sequence):
    """Encode bases into 2-bit packed bytes, padding the final byte with T."""
    padded = sequence.upper().replace("N", "T")  # N is carried by N blocks
    while len(padded) % 4:
        padded += "T"
    out = bytearray()
    for i in range(0, len(padded), 4):
        byte = 0
        for base in padded[i:i + 4]:
            byte = (byte << 2) | CODE[base]
        out.append(byte)
    return bytes(out)


def n_blocks_of(sequence):
    blocks, start = [], None
    for i, base in enumerate(sequence):
        if base.upper() == "N" and start is None:
            start = i
        elif base.upper() != "N" and start is not None:
            blocks.append((start, i - start))
            start = None
    if start is not None:
        blocks.append((start, len(sequence) - start))
    return blocks


def mask_blocks_of(sequence):
    blocks, start = [], None
    for i, base in enumerate(sequence):
        lower = base.islower()
        if lower and start is None:
            start = i
        elif not lower and start is not None:
            blocks.append((start, i - start))
            start = None
    if start is not None:
        blocks.append((start, len(sequence) - start))
    return blocks


def write_twobit(path, sequences, version=0):
    """Build a minimal but spec-correct .2bit file."""
    names = list(sequences)
    header = struct.pack("<IIII", SIGNATURE, version, len(names), 0)

    index_size = sum(1 + len(n.encode()) + 4 for n in names)
    records, offsets, cursor = [], {}, len(header) + index_size

    for name in names:
        seq = sequences[name]
        nb = n_blocks_of(seq)
        mb = mask_blocks_of(seq)
        body = struct.pack("<I", len(seq))
        body += struct.pack("<I", len(nb))
        body += b"".join(struct.pack("<I", s) for s, _ in nb)
        body += b"".join(struct.pack("<I", z) for _, z in nb)
        body += struct.pack("<I", len(mb))
        body += b"".join(struct.pack("<I", s) for s, _ in mb)
        body += b"".join(struct.pack("<I", z) for _, z in mb)
        body += struct.pack("<I", 0)          # reserved
        body += pack_dna(seq)
        offsets[name] = cursor
        cursor += len(body)
        records.append(body)

    index = b"".join(
        bytes([len(n.encode())]) + n.encode() + struct.pack("<I", offsets[n])
        for n in names
    )
    path.write_bytes(header + index + b"".join(records))
    return path


@pytest.fixture
def simple(tmp_path):
    return write_twobit(tmp_path / "g.2bit", {
        "chr1": "ACGTACGTAC",
        "chr2": "TTTTCCCCGGGGAAAA",
    })


class TestDetection:
    def test_recognises_a_2bit_file(self, simple):
        assert is_twobit(simple) is True

    def test_rejects_a_fasta(self, tmp_path):
        fasta = tmp_path / "g.fa"
        fasta.write_text(">chr1\nACGT\n")
        assert is_twobit(fasta) is False

    def test_missing_file_is_not_2bit(self, tmp_path):
        assert is_twobit(tmp_path / "nope.2bit") is False

    def test_detection_is_by_content_not_extension(self, tmp_path, simple):
        renamed = tmp_path / "genome.bin"
        renamed.write_bytes(simple.read_bytes())
        assert is_twobit(renamed) is True


class TestDecoding:
    def test_sequences_round_trip(self, simple):
        with TwoBitFile(simple) as handle:
            assert handle.names == ["chr1", "chr2"]
            assert handle.read("chr1") == "ACGTACGTAC"
            assert handle.read("chr2") == "TTTTCCCCGGGGAAAA"

    def test_lengths_without_decoding(self, simple):
        with TwoBitFile(simple) as handle:
            assert handle.sequence_sizes() == {"chr1": 10, "chr2": 16}

    def test_n_blocks_are_restored(self, tmp_path):
        """The packed stream holds junk at N positions; skipping N blocks
        yields plausible-looking but wrong sequence."""
        seq = "ACGTNNNNNNACGT"
        path = write_twobit(tmp_path / "n.2bit", {"chr1": seq})
        with TwoBitFile(path) as handle:
            assert handle.read("chr1") == seq

    def test_leading_and_trailing_n_runs(self, tmp_path):
        seq = "NNNNACGTNN"
        path = write_twobit(tmp_path / "n2.2bit", {"c": seq})
        with TwoBitFile(path) as handle:
            assert handle.read("c") == seq

    def test_soft_masking_is_preserved(self, tmp_path):
        seq = "ACGTacgtACGT"
        path = write_twobit(tmp_path / "m.2bit", {"c": seq})
        with TwoBitFile(path) as handle:
            assert handle.read("c") == seq

    def test_soft_masking_can_be_dropped(self, tmp_path):
        path = write_twobit(tmp_path / "m2.2bit", {"c": "ACGTacgt"})
        with TwoBitFile(path) as handle:
            assert handle.read("c", soft_mask=False) == "ACGTACGT"

    def test_length_not_a_multiple_of_four(self, tmp_path):
        """The last byte is padded; the padding must not leak into the output."""
        for seq in ("A", "AC", "ACG", "ACGTA", "ACGTAC", "ACGTACG"):
            path = write_twobit(tmp_path / f"p{len(seq)}.2bit", {"c": seq})
            with TwoBitFile(path) as handle:
                assert handle.read("c") == seq, f"failed for length {len(seq)}"

    def test_empty_sequence(self, tmp_path):
        path = write_twobit(tmp_path / "e.2bit", {"c": ""})
        with TwoBitFile(path) as handle:
            assert handle.read("c") == ""
            assert handle.sequence_sizes() == {"c": 0}

    def test_chunked_reads_match_the_whole_sequence(self, tmp_path):
        seq = "ACGTNNAcgtTTGGCA" * 40
        path = write_twobit(tmp_path / "c.2bit", {"c": seq})
        with TwoBitFile(path) as handle:
            for chunk_size in (1, 3, 4, 5, 7, 64, 10_000):
                joined = "".join(handle.iter_chunks("c", chunk_size=chunk_size))
                assert joined == seq, f"chunk_size={chunk_size} decoded differently"

    def test_chunk_boundaries_inside_an_n_block(self, tmp_path):
        """A chunk edge falling inside a block is where off-by-one lives."""
        seq = "AC" + "N" * 10 + "GT"
        path = write_twobit(tmp_path / "cb.2bit", {"c": seq})
        with TwoBitFile(path) as handle:
            for chunk_size in range(1, len(seq) + 2):
                assert "".join(handle.iter_chunks("c", chunk_size=chunk_size)) == seq


class TestIteration:
    def test_sc_iter_twobit_matches_the_fasta_reader_shape(self, simple):
        records = list(sc_iter_twobit(simple))
        assert records == [(">chr1", "ACGTACGTAC"), (">chr2", "TTTTCCCCGGGGAAAA")]


class TestErrors:
    def test_not_a_2bit_file_is_named_as_such(self, tmp_path):
        path = tmp_path / "bad.2bit"
        path.write_bytes(b">chr1\nACGT\n")
        with pytest.raises(TwoBitError, match="not a .2bit file"):
            TwoBitFile(path)

    def test_truncated_header(self, tmp_path):
        path = tmp_path / "short.2bit"
        path.write_bytes(struct.pack("<I", SIGNATURE))
        with pytest.raises(TwoBitError, match="truncated header"):
            TwoBitFile(path)

    def test_unsupported_version_is_refused_not_misread(self, tmp_path):
        path = tmp_path / "v9.2bit"
        path.write_bytes(struct.pack("<IIII", SIGNATURE, 9, 0, 0))
        with pytest.raises(TwoBitError, match="unsupported .2bit version"):
            TwoBitFile(path)

    def test_unknown_sequence_name(self, simple):
        with TwoBitFile(simple) as handle:
            with pytest.raises(KeyError, match="no sequence named"):
                handle.read("chrX")

    def test_truncated_dna_is_reported(self, tmp_path, simple):
        path = tmp_path / "cut.2bit"
        path.write_bytes(simple.read_bytes()[:-2])
        with TwoBitFile(path) as handle:
            with pytest.raises(TwoBitError, match="truncated DNA"):
                handle.read("chr2")


class TestFastaConversion:
    def test_writes_wrapped_fasta(self, tmp_path):
        seq = "ACGT" * 25  # 100 bases
        path = write_twobit(tmp_path / "w.2bit", {"chr1": seq})
        out = tmp_path / "out.fa"

        assert twobit_to_fasta(path, out, line_width=60) == 1

        lines = out.read_text().splitlines()
        assert lines[0] == ">chr1"
        assert [len(x) for x in lines[1:]] == [60, 40]
        assert "".join(lines[1:]) == seq

    def test_conversion_preserves_content_exactly(self, tmp_path):
        sequences = {
            "chr1": "ACGTNNNNacgtACGT" * 13,
            "chr2": "N" * 7 + "ACG",
            "chr3": "A",
        }
        path = write_twobit(tmp_path / "x.2bit", sequences)
        out = tmp_path / "x.fa"

        assert twobit_to_fasta(path, out, line_width=17) == 3

        from satellome.core_functions.io.fasta_file import sc_iter_fasta_brute

        got = {h[1:]: s for h, s in sc_iter_fasta_brute(str(out))}
        assert got == sequences, "round trip through FASTA must not alter a base"

    def test_no_partial_file_is_left_behind(self, tmp_path):
        path = write_twobit(tmp_path / "a.2bit", {"c": "ACGT"})
        out = tmp_path / "a.fa"
        twobit_to_fasta(path, out)
        assert not (tmp_path / "a.fa.partial").exists()
