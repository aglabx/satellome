"""Tests for accepting compressed / .2bit genomes.

The pipeline hands a *path* to FasTAN, TRF and the Rust helpers, and those read
plain FASTA. So the contract under test is: whatever the user passes, what
comes back is something those tools can open — or a clear error, never a
silent pass-through that yields a genome with zero tandem repeats.
"""

import bz2
import gzip
import lzma
import os
import zipfile

import pytest

from satellome.core_functions.io.fasta_file import sc_iter_fasta_brute
from satellome.core_functions.tools.input_prep import (
    InputFormatError,
    ensure_plain_fasta,
    needs_conversion,
    sniff_format,
)
from tests.unit.test_twobit_file import write_twobit

GENOME = ">chr1\nACGTACGTAC\n>chr2\nTTTTCCCCGG\n"


@pytest.fixture
def plain(tmp_path):
    path = tmp_path / "g.fasta"
    path.write_text(GENOME)
    return path


class TestSniffing:
    def test_plain_fasta(self, plain):
        assert sniff_format(plain) == "fasta"
        assert needs_conversion(plain) is False

    @pytest.mark.parametrize("kind,writer", [
        ("gzip", lambda p, d: p.write_bytes(gzip.compress(d))),
        ("bzip2", lambda p, d: p.write_bytes(bz2.compress(d))),
        ("xz", lambda p, d: p.write_bytes(lzma.compress(d))),
    ])
    def test_compressed_forms(self, tmp_path, kind, writer):
        path = tmp_path / "g.bin"
        writer(path, GENOME.encode())
        assert sniff_format(path) == kind
        assert needs_conversion(path) is True

    def test_zip(self, tmp_path):
        path = tmp_path / "g.zip"
        with zipfile.ZipFile(path, "w") as z:
            z.writestr("g.fasta", GENOME)
        assert sniff_format(path) == "zip"

    def test_2bit(self, tmp_path):
        path = write_twobit(tmp_path / "g.2bit", {"chr1": "ACGT"})
        assert sniff_format(path) == "2bit"

    def test_content_wins_over_extension(self, tmp_path):
        """A gzip named .fasta is common; the bytes decide, not the name."""
        path = tmp_path / "genome.fasta"
        path.write_bytes(gzip.compress(GENOME.encode()))
        assert sniff_format(path) == "gzip"
        assert needs_conversion(path) is True

    def test_empty_file_is_an_error(self, tmp_path):
        path = tmp_path / "e.fasta"
        path.write_bytes(b"")
        with pytest.raises(InputFormatError, match="is empty"):
            sniff_format(path)

    def test_missing_file_is_an_error(self, tmp_path):
        with pytest.raises(InputFormatError, match="cannot read"):
            sniff_format(tmp_path / "nope.fa")


class TestConversion:
    def _check_content(self, path):
        got = {h[1:]: s for h, s in sc_iter_fasta_brute(str(path))}
        assert got == {"chr1": "ACGTACGTAC", "chr2": "TTTTCCCCGG"}

    def test_plain_fasta_is_passed_through_untouched(self, plain, tmp_path):
        work = tmp_path / "work"
        usable, kind = ensure_plain_fasta(plain, work)
        assert usable == str(plain)
        assert kind is None
        assert not work.exists(), "nothing should be written for plain input"

    @pytest.mark.parametrize("suffix,compress", [
        (".gz", gzip.compress),
        (".bz2", bz2.compress),
        (".xz", lzma.compress),
    ])
    def test_compressed_inputs_become_readable_fasta(self, tmp_path, suffix, compress):
        source = tmp_path / f"g.fasta{suffix}"
        source.write_bytes(compress(GENOME.encode()))

        usable, kind = ensure_plain_fasta(source, tmp_path / "work")

        assert usable != str(source)
        self._check_content(usable)
        assert not str(usable).endswith(suffix)

    def test_zip_with_one_fasta(self, tmp_path):
        source = tmp_path / "g.zip"
        with zipfile.ZipFile(source, "w") as z:
            z.writestr("genome.fasta", GENOME)

        usable, kind = ensure_plain_fasta(source, tmp_path / "work")

        assert kind == "zip"
        self._check_content(usable)

    def test_zip_with_several_fastas_refuses_to_guess(self, tmp_path):
        source = tmp_path / "g.zip"
        with zipfile.ZipFile(source, "w") as z:
            z.writestr("a.fasta", GENOME)
            z.writestr("b.fasta", GENOME)

        with pytest.raises(InputFormatError, match="holds 2 FASTA files"):
            ensure_plain_fasta(source, tmp_path / "work")

    def test_zip_without_fasta_names_what_it_found(self, tmp_path):
        source = tmp_path / "g.zip"
        with zipfile.ZipFile(source, "w") as z:
            z.writestr("readme.txt", "hello")
            z.writestr("notes.txt", "hi")

        with pytest.raises(InputFormatError, match="no FASTA inside"):
            ensure_plain_fasta(source, tmp_path / "work")

    def test_2bit_becomes_fasta(self, tmp_path):
        source = write_twobit(tmp_path / "g.2bit", {
            "chr1": "ACGTACGTAC", "chr2": "TTTTCCCCGG",
        })

        usable, kind = ensure_plain_fasta(source, tmp_path / "work")

        assert kind == "2bit"
        self._check_content(usable)

    def test_2bit_n_and_masking_survive_the_conversion(self, tmp_path):
        source = write_twobit(tmp_path / "g.2bit", {"c": "ACGTNNNNacgt"})
        usable, _ = ensure_plain_fasta(source, tmp_path / "work")
        got = {h[1:]: s for h, s in sc_iter_fasta_brute(str(usable))}
        assert got == {"c": "ACGTNNNNacgt"}

    def test_corrupt_archive_stops_the_run(self, tmp_path):
        """A silent fallback here yields a genome with no repeats and no reason."""
        source = tmp_path / "g.fasta.gz"
        source.write_bytes(gzip.compress(GENOME.encode())[:20])

        with pytest.raises(InputFormatError, match="cannot decompress"):
            ensure_plain_fasta(source, tmp_path / "work")

    def test_compressed_non_fasta_is_rejected(self, tmp_path):
        source = tmp_path / "notes.gz"
        source.write_bytes(gzip.compress(b"just some text\n"))

        with pytest.raises(InputFormatError, match="does not look like FASTA"):
            ensure_plain_fasta(source, tmp_path / "work")

    def test_no_partial_file_survives_a_failure(self, tmp_path):
        source = tmp_path / "g.fasta.gz"
        source.write_bytes(gzip.compress(GENOME.encode())[:20])
        work = tmp_path / "work"
        with pytest.raises(InputFormatError):
            ensure_plain_fasta(source, work)
        assert not list(work.glob("*.partial")), "a half-written genome must not be left"


class TestCaching:
    def test_second_call_reuses_the_conversion(self, tmp_path):
        source = tmp_path / "g.fasta.gz"
        source.write_bytes(gzip.compress(GENOME.encode()))
        work = tmp_path / "work"

        first, _ = ensure_plain_fasta(source, work)
        stamp = os.path.getmtime(first)
        os.utime(first, (stamp - 100, stamp))  # make "unchanged" detectable

        second, _ = ensure_plain_fasta(source, work)

        assert second == first
        assert os.path.getmtime(second) == stamp, "must not rewrite an up-to-date cache"

    def test_a_newer_source_invalidates_the_cache(self, tmp_path):
        source = tmp_path / "g.fasta.gz"
        source.write_bytes(gzip.compress(GENOME.encode()))
        work = tmp_path / "work"
        first, _ = ensure_plain_fasta(source, work)

        # rewrite the source with different content, stamped in the future
        source.write_bytes(gzip.compress(">chrZ\nTTTT\n".encode()))
        future = os.path.getmtime(first) + 1000
        os.utime(source, (future, future))

        again, _ = ensure_plain_fasta(source, work)

        got = {h[1:]: s for h, s in sc_iter_fasta_brute(str(again))}
        assert got == {"chrZ": "TTTT"}, "a changed input must not reuse the old genome"

    def test_force_reconverts(self, tmp_path):
        source = tmp_path / "g.fasta.gz"
        source.write_bytes(gzip.compress(GENOME.encode()))
        work = tmp_path / "work"
        first, _ = ensure_plain_fasta(source, work)
        os.utime(first, (0, 0))

        again, _ = ensure_plain_fasta(source, work, force=True)

        assert os.path.getmtime(again) > 0


class TestConvertedFileName:
    """The converted name must carry exactly one extension.

    Downstream tools derive their output names by stripping the last extension,
    and FasTAN strips a further one from the -o root it is given. A converted
    name like `hs1.fa.satellome.fasta` made FasTAN write `hs1.fa.1aln` while
    satellome looked for `hs1.fa.satellome.1aln` — killing the run *after* the
    genome-wide search had already been paid for. Caught on CHM13.
    """

    @pytest.mark.parametrize("source,expected", [
        ("hs1.fa.gz", "hs1.fasta"),
        ("hs1.fasta.gz", "hs1.fasta"),
        ("genome.fna.bz2", "genome.fasta"),
        ("genome.fa.xz", "genome.fasta"),
        ("g.2bit", "g.fasta"),
        ("x.zip", "x.fasta"),
        ("HS1.FA.GZ", "HS1.fasta"),
    ])
    def test_exactly_one_extension(self, source, expected, tmp_path):
        from satellome.core_functions.tools.input_prep import converted_path

        name = converted_path(source, tmp_path).name
        assert name == expected
        assert name.count(".") == 1, f"{name} has more than one extension"

    def test_dots_inside_the_stem_are_preserved(self, tmp_path):
        """An NCBI name keeps its dots; that case is handled at the FasTAN root."""
        from satellome.core_functions.tools.input_prep import converted_path

        name = converted_path("GCA_009914755.4_T2T-CHM13v2.0_genomic.fna.gz", tmp_path).name
        assert name == "GCA_009914755.4_T2T-CHM13v2.0_genomic.fasta"

    def test_converted_file_is_actually_written_under_that_name(self, tmp_path):
        import gzip
        from satellome.core_functions.tools.input_prep import ensure_plain_fasta

        source = tmp_path / "hs1.fa.gz"
        source.write_bytes(gzip.compress(b">chr1\nACGT\n"))

        usable, _ = ensure_plain_fasta(source, tmp_path / "work")

        assert os.path.basename(usable) == "hs1.fasta"
