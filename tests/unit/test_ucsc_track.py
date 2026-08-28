"""Tests for the UCSC Genome Browser track output."""

import os
import subprocess

import pytest

from satellome.core_functions.models.trf_model import TRModel
from satellome.core_functions.tools import ucsc_track
from satellome.core_functions.tools.ucsc_track import (
    PERIOD_CLASSES,
    build_ucsc_track,
    class_label,
    colour_for_period,
    make_bigbed,
    write_chrom_sizes,
    write_ucsc_bed,
)


def write_sat(path, rows):
    """Write a .sat file with the real column layout of TRModel.

    Numeric columns are filled with zeros rather than left blank: TRModel
    refuses to parse an empty float, so a blank fixture would be testing the
    parser rather than the track writer.
    """
    columns = TRModel.dumpable_attributes
    numeric = set(TRModel.int_attributes) | set(TRModel.float_attributes)
    lines = ["\t".join(columns)]
    for row in rows:
        values = []
        for column in columns:
            if column in row:
                values.append(str(row[column]))
            elif column in numeric:
                values.append("0")
            else:
                values.append(".")
        lines.append("\t".join(values))
    path.write_text("\n".join(lines) + "\n")
    return path


@pytest.fixture
def sat_file(tmp_path):
    return write_sat(tmp_path / "p.sat", [
        dict(trf_id=1, trf_head="chr1 some description", trf_l_ind=1000,
             trf_r_ind=1200, trf_period=3, trf_family="(AT)n", trf_pmatch=98.5),
        dict(trf_id=2, trf_head="chr1", trf_l_ind=50, trf_r_ind=90,
             trf_period=171, trf_family="alpha", trf_pmatch=85.0),
        dict(trf_id=3, trf_head="chr2", trf_l_ind=500, trf_r_ind=560,
             trf_period=20, trf_family="mini", trf_pmatch=100.0),
    ])


class TestColouring:
    def test_period_classes_cover_every_period_without_gaps(self):
        bounds = [(low, high) for low, high, _, _ in PERIOD_CLASSES]
        for (_, prev_high), (next_low, _) in zip(bounds, bounds[1:]):
            assert next_low == prev_high + 1, "a period would fall between classes"
        assert bounds[0][0] == 1

    def test_each_class_has_a_distinct_colour(self):
        colours = [c for _, _, c, _ in PERIOD_CLASSES]
        assert len(set(colours)) == len(colours)

    def test_colour_is_valid_rgb(self):
        for period in (1, 3, 8, 50, 171, 5000):
            parts = colour_for_period(period).split(",")
            assert len(parts) == 3
            assert all(0 <= int(p) <= 255 for p in parts)

    def test_microsatellite_and_satellite_differ(self):
        assert colour_for_period(3) != colour_for_period(171)
        assert "micro" in class_label(3)
        assert "satellite" in class_label(171)

    def test_nonsense_period_falls_back_rather_than_crashing(self):
        assert colour_for_period(0) == "128,128,128"
        assert class_label(-5) == "unclassified"


class TestBedOutput:
    def test_rows_carry_coordinates_name_score_and_colour(self, sat_file, tmp_path):
        out = tmp_path / "t.bed"
        assert write_ucsc_bed(str(sat_file), str(out)) == 3

        lines = out.read_text().splitlines()
        assert lines[0].startswith("track name=")
        assert 'itemRgb="On"' in lines[0], "without this the browser ignores the colours"

        fields = lines[1].split("\t")
        assert len(fields) == 9, "BED9 is required for per-item colour"
        chrom, start, end, name, score, strand = fields[:6]
        assert chrom == "chr1"
        assert (int(start), int(end)) == (50, 90)
        assert 0 <= int(score) <= 1000

    def test_description_after_a_space_is_dropped_from_the_name(self, sat_file, tmp_path):
        """UCSC matches the assembly's sequence name, not the FASTA comment."""
        out = tmp_path / "t.bed"
        write_ucsc_bed(str(sat_file), str(out))
        for line in out.read_text().splitlines()[1:]:
            assert " " not in line.split("\t")[0]

    def test_output_is_sorted_by_chrom_then_start(self, sat_file, tmp_path):
        """bedToBigBed refuses unsorted input; the browser merely tolerates it."""
        out = tmp_path / "t.bed"
        write_ucsc_bed(str(sat_file), str(out))
        rows = [l.split("\t") for l in out.read_text().splitlines()[1:]]
        keys = [(r[0], int(r[1])) for r in rows]
        assert keys == sorted(keys)

    def test_coordinates_pass_through_unshifted(self, sat_file, tmp_path):
        """.sat is already 0-based half-open, same as BED - shifting would be a bug."""
        out = tmp_path / "t.bed"
        write_ucsc_bed(str(sat_file), str(out))
        spans = {(r.split("\t")[0], int(r.split("\t")[1]), int(r.split("\t")[2]))
                 for r in out.read_text().splitlines()[1:]}
        assert ("chr1", 1000, 1200) in spans
        assert ("chr2", 500, 560) in spans

    def test_min_length_filters(self, sat_file, tmp_path):
        out = tmp_path / "t.bed"
        assert write_ucsc_bed(str(sat_file), str(out), min_length=100) == 1

    def test_track_line_can_be_omitted(self, sat_file, tmp_path):
        out = tmp_path / "t.bed"
        write_ucsc_bed(str(sat_file), str(out), with_track_line=False)
        assert not out.read_text().startswith("track")

    def test_reversed_coordinates_are_normalised(self, tmp_path):
        sat = write_sat(tmp_path / "r.sat", [
            dict(trf_id=1, trf_head="chr1", trf_l_ind=900, trf_r_ind=100,
                 trf_period=5, trf_family="x", trf_pmatch=90.0),
        ])
        out = tmp_path / "r.bed"
        write_ucsc_bed(str(sat), str(out))
        fields = out.read_text().splitlines()[1].split("\t")
        assert (int(fields[1]), int(fields[2])) == (100, 900), "BED start must precede end"

    def test_no_partial_file_left(self, sat_file, tmp_path):
        out = tmp_path / "t.bed"
        write_ucsc_bed(str(sat_file), str(out))
        assert not (tmp_path / "t.bed.partial").exists()


class TestChromSizes:
    def test_written_from_the_genome(self, tmp_path):
        fasta = tmp_path / "g.fa"
        fasta.write_text(">chr1 with description\nACGTACGTAC\n>chr2\nTTTT\n")
        out = tmp_path / "g.chrom.sizes"

        sizes = write_chrom_sizes(str(fasta), str(out))

        assert sizes == {"chr1": 10, "chr2": 4}
        assert out.read_text() == "chr1\t10\nchr2\t4\n"

    def test_names_match_the_bed(self, tmp_path, sat_file):
        """A chrom.sizes whose names differ from the BED makes bigBed fail."""
        fasta = tmp_path / "g.fa"
        fasta.write_text(">chr1 desc\n" + "A" * 2000 + "\n>chr2\n" + "C" * 1000 + "\n")
        sizes = write_chrom_sizes(str(fasta), str(tmp_path / "s.txt"))

        bed = tmp_path / "t.bed"
        write_ucsc_bed(str(sat_file), str(bed))
        used = {l.split("\t")[0] for l in bed.read_text().splitlines()[1:]}

        assert used <= set(sizes), f"{used - set(sizes)} missing from chrom.sizes"


class TestBigBed:
    def test_absent_tool_is_reported_as_optional_not_an_error(self, tmp_path, monkeypatch, caplog):
        monkeypatch.setattr(ucsc_track.shutil, "which", lambda name: None)
        with caplog.at_level("INFO"):
            assert make_bigbed("a.bed", "s.txt", "o.bb") is False
        assert "bedToBigBed not found" in caplog.text
        assert "ERROR" not in caplog.text

    def test_track_line_is_stripped_before_conversion(self, tmp_path, monkeypatch):
        """bedToBigBed rejects a file with a track line."""
        bed = tmp_path / "t.bed"
        bed.write_text('track name="x"\nchr1\t1\t2\tn\t0\t+\t1\t2\t0,0,0\n')
        seen = {}

        def fake_run(cmd, **kwargs):
            seen["content"] = open(cmd[-3]).read()
            return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

        monkeypatch.setattr(ucsc_track.shutil, "which", lambda n: "/usr/bin/bedToBigBed")
        monkeypatch.setattr(ucsc_track.subprocess, "run", fake_run)

        assert make_bigbed(str(bed), "s.txt", str(tmp_path / "o.bb")) is True
        assert "track" not in seen["content"]
        assert not (tmp_path / "t.bed.notrack").exists(), "scratch file must be cleaned up"

    def test_failure_surfaces_the_tool_output(self, tmp_path, monkeypatch, caplog):
        bed = tmp_path / "t.bed"
        bed.write_text("chr1\t1\t2\tn\t0\t+\t1\t2\t0,0,0\n")
        monkeypatch.setattr(ucsc_track.shutil, "which", lambda n: "/usr/bin/bedToBigBed")
        monkeypatch.setattr(
            ucsc_track.subprocess, "run",
            lambda cmd, **k: subprocess.CompletedProcess(cmd, 255, "", "chrom chr9 not found"),
        )

        with caplog.at_level("ERROR"):
            assert make_bigbed(str(bed), "s.txt", str(tmp_path / "o.bb")) is False
        assert "chrom chr9 not found" in caplog.text


class TestBuild:
    def test_produces_bed_and_chrom_sizes(self, tmp_path, sat_file, monkeypatch):
        fasta = tmp_path / "g.fa"
        fasta.write_text(">chr1\n" + "A" * 2000 + "\n>chr2\n" + "C" * 1000 + "\n")
        monkeypatch.setattr(ucsc_track.shutil, "which", lambda n: None)

        produced = build_ucsc_track(str(sat_file), str(fasta), str(tmp_path), "proj")

        assert set(produced) == {"bed", "chrom_sizes"}
        assert os.path.exists(produced["bed"])
        assert os.path.exists(produced["chrom_sizes"])

    def test_empty_track_is_warned_about(self, tmp_path, monkeypatch, caplog):
        sat = write_sat(tmp_path / "e.sat", [])
        fasta = tmp_path / "g.fa"
        fasta.write_text(">chr1\nACGT\n")
        monkeypatch.setattr(ucsc_track.shutil, "which", lambda n: None)

        with caplog.at_level("WARNING"):
            build_ucsc_track(str(sat), str(fasta), str(tmp_path), "proj")

        assert "UCSC track is empty" in caplog.text, (
            "an empty track must not look like a successful full one"
        )
