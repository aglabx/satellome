"""Tests for the karyotype output prefix (the slash-in-organism-name crash).

Five genomes failed in trf_draw.py, one-to-one with the organisms whose strain
designation contains a slash: MHOM/BR/75/M2904, CCAP 1055/1, MHOM/GT/2001/U1103,
MF3/22, NIH/UT8656. The prefix was built as::

    os.path.join(output_folder, f"{taxon}.karyo")

so the slashes turned the file name into ``images/MHOM/BR/75/M2904.karyo.*`` —
directories that do not exist. Plotly's write then raised FileNotFoundError at the
last step of the pipeline, with all the data already on disk.

These tests assert the property that broke: everything the drawing step writes
must land *directly* in the folder it was given.
"""

import os

import pytest

from satellome.core_functions.trf_clusters import draw_karyotypes, karyotype_prefix


FAILING_TAXA = [
    "Leishmania_braziliensis_MHOM/BR/75/M2904",
    "Chlorella_vulgaris_CCAP_1055/1",
    "Leishmania_panamensis_MHOM/GT/2001/U1103",
    "Fusarium_MF3/22",
    "Streptomyces_NIH/UT8656",
]

# Every suffix draw_karyotypes() appends to the prefix.
KARYOTYPE_SUFFIXES = [
    ".gaps.svg",
    ".gaps.1000bp.enhanced.svg",
    ".repeats.with.gaps.enhanced.svg",
    ".repeats.nogaps.enhanced.svg",
    ".raw.svg",
    ".nosing.svg",
    ".sing.svg",
    ".raw.enhanced.svg",
    ".nosing.enchanced.svg",
    ".sing.enchanced.svg",
]


@pytest.mark.parametrize("taxon", FAILING_TAXA)
def test_every_karyotype_file_lands_in_the_output_folder(tmp_path, taxon):
    folder = str(tmp_path / "images")
    prefix = karyotype_prefix(folder, taxon)

    for suffix in KARYOTYPE_SUFFIXES:
        path = prefix + suffix
        assert os.path.dirname(path) == folder, f"{path} escapes {folder}"
        # The parent directory exists as soon as the images folder does, so the
        # write cannot fail with FileNotFoundError.
        assert os.path.dirname(path) == os.path.dirname(prefix)


@pytest.mark.parametrize("taxon", FAILING_TAXA)
def test_the_naive_prefix_would_have_escaped_the_output_folder(tmp_path, taxon):
    """Shows the defect these names triggered, so the guard above has teeth."""
    folder = str(tmp_path / "images")
    naive = os.path.join(folder, f"{taxon}.karyo")

    assert os.path.dirname(naive) != folder
    assert not os.path.isdir(os.path.dirname(naive))


def test_karyotype_files_are_actually_writable(tmp_path):
    folder = tmp_path / "images"
    folder.mkdir()
    prefix = karyotype_prefix(str(folder), "Leishmania_braziliensis_MHOM/BR/75/M2904")

    for suffix in KARYOTYPE_SUFFIXES:
        with open(prefix + suffix, "w") as fh:
            fh.write("<svg/>")

    written = sorted(p.name for p in folder.iterdir())
    assert len(written) == len(KARYOTYPE_SUFFIXES)
    assert written[0].startswith("Leishmania_braziliensis_MHOM_BR_75_M2904.karyo")
    # Nothing was created below the images folder.
    assert not any(p.is_dir() for p in folder.iterdir())


def test_a_safe_taxon_name_keeps_its_file_names(tmp_path):
    folder = str(tmp_path)
    assert karyotype_prefix(folder, "Homo_sapiens") == os.path.join(folder, "Homo_sapiens.karyo")


def test_a_renamed_taxon_is_reported_not_silently_renamed(tmp_path, caplog):
    with caplog.at_level("WARNING"):
        prefix = karyotype_prefix(str(tmp_path), "Fusarium_MF3/22")

    assert prefix.endswith("Fusarium_MF3_22.karyo")
    messages = [r.getMessage() for r in caplog.records]
    assert any("Fusarium_MF3/22" in m and "Fusarium_MF3_22" in m for m in messages), messages


def test_a_safe_taxon_name_produces_no_warning(tmp_path, caplog):
    with caplog.at_level("WARNING"):
        karyotype_prefix(str(tmp_path), "Homo_sapiens")

    assert [r.getMessage() for r in caplog.records] == []


def test_drawing_a_slashed_taxon_writes_all_ten_charts(tmp_path):
    """The whole karyotype step, as draw_all() calls it, for a slashed organism."""
    taxon = "Leishmania_braziliensis_MHOM/BR/75/M2904"
    folder = tmp_path / "images"
    folder.mkdir()

    scaffold_for_plot = {"scaffold": ["chr1", "chr2"], "end": [1_000_000, 500_000]}
    gaps_df = {
        "scaffold": ["chr1"],
        "start": [400_000],
        "end": [402_000],
        "length": [2_000],
    }
    df_trs = [
        {"chrm": "chr1", "start": 10_000, "end": 12_000, "length": 2_000, "family_name": "FAM1"},
        {"chrm": "chr2", "start": 20_000, "end": 20_500, "length": 500, "family_name": "SING"},
    ]
    repeats_with_gap = [["chr1", 399_000, 402_000, "FAM1", "aN", 3_000]]
    repeats_without_gaps = [
        {"chrm": "chr2", "start": 20_000, "end": 20_500, "length": 500, "family_name": "SING"}
    ]

    draw_karyotypes(
        karyotype_prefix(str(folder), taxon),
        taxon,
        df_trs,
        scaffold_for_plot,
        gaps_df,
        repeats_with_gap,
        repeats_without_gaps,
        enhance=10_000,
        gap_cutoff=1_000,
    )

    charts = sorted(p.name for p in folder.iterdir())
    assert len(charts) == len(KARYOTYPE_SUFFIXES), charts
    # Everything is a file in `folder` — no directory was invented from the name.
    assert all(p.is_file() for p in folder.iterdir())
    assert all(name.startswith("Leishmania_braziliensis_MHOM_BR_75_M2904.karyo") for name in charts)
    assert all(name.endswith(".html") for name in charts)
    # No leftover in-progress files: the charts are published atomically.
    assert not any(name.endswith(".partial") for name in charts)
    # The real name, slashes included, stays in the chart titles.
    assert taxon in (folder / charts[0]).read_text()


def test_a_missing_taxon_name_still_yields_a_usable_prefix(tmp_path):
    folder = str(tmp_path)
    for taxon in [None, "", "   "]:
        prefix = karyotype_prefix(folder, taxon)
        assert os.path.dirname(prefix) == folder
        assert os.path.basename(prefix).endswith(".karyo")
        assert not os.path.basename(prefix).startswith(".")
