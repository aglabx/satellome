"""Classification is by trailing suffix, and the corpus punishes anything else."""

import pytest

from satellome.compact.policy import (
    DROPPED_CLASSES,
    KINDS,
    FileClass,
    classify_path,
    pruned_kind_names,
    strip_compression,
)


def kind_of(rel_path):
    match = classify_path(rel_path)
    return match.kind.name if match.kind else None


def test_basename_need_not_match_the_directory():
    # The live case: directory GCA_029289425.3 holds GCF_029289425.2_* files.
    # Classifying by prefix would skip the whole assembly and call it a success.
    rel = "fastan/GCF_029289425.2_NHGRI_mPanPan1-v2.0_pri_genomic.monomers.tsv.gz"
    match = classify_path(rel)
    assert match.kind.name == "monomers"
    assert match.kind.file_class == FileClass.PRIMARY_FILTERED
    assert match.compression == "gzip"


def test_a_non_accession_basename_classifies_the_same():
    assert kind_of("fastan/HLmelFor1.monomers.tsv") == "monomers"
    assert kind_of("HLmelFor1.sat.gz") == "sat_master"


@pytest.mark.parametrize(
    "rel_path,expected",
    [
        ("X.sat.gz", "sat_master"),
        ("X.1kb.sat.gz", "sat_1kb"),
        ("X.10kb.sat.gz", "sat_10kb"),
        ("X.100kb.sat.gz", "sat_100kb"),
        ("X.1000kb.sat.gz", "sat_1000kb"),
        ("X.micro.sat.gz", "sat_micro"),
        ("X.pmicro.sat.gz", "sat_pmicro"),
        ("X.tssr.sat.gz", "sat_tssr"),
        ("X.complex.sat.gz", "sat_complex"),
    ],
)
def test_longest_suffix_wins_over_the_bare_master(rel_path, expected):
    assert kind_of(rel_path) == expected


def test_compression_is_stripped_before_matching():
    assert strip_compression("X.monomers.tsv.gz") == ("X.monomers.tsv", "gzip")
    assert strip_compression("X.monomers.tsv") == ("X.monomers.tsv", "none")
    # Both forms exist in the corpus and must classify identically.
    assert kind_of("fastan/X.monomers.tsv") == kind_of("fastan/X.monomers.tsv.gz")


def test_subdirectory_matters_where_the_table_says_so():
    assert kind_of("gff3/X.micro.gff") == "gff"
    assert kind_of("fastan/X.bed") == "bed"
    # A .bed outside fastan/ is not the SINE-masking bed the policy means, so it
    # falls to the generic row rather than inheriting that one's meaning.
    assert kind_of("X.bed") == "bed_other"


def test_an_unknown_file_is_reported_not_guessed():
    assert classify_path("something_new.parquet").kind is None


def test_the_master_is_primary_and_never_dropped():
    match = classify_path("X.sat.gz")
    assert match.kind.file_class == FileClass.PRIMARY
    assert match.kind.file_class not in DROPPED_CLASSES


def test_the_10kb_view_is_kept_and_the_1kb_view_is_not():
    assert classify_path("X.10kb.sat.gz").kind.file_class == FileClass.PRIMARY
    assert classify_path("X.1kb.sat.gz").kind.file_class == FileClass.VIEW


def test_pruned_kinds_exclude_what_the_pipeline_still_reads():
    pruned = set(pruned_kind_names())
    # --large_file defaults to 1kb, and drawing/annotation read that file.
    assert "sat_1kb" not in pruned
    assert {"gff", "lengths", "aln", "sat_micro", "sat_100kb"} <= pruned


def test_every_kind_has_a_reason_recorded():
    # --explain prints these; a policy row with no justification is a row
    # nobody can argue with later.
    assert all(kind.why for kind in KINDS)


def test_every_dropped_kind_names_a_recipe():
    for kind in KINDS:
        if kind.file_class in DROPPED_CLASSES:
            assert kind.recipe, f"{kind.name} is dropped with no recipe"


@pytest.mark.parametrize(
    "rel_path,expected",
    [
        # Found unclassified on the live corpus: ~1.9 MB per assembly, ~30 GB
        # across the roster, silently left alone because no row matched.
        ("GCF_029289425.2_NHGRI_mPanPan1-v2.0_pri_genomic.gaps.bed", "gaps_bed"),
        ("reports/X.its.bed", "bed_other"),
        ("reports/telomeres.tsv", "tsv_other"),
        ("reports/sat_families.tsv", "tsv_other"),
        ("X.fna.fai", "fai"),
        ("fasta/X.arrays.fasta", "fasta_other"),
        ("reports/annotation_report.txt", "text"),
    ],
)
def test_every_kind_the_corpus_actually_contains_is_classified(rel_path, expected):
    assert kind_of(rel_path) == expected


def test_the_gaps_bed_is_kept_because_only_the_genome_could_rebuild_it():
    match = classify_path("X.gaps.bed")
    assert match.kind.file_class == FileClass.PRIMARY
    assert match.kind.file_class not in DROPPED_CLASSES


def test_beds_are_declared_headerless_rather_than_special_cased():
    for rel in ["fastan/X.bed", "X.gaps.bed", "reports/X.its.bed"]:
        assert classify_path(rel).kind.has_header is False
    assert classify_path("fastan/X.monomers.tsv").kind.has_header is True


def test_the_specific_summary_row_still_beats_the_generic_tsv_row():
    assert kind_of("reports/microsatellites.summary.tsv") == "micro_summary"
    assert kind_of("reports/anything_else.tsv") == "tsv_other"


def test_the_manifest_row_still_beats_the_generic_json_row():
    assert kind_of("run_manifest.json") == "run_manifest"
    assert kind_of("something.json") == "json"
