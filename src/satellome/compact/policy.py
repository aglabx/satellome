"""What each kind of satellome output is, and what compaction does with it.

The table below is the artifact; the engine only applies it. Adding a new output
kind means adding a row, not a branch.

Matching is on the **trailing suffix**, never on a prefix. The basename of every
file in a run is derived from the input FASTA, not from the accession or the
directory holding it — live case: directory ``GCA_029289425.3`` containing
``GCF_029289425.2_NHGRI_mPanPan1-v2.0_pri_genomic.*``. A run that classified by
prefix would skip that directory entirely and call it success.

Classes
-------
``PRIMARY``           kept; re-encoded columnar. Cannot be recomputed.
``PRIMARY_FILTERED``  kept, minus per-copy rows belonging to arrays below the
                      threshold; re-encoded columnar.
``PRIMARY_BLOB``      kept; whole-file re-encoded. No columns to split.
``PROJECTION``        dropped; a bit-exact function of a file that stays.
``VIEW``              dropped; recomputed by re-running the classifier over the
                      master table.
``NOT_NEEDED``        dropped by decision, not by derivability.
``KEEP``              left exactly as it is — the run's own record of itself.
``UNKNOWN``           not in the table. Kept and *named in the report*: a file
                      nobody classified must not disappear because nobody
                      classified it.
"""

import os
from typing import NamedTuple, Optional

#: Bumped whenever an action in the table changes. Recorded in every
#: ``.compact.json`` so a directory says which policy produced it.
POLICY_VERSION = 1


class FileClass:
    PRIMARY = "primary"
    PRIMARY_FILTERED = "primary_filtered"
    PRIMARY_BLOB = "primary_blob"
    PROJECTION = "projection"
    VIEW = "view"
    NOT_NEEDED = "not_needed"
    KEEP = "keep"
    UNKNOWN = "unknown"


#: Classes whose files disappear from the directory.
DROPPED_CLASSES = frozenset(
    {FileClass.PROJECTION, FileClass.VIEW, FileClass.NOT_NEEDED}
)

#: Classes whose files are rewritten into a container.
REENCODED_CLASSES = frozenset(
    {FileClass.PRIMARY, FileClass.PRIMARY_FILTERED, FileClass.PRIMARY_BLOB}
)


class Kind(NamedTuple):
    """One row of the policy table."""

    #: Trailing suffix that identifies the kind, e.g. ``.monomers.tsv``.
    suffix: str
    #: Stable name used in reports, recipes and the ledger.
    name: str
    #: One of :class:`FileClass`.
    file_class: str
    #: Container layout to write: ``tsv``, ``fasta``, ``blob`` or None.
    layout: Optional[str]
    #: Recipe kind used to put the file back, or None when it is kept.
    recipe: Optional[str]
    #: Directory the kind lives in relative to the run root, or None for any.
    subdir: Optional[str] = None
    #: True when the per-copy row filter applies.
    filtered: bool = False
    #: False for a table with no header row (every BED the pipeline writes).
    has_header: bool = True
    #: Human-readable justification, printed by ``--explain``.
    why: str = ""


#: The ``.gz`` variants are handled by stripping the suffix before matching:
#: compression across the corpus is inconsistent (14,606 assemblies carry
#: ``.tsv.gz`` for monomers, 1,012 carry the same kinds uncompressed), and both
#: forms mean the same thing.
CONTAINER_SUFFIXES = (".gz",)

#: Extension written for a re-encoded file.
SATZ_SUFFIX = ".satz"

KINDS = (
    # --- per-copy tables: the bulk of the corpus -------------------------
    Kind(".monomers.tsv", "monomers", FileClass.PRIMARY_FILTERED, "tsv",
         "decompose_monomers", subdir="fastan", filtered=True,
         why="one row per monomer copy; copies of sub-threshold arrays are "
             "recomputable from the master's trf_array"),
    Kind(".hors.tsv", "hors", FileClass.PRIMARY_FILTERED, "tsv",
         "decompose_hors", subdir="fastan", filtered=True,
         why="same schema as monomers minus parent_idx (16 columns, not 17)"),
    Kind(".decomposed.fasta", "decomposed", FileClass.PRIMARY_FILTERED, "fasta",
         "decompose_fasta", subdir="fastan", filtered=True,
         why="per-copy sequences; header lines are highly repetitive and "
             "compress apart from the sequence"),
    Kind(".lengths", "lengths", FileClass.PROJECTION, None,
         "lengths_from_decomposed", subdir="fastan",
         why="decomposed with every monomer sequence replaced by its length; "
             "verified byte-identical regeneration on 8 assemblies"),
    # --- per-array tables -------------------------------------------------
    Kind(".summary.tsv", "summary", FileClass.PRIMARY, "tsv", None, subdir="fastan",
         why="24 columns, one row per array; consensus and HOR structure that "
             "nothing else carries"),
    Kind(".bed", "bed", FileClass.PRIMARY, "tsv", None, subdir="fastan",
         has_header=False,
         why="the SINE scan masks with it"),
    Kind(".gaps.bed", "gaps_bed", FileClass.PRIMARY, "tsv", None, has_header=False,
         why="assembly gaps; found by scanning the genome, so nothing in the "
             "output directory can recompute it"),
    Kind(".bed", "bed_other", FileClass.PRIMARY, "tsv", None, has_header=False,
         why="any other BED the pipeline writes (reports/*.its.bed and the "
             "like): kept, and compressed rather than left alone"),
    Kind(".1aln", "aln", FileClass.NOT_NEEDED, None, "rerun_fastan", subdir="fastan",
         why="FasTAN's alignment intermediate; dropped by decision, and it is "
             "the one kind expand cannot rebuild without the genome"),
    # --- the master and its views ----------------------------------------
    Kind(".10kb.sat", "sat_10kb", FileClass.PRIMARY, "tsv", None,
         why="working view kept by decision, not by derivability"),
    Kind(".1kb.sat", "sat_1kb", FileClass.VIEW, None, "classify_from_master",
         why="trf_array_length >= 1000 over the master"),
    Kind(".100kb.sat", "sat_100kb", FileClass.VIEW, None, "classify_from_master",
         why="trf_array_length >= 100000 over the master"),
    Kind(".1000kb.sat", "sat_1000kb", FileClass.VIEW, None, "classify_from_master",
         why="trf_array_length >= 1000000 over the master"),
    Kind(".micro.sat", "sat_micro", FileClass.VIEW, None, "classify_from_master",
         why="every row matches a master row on all substantive columns; the "
             "only addition is trf_family = fSSR_*, declared dead"),
    Kind(".pmicro.sat", "sat_pmicro", FileClass.VIEW, None, "classify_from_master",
         why="adds trf_family = the motif, which equals (trf_consensus)n for "
             "3,733 of 3,733 rows checked"),
    Kind(".tssr.sat", "sat_tssr", FileClass.VIEW, None, "classify_from_master",
         why="adds trf_family = tSSR_*"),
    Kind(".complex.sat", "sat_complex", FileClass.VIEW, None, "classify_from_master",
         why="adds nothing at all to the master rows it selects"),
    Kind(".sat", "sat_master", FileClass.PRIMARY, "tsv", None,
         why="carries every array's complete sequence in trf_array, which is "
             "what makes the per-copy layer recomputable without the genome"),
    # --- annotation -------------------------------------------------------
    Kind(".gff", "gff", FileClass.NOT_NEEDED, None, "classify_from_master",
         subdir="gff3",
         why="dropped by decision; fSSR is an early-version atavism"),
    # --- kept as-is -------------------------------------------------------
    Kind(".tsv", "tsv_other", FileClass.PRIMARY_BLOB, "blob", None,
         why="reports/telomeres.tsv, reports/sat_families.tsv and any other "
             "small table; encoded whole because their header conventions "
             "differ and a wrong guess about a header is a wrong file"),
    Kind(".fasta", "fasta_other", FileClass.PRIMARY_BLOB, "blob", None,
         why="leftover extracted-array FASTA from the earlier campaign; kept "
             "because dropping a file on the assumption it is a leftover is "
             "exactly the assumption that loses data"),
    Kind(".txt", "text", FileClass.PRIMARY_BLOB, "blob", None,
         why="annotation reports and other free text"),
    Kind(".fai", "fai", FileClass.KEEP, None, None,
         why="a samtools index: tiny, and it has to stay readable beside the "
             "file it indexes"),
    Kind(".json", "json", FileClass.KEEP, None, None,
         why="small records that must stay readable without satellome"),
    Kind(".report", "report", FileClass.PRIMARY_BLOB, "blob", None, subdir="reports",
         why="per-family statistics; small, but 15,620 of them"),
    Kind("microsatellites.summary.tsv", "micro_summary", FileClass.PRIMARY_BLOB,
         "blob", None, subdir="reports",
         why="per-motif counts; kept rather than derived because the file's row "
             "order among ties depends on set iteration and is not reproducible"),
    Kind(".html", "html", FileClass.PRIMARY_BLOB, "blob", None,
         why="the run's report page"),
    Kind(".svg", "svg", FileClass.PRIMARY_BLOB, "blob", None,
         why="karyotype and distribution figures"),
    Kind(".png", "png", FileClass.PRIMARY, None, None,
         why="already compressed; re-encoding buys nothing"),
    Kind(".log", "log", FileClass.PRIMARY_BLOB, "blob", None,
         why="the run's own account of itself"),
    Kind("results.yaml", "results_yaml", FileClass.KEEP, None, None,
         why="the run's statistics; read by the panel's aggregation directly"),
    Kind("run_manifest.json", "run_manifest", FileClass.KEEP, None, None,
         why="the record --verify-run checks against"),
    Kind(".compact.json", "compact_record", FileClass.KEEP, None, None,
         why="this compaction's own record"),
    Kind(".rc_complete", "rc_complete", FileClass.KEEP, None, None,
         why="the driver's marker"),
)

#: Longest suffix wins, so ``.10kb.sat`` is matched before ``.sat`` and
#: ``.monomers.tsv`` before a bare ``.tsv``.
_ORDERED_KINDS = tuple(sorted(KINDS, key=lambda k: len(k.suffix), reverse=True))


class Match(NamedTuple):
    kind: Optional[Kind]
    #: Relative path as given.
    rel_path: str
    #: The file's basename with the compression suffix removed.
    stem: str
    #: "gzip" or "none".
    compression: str


def strip_compression(name):
    """Return ``(name_without_compression_suffix, "gzip" | "none")``."""
    for suffix in CONTAINER_SUFFIXES:
        if name.endswith(suffix):
            return name[: -len(suffix)], "gzip"
    return name, "none"


def classify_path(rel_path):
    """Classify one file of a run directory by its trailing suffix.

    Args:
        rel_path: path relative to the run directory, e.g.
            ``fastan/GCF_029289425.2_..._genomic.monomers.tsv.gz``.

    Returns:
        Match: with ``kind`` None when nothing in the table applies.
    """
    rel_path = rel_path.replace(os.sep, "/") if os.sep != "/" else rel_path
    parts = rel_path.split("/")
    name = parts[-1]
    subdir = parts[-2] if len(parts) > 1 else None
    stem, compression = strip_compression(name)

    for kind in _ORDERED_KINDS:
        if not stem.endswith(kind.suffix):
            continue
        if kind.subdir is not None and subdir != kind.subdir:
            continue
        return Match(kind, rel_path, stem, compression)
    return Match(None, rel_path, stem, compression)


def prefix_of(stem, kind):
    """The run's basename prefix, i.e. *stem* with the kind's suffix removed."""
    return stem[: -len(kind.suffix)] if kind.suffix.startswith(".") else stem


class CompactConfig(NamedTuple):
    """Everything about a compaction that is a decision rather than a fact.

    ``min_array_length`` is not derived from storage economics and must not be
    presented as if it were: the fraction it frees ranges from 19.0% to 91.9%
    across assemblies, so no storage argument can fix a value that moves by 5x.
    It comes from scope — the tool is for satellites, not for tandem repeats in
    general — and it matches the boundary at which ``1kb.sat`` was already cut.
    Revisit it when the panel's scientific scope changes, not when the volume
    fills up.
    """

    min_array_length: int = 1000
    level: int = 12
    verify_drops: bool = True
    keep_unknown: bool = True


#: Kinds a normal run still needs on disk when it finishes, even though
#: compaction later removes them. ``.1kb.sat`` is the drawing and annotation
#: steps' input (``--large_file`` defaults to ``1kb``) and what ``--rerun
#: drawing`` reads, so a run that deleted it would break its own amendment path.
PIPELINE_KEEPS_ANYWAY = frozenset({"sat_1kb"})


def pruned_kind_names():
    """Kinds a run without ``--extended-output`` does not leave behind.

    Derived from the table rather than listed again, so the pipeline and the
    compactor cannot drift apart: a kind reclassified in one place changes both.
    """
    return tuple(
        kind.name
        for kind in KINDS
        if kind.file_class in DROPPED_CLASSES and kind.name not in PIPELINE_KEEPS_ANYWAY
    )
