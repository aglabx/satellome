"""Putting back what compaction removed.

Each recipe is a function of files that stay, plus — for the per-copy layer —
the decomposer that produced them in the first place. Nothing here consults the
genome; that is the property that lets compaction run over a corpus whose source
assemblies live on a different volume.

What was verified on real output before any of this was written (assembly
``GCA_022385595.1``, 13,837 arrays):

* ``classify_from_master`` reproduces all **8** derived ``.sat`` views and all
  **9** ``.gff`` files byte-identically from the master alone.
* ``lengths_from_decomposed`` reproduces ``fastan/*.lengths`` byte-identically.
* the decomposer is **deterministic** — three runs at ``-t 4`` and two at
  ``-t 1`` give byte-identical output — and **per-array independent**: run over
  only the 13,233 sub-kilobase arrays it reproduces their monomer rows
  byte-identically, 13,233 of 13,233, against the full-corpus run.
* the corpus files are in **input order**, i.e. the master's order, with each
  array's rows contiguous. That is why re-inserting regenerated rows needs no
  stored positions: the master already says where they go.

The one thing no recipe can do offline is rebuild ``fastan/*.1aln``; it needs the
genome. That is recorded as unrestorable rather than quietly omitted.
"""

import logging
import os
import shutil
import subprocess
import tempfile

import yaml

from satellome.compact import formats
from satellome.core_functions.tools.env_check import find_console_script

logger = logging.getLogger("satellome")

#: Files ``classify_trf_data`` produces, keyed by the policy kind name.
CLASSIFY_OUTPUTS = {
    "sat_1kb": "{prefix}.1kb.sat",
    "sat_10kb": "{prefix}.10kb.sat",
    "sat_100kb": "{prefix}.100kb.sat",
    "sat_1000kb": "{prefix}.1000kb.sat",
    "sat_micro": "{prefix}.micro.sat",
    "sat_pmicro": "{prefix}.pmicro.sat",
    "sat_tssr": "{prefix}.tssr.sat",
    "sat_complex": "{prefix}.complex.sat",
}


class RecipeError(Exception):
    """A recipe cannot run, with a reason the operator can act on."""


class Unrestorable(RecipeError):
    """This kind cannot be rebuilt from what compaction keeps. Stated, not hidden."""


def genome_size_from_results(run_dir):
    """Read the assembly length out of ``results.yaml``.

    ``classify_trf_data`` needs it to compute the ``pgenome`` figures in
    ``results.yaml``; it does not affect the ``.sat``/``.gff`` bytes, but a wrong
    value would silently change the statistics a re-run writes.
    """
    path = os.path.join(run_dir, "results.yaml")
    if not os.path.exists(path):
        raise RecipeError(f"no results.yaml in {run_dir}: cannot recover the genome size")
    try:
        with open(path, "r") as fh:
            data = yaml.safe_load(fh)
    except yaml.YAMLError as e:
        raise RecipeError(f"results.yaml in {run_dir} is not readable YAML: {e}") from e
    try:
        work = data["work_files"]
        dataset = work.get("ref_assembly_name_for_trf", "dataset")
        return int(work["assembly_stats"][dataset]["genome_size"])
    except (KeyError, TypeError, ValueError) as e:
        raise RecipeError(
            f"results.yaml in {run_dir} does not carry work_files."
            f"assembly_stats.<dataset>.genome_size ({e})"
        ) from e


class ClassifyWorkspace:
    """One run of the classifier over a restored master, reused by every view.

    Regenerating eight ``.sat`` views and nine ``.gff`` files means running the
    classifier once, not seventeen times — and running the *same* classifier the
    original run used, so fidelity is a property of code reuse rather than of a
    reimplementation that has to be kept in step.
    """

    def __init__(self, write_master, prefix, genome_size, workdir):
        self._write_master = write_master
        self.prefix = prefix
        self.genome_size = genome_size
        self.workdir = workdir
        self._ready = False

    def build(self):
        if self._ready:
            return self
        from satellome.steps.trf_classify import classify_trf_data

        os.makedirs(self.workdir, exist_ok=True)
        master = os.path.join(self.workdir, f"{self.prefix}.sat")
        with open(master, "wb") as fh:
            self._write_master(fh)
        classify_trf_data(
            os.path.join(self.workdir, self.prefix), self.workdir, self.genome_size
        )
        self._ready = True
        return self

    def path_for(self, kind_name, basename):
        """Absolute path of the regenerated file for *kind_name*."""
        self.build()
        if kind_name in CLASSIFY_OUTPUTS:
            return os.path.join(
                self.workdir, CLASSIFY_OUTPUTS[kind_name].format(prefix=self.prefix)
            )
        if kind_name == "gff":
            return os.path.join(self.workdir, "gff3", basename)
        raise RecipeError(f"classify_from_master does not produce {kind_name!r}")


def lengths_from_decomposed(decomposed_records, out_fh):
    """Write ``.lengths`` from decomposed FASTA records."""
    return formats.lengths_from_decomposed_records(decomposed_records, out_fh)


# --------------------------------------------------------------------------
# per-copy layer
# --------------------------------------------------------------------------

#: Points the recipes at a specific decomposer build. The corpus was produced
#: by one version of arraysplitter and a different version does not reproduce
#: it, so an operator restoring an old catalogue needs to be able to say which
#: binary to use rather than relying on whatever is first on PATH.
ARRAYSPLITTER_ENV = "SATELLOME_ARRAYSPLITTER"


def arraysplitter_binary():
    """Locate ``arraysplitter`` or say what is missing.

    ``$SATELLOME_ARRAYSPLITTER`` wins when set, so the exact build that produced
    a catalogue can be named explicitly.

    Resolved through :func:`find_console_script`, never bare ``shutil.which``: a
    ``pip install --user`` puts the launcher in ``~/.local/bin``, which is
    routinely absent from ``PATH``, and reporting an installed tool as missing
    here would turn a reversible compaction into a refusal for no reason.
    """
    override = os.environ.get(ARRAYSPLITTER_ENV)
    if override:
        if not os.path.exists(override):
            raise Unrestorable(
                f"{ARRAYSPLITTER_ENV} points at {override}, which does not exist"
            )
        return override

    location = find_console_script("arraysplitter")
    if not location.found:
        raise Unrestorable(
            "arraysplitter is not installed, so the per-copy rows of "
            "sub-threshold arrays cannot be recomputed. Install it with "
            "'pip install arraysplitter'."
        )
    if location.hidden:
        logger.warning(
            f"arraysplitter: found at {os.path.dirname(location.path)}, but that "
            f"directory is not on your PATH (satellome can still use it)"
        )
    return location.path


def write_arrays_fasta(master_rows, out_fh, keep=None):
    """Write the arrays FASTA the decomposer consumes, from master rows.

    The header is ``>{chrom}_{start}_{end}_{array_length}_{period}`` with a
    0-based start, matching what ``bed-extract`` writes during a run; the body is
    the master's ``trf_array`` column. This is the whole basis of the policy: the
    master carries every array's complete sequence, so the per-copy layer is
    recomputable from the master alone.

    Args:
        master_rows: iterable of raw ``.sat`` data lines (bytes, no terminator).
        keep: ``callable(array_length) -> bool``, or None for every array.

    Returns:
        list of array ids written, in master order.
    """
    ids = []
    for line in master_rows:
        fields = line.split(b"\t")
        if len(fields) < 15:
            continue
        try:
            l_ind = int(fields[3])
            r_ind = int(fields[4])
            array_length = int(fields[14])
        except ValueError:
            continue
        if keep is not None and not keep(array_length):
            continue
        head, period, array = fields[2], fields[5], fields[11]
        array_id = b"%s_%d_%d_%d_%s" % (head, l_ind - 1, r_ind, array_length, period)
        out_fh.write(b">%s\n%s\n" % (array_id, array))
        ids.append(array_id)
    return ids


def run_decomposer(arrays_fasta, out_prefix, threads=1, binary=None):
    """Run the decomposer over *arrays_fasta*, writing ``<out_prefix>.*``.

    Single-threaded by default: the content is thread-count independent (checked
    across ``-t 1`` and ``-t 4``), but the order rows come out in is not, and
    ``expand`` has to reproduce a file, not merely its contents.
    """
    binary = binary or arraysplitter_binary()
    command = [binary, "-i", arrays_fasta, "-o", out_prefix, "-t", str(threads)]
    completed = subprocess.run(command, capture_output=True, text=True)
    if completed.returncode != 0:
        raise RecipeError(
            f"arraysplitter failed with exit code {completed.returncode}: "
            f"{(completed.stderr or completed.stdout or '').strip()[:400]}"
        )
    return out_prefix


def group_rows_by_array(path):
    """Read a per-copy TSV into ``(header_line, {array_id: [raw lines]})``."""
    groups = {}
    with open(path, "rb") as fh:
        header = fh.readline()
        for raw in fh:
            key = raw.split(b"\t", 1)[0]
            groups.setdefault(key, []).append(raw)
    return header, groups


def group_fasta_by_array(path):
    """Read a decomposed FASTA into ``{array_id: [(header, body)]}``."""
    groups = {}
    with open(path, "rb") as fh:
        for head, body in formats.iter_fasta_records(fh):
            key = head[1:].split(None, 1)[0] if head.startswith(b">") else b""
            groups.setdefault(key, []).append((head, body))
    return groups


class DecomposerWorkspace:
    """One decomposer run over the sub-threshold arrays, reused by every kind.

    monomers, hors and decomposed all come out of the same invocation, so an
    expand that needs all three pays for one run, not three.
    """

    def __init__(self, master_rows_factory, prefix, workdir, min_array_length,
                 threads=1):
        self._rows = master_rows_factory
        self.prefix = prefix
        self.workdir = workdir
        self.min_array_length = min_array_length
        self.threads = threads
        self._built = False
        self.array_ids = []

    def build(self):
        if self._built:
            return self
        os.makedirs(self.workdir, exist_ok=True)
        arrays_fasta = os.path.join(self.workdir, "sub_threshold.fasta")
        with open(arrays_fasta, "wb") as fh:
            self.array_ids = write_arrays_fasta(
                self._rows(), fh, keep=lambda n: n < self.min_array_length
            )
        if self.array_ids:
            run_decomposer(
                arrays_fasta,
                os.path.join(self.workdir, self.prefix),
                threads=self.threads,
            )
        self._built = True
        return self

    def produced(self, suffix):
        """Path of a decomposer output, or None when no array was below threshold."""
        self.build()
        if not self.array_ids:
            return None
        path = os.path.join(self.workdir, f"{self.prefix}{suffix}")
        if not os.path.exists(path):
            raise RecipeError(
                f"the decomposer did not produce {os.path.basename(path)}; this "
                f"build of arraysplitter wrote {sorted(os.listdir(self.workdir))}. "
                f"A different arraysplitter version cannot reproduce this corpus."
            )
        return path


def merge_per_copy_rows(kept_lines, regenerated_groups, array_order, out_fh,
                        header_line):
    """Interleave kept and regenerated rows back into the master's array order.

    Args:
        kept_lines: iterable of raw kept rows (bytes, no terminator), in order.
        regenerated_groups: ``{array_id: [raw lines, terminator included]}``.
        array_order: array ids in the master's order — where each block goes.
        out_fh: binary handle.
        header_line: the original header row, terminator included.

    Returns:
        A digest dict of what was written, so the caller can compare it with the
        md5 recorded before compaction instead of trusting the merge.
    """
    writer = formats.DigestWriter(out_fh)
    writer.write(header_line)

    kept_by_array = {}
    for line in kept_lines:
        kept_by_array.setdefault(line.split(b"\t", 1)[0], []).append(line + b"\n")

    for array_id in array_order:
        block = kept_by_array.get(array_id) or regenerated_groups.get(array_id)
        if not block:
            continue
        for raw in block:
            writer.write(raw)
    return writer.result()


def merge_fasta_records(kept_records, regenerated_groups, array_order, out_fh):
    """Same as :func:`merge_per_copy_rows` for the decomposed FASTA."""
    writer = formats.DigestWriter(out_fh)
    kept_by_array = {}
    for head, body in kept_records:
        key = head[1:].split(None, 1)[0] if head.startswith(b">") else b""
        kept_by_array.setdefault(key, []).append((head, body))

    for array_id in array_order:
        block = kept_by_array.get(array_id) or regenerated_groups.get(array_id)
        if not block:
            continue
        for head, body in block:
            writer.write(head)
            writer.write(body)
    return writer.result()


def master_array_order(master_rows):
    """Array ids in the master's own order — where every regenerated block goes."""
    order = []
    for line in master_rows:
        fields = line.split(b"\t")
        if len(fields) < 15:
            continue
        try:
            l_ind = int(fields[3])
            r_ind = int(fields[4])
            array_length = int(fields[14])
        except ValueError:
            continue
        order.append(
            b"%s_%d_%d_%d_%s" % (fields[2], l_ind - 1, r_ind, array_length, fields[5])
        )
    return order


def temp_workspace(run_dir, name):
    """A scratch directory beside the run, never in the system temp.

    Beside the run so that a 3 GB restore does not land on a small ``/tmp``, and
    so an interrupted expand leaves its debris where the operator is looking.
    """
    parent = os.path.join(run_dir, ".compact_work")
    os.makedirs(parent, exist_ok=True)
    return tempfile.mkdtemp(prefix=f"{name}_", dir=parent)


def drop_workspace(run_dir):
    parent = os.path.join(run_dir, ".compact_work")
    if os.path.isdir(parent):
        shutil.rmtree(parent, ignore_errors=True)
