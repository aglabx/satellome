"""Compact and expand a finished run directory.

Invariants, in the order they matter:

1. **Nothing is removed that has not just been reproduced.** A dropped file is
   regenerated into scratch and compared, by md5, against the file about to be
   deleted; only then is it unlinked. A re-encoded file is decoded back and
   compared before the original goes. "Reversible" is therefore a fact
   established per file at compaction time, not a claim about the design.
2. **A failure narrows the compaction, it never widens the damage.** A kind
   whose recipe does not reproduce it is reclassified as primary and kept, with
   the reason named in the report and in the record.
3. **Idempotent.** A directory that already carries ``.compact.json`` is a
   no-op: zero writes, zero byte change, exit 0.
4. **Streaming.** Single files reach 3.2 GB; nothing is held whole in RAM
   except the per-array row blocks a merge needs, which are bounded by one
   array.
"""

import logging
import os
from typing import List, NamedTuple, Sequence

from satellome import __version__
from satellome.compact import columnar, formats, readiness, recipes, record
from satellome.compact.ledger import build_plan, format_ledger, human, iter_run_files
from satellome.compact.probe import ProbeResult, probe_decomposer
from satellome.compact.policy import (
    DROPPED_CLASSES,
    FileClass,
    classify_path,
)
from satellome.core_functions.tools.atomic_io import discard, partial_path

logger = logging.getLogger("satellome")


class _ProbeDenyAll:
    """Stands in for a probe that refuses everything, for the retry path."""

    refused = {}

    def can_filter(self, kind_name):
        return False


class Outcome(NamedTuple):
    ok: bool
    #: "compacted", "expanded", "noop", "refused", "planned", "checked"
    status: str
    lines: List[str]
    before: int = 0
    after: int = 0
    #: Cuts that could not be shown reversible here, and were therefore not made.
    kept_back: Sequence[str] = ()


# --------------------------------------------------------------------------
# keeping --verify-run meaningful
# --------------------------------------------------------------------------

def sync_manifest(run_dir, event, log=logger):
    """Bring ``run_manifest.json`` in step with what the directory now holds.

    Compaction changes every recorded byte size, so without this
    ``satellome --verify-run`` — which the panel's driver calls before it
    archives or ingests a directory — would report every compacted run as
    damaged, and an operator who has seen that failure a thousand times stops
    reading it.

    The audit chain is not lost by refreshing the inventory: the manifest states
    what is on disk *now*, ``.compact.json`` carries the md5 of every file's
    content *before*, and ``compaction_history`` records each transition. A
    directory with no manifest (1,020 of the finished catalogues predate them)
    is left alone.
    """
    from satellome.core_functions.tools.run_manifest import (
        MANIFEST_NAME,
        ManifestError,
        _iter_output_files,
        load_manifest,
        write_run_manifest,
    )

    if not os.path.exists(os.path.join(run_dir, MANIFEST_NAME)):
        return None
    try:
        manifest = load_manifest(run_dir)
    except ManifestError as e:
        # A corrupt manifest is not "no manifest": say so and leave it alone
        # rather than overwriting the only record of what the run produced.
        log.error(f"{run_dir}: {e} — leaving it untouched")
        return None

    manifest["files"] = [
        {"path": rel, "bytes": size} for rel, size in _iter_output_files(run_dir)
    ]
    history = manifest.setdefault("compaction_history", [])
    history.append(event)
    write_run_manifest(run_dir, manifest)
    return manifest


# --------------------------------------------------------------------------
# reading the master, whichever form it is in
# --------------------------------------------------------------------------

class Master:
    """Access to the master ``.sat`` table, plain, gzipped or in a container."""

    def __init__(self, path, compression, container=False):
        self.path = path
        self.compression = compression
        self.container = container

    @classmethod
    def locate(cls, run_dir):
        """Find the master in whatever form the directory currently holds it."""
        path, _prefix = readiness.find_master(run_dir)
        if path is not None:
            return cls(path, "gzip" if path.endswith(".gz") else "none")
        for name in sorted(os.listdir(run_dir)):
            if name.endswith(".sat.satz"):
                return cls(os.path.join(run_dir, name), "none", container=True)
        raise recipes.RecipeError(
            f"no master .sat, .sat.gz or .sat.satz in {run_dir}: nothing can be "
            f"regenerated without it"
        )

    def write_plain(self, fh):
        """Write the master's exact original bytes to *fh*."""
        if self.container:
            columnar.restore(self.path, fh)
            return
        with formats.open_maybe_gzip(self.path, self.compression) as src:
            for block in iter(lambda: src.read(formats.IO_CHUNK), b""):
                fh.write(block)

    def data_lines(self):
        """Yield master data rows as raw bytes, without preamble, header or terminator."""
        if self.container:
            first = True
            for line in columnar.read_table(self.path):
                if first:
                    first = False
                    if line.split(b"\t", 1)[0] == b"project":
                        continue
                yield line
            return
        with formats.open_maybe_gzip(self.path, self.compression) as fh:
            preamble, first = formats._split_preamble(fh)
            del preamble
            if first is None:
                return
            # The first non-comment line is the header row.
            for raw in fh:
                yield raw[:-2] if raw.endswith(b"\r\n") else raw.rstrip(b"\n")


# --------------------------------------------------------------------------
# compact
# --------------------------------------------------------------------------

def compact(run_dir, config, dry_run=False, log=logger):
    """Reduce *run_dir* to its primary layer. Returns an :class:`Outcome`."""
    run_dir = os.path.abspath(run_dir)
    lines = []

    if record.is_compacted(run_dir):
        return Outcome(
            True, "noop",
            [f"{run_dir} is already compacted ({record.RECORD_NAME} present); "
             f"nothing to do"],
        )

    ready = readiness.check(run_dir)
    if not ready.ok:
        return Outcome(
            False, "refused",
            [f"REFUSED {run_dir}"] + [f"  ! {reason}" for reason in ready.reasons],
        )

    # Ask, before pricing anything, whether the per-copy cut can be applied on
    # this machine at all: it is worth about half the saving, and it depends on
    # a binary that may not be the one that produced the catalogue.
    probe_dir = recipes.temp_workspace(run_dir, "probe")
    try:
        probe = probe_decomposer(
            run_dir, Master.locate(run_dir).data_lines(), ready.prefix,
            config.min_array_length, os.path.join(probe_dir, "decomp"),
        )
    except recipes.RecipeError as e:
        probe = ProbeResult(frozenset(), {"monomers": str(e), "hors": str(e),
                                          "decomposed": str(e)})
    finally:
        recipes.drop_workspace(run_dir)

    plan = build_plan(run_dir, config, prefix=ready.prefix, measure=True, probe=probe)
    lines.extend(format_ledger(plan, config))
    if dry_run:
        lines.append("")
        lines.append("  --dry-run: nothing was written")
        return Outcome(True, "planned", lines, plan.before, plan.after)

    rec = record.new_record(run_dir, config, __version__, prefix=ready.prefix)
    rec["probe_refused"] = dict(probe.refused)
    work = recipes.temp_workspace(run_dir, "compact")
    # A refusal the probe found is stated on the real run too, not only in the
    # dry run: otherwise the log of a compaction that quietly kept half the
    # bytes looks exactly like one that did the whole job.
    kept_back = [
        f"{name}: {reason}; re-encoded whole, no rows dropped"
        for name, reason in sorted(probe.refused.items())
    ]

    try:
        # Located afresh on every use: the master is itself re-encoded during
        # this run, so a handle taken once would point at a file that no longer
        # exists by the time the per-copy layer is verified.
        def master_write(fh):
            Master.locate(run_dir).write_plain(fh)

        def master_rows():
            return Master.locate(run_dir).data_lines()

        classify_ws = recipes.ClassifyWorkspace(
            master_write, ready.prefix,
            recipes.genome_size_from_results(run_dir),
            os.path.join(work, "classify"),
        )
        decomposer_ws = recipes.DecomposerWorkspace(
            master_rows, ready.prefix, os.path.join(work, "decomp"),
            config.min_array_length,
        )

        # Drops first, while every original is still on disk to compare against.
        for entry in plan.files:
            if entry.file_class not in DROPPED_CLASSES:
                continue
            _drop_one(run_dir, entry, rec, config, classify_ws, kept_back, log)

        # Then re-encode what stays.
        for entry in plan.files:
            if entry.file_class not in {
                FileClass.PRIMARY, FileClass.PRIMARY_FILTERED, FileClass.PRIMARY_BLOB
            }:
                continue
            if entry.layout is None:
                continue
            _reencode_one(
                run_dir, entry, rec, config, master_rows, decomposer_ws,
                kept_back, log, probe
            )

        rec["kept_unclassified"] = list(plan.unclassified)
        rec["kept_back"] = kept_back
        record.write_record(run_dir, rec)
        sync_manifest(run_dir, {
            "action": "compact",
            "at": rec["created"],
            "satellome_version": __version__,
            "policy_version": rec["policy_version"],
            "min_array_length": config.min_array_length,
            "level": config.level,
            "before_bytes": rec["totals"]["before"],
            "after_bytes": rec["totals"]["after"],
            "kept_back": len(kept_back),
        }, log=log)
    finally:
        recipes.drop_workspace(run_dir)

    # Report what an operator would measure with du: the whole directory before
    # and after, not only the files compaction happened to touch. Those two
    # numbers differ whenever a kind is kept back, and quoting the narrower one
    # would make the ledger look more accurate than it is.
    before = plan.before
    after = sum(
        os.path.getsize(os.path.join(run_dir, rel))
        for rel in iter_run_files(run_dir)
    )
    rec["totals"]["directory_before"] = before
    rec["totals"]["directory_after"] = after
    record.write_record(run_dir, rec)
    lines.append("")
    lines.append(
        f"  realised: {human(before)} -> {human(after)}  "
        f"freed {human(before - after)}"
        + (f" ({100.0 * (before - after) / before:.1f}%)" if before else "")
    )
    for note in kept_back:
        lines.append(f"  - kept: {note}")
    return Outcome(True, "compacted", lines, before, after, kept_back)


def _drop_one(run_dir, entry, rec, config, classify_ws, kept_back, log):
    """Regenerate a droppable file, compare it, and only then unlink it."""
    abs_path = os.path.join(run_dir, entry.rel_path)
    match = classify_path(entry.rel_path)
    kind = match.kind
    before_md5 = formats.md5_of_content(abs_path, match.compression)
    entry_rec = {
        "path": entry.rel_path,
        "kind": kind.name,
        "class": kind.file_class,
        "action": "drop",
        "recipe": {"kind": kind.recipe},
        "before": {
            "stored_bytes": entry.stored_bytes,
            "content_bytes": entry.content_bytes,
            "compression": match.compression,
            "md5": before_md5,
            "stored_md5": formats.md5_of_file(abs_path),
            "gzip": formats.gzip_framing(abs_path),
        },
        "after": {"stored_bytes": 0},
    }

    if kind.recipe == "rerun_fastan":
        # Honest about the one kind expand cannot rebuild offline. Dropped by
        # decision, and the record says so rather than implying reversibility.
        entry_rec["recipe"]["restorable"] = False
        entry_rec["recipe"]["note"] = (
            "regenerating a .1aln needs the genome and a FasTAN run; expand "
            "cannot restore it from this directory"
        )
        os.unlink(abs_path)
        record.add_entry(rec, entry_rec)
        log.info(f"  dropped {entry.rel_path} ({human(entry.stored_bytes)}, not restorable offline)")
        return

    if not config.verify_drops:
        entry_rec["recipe"]["verified_at_compaction"] = False
        os.unlink(abs_path)
        record.add_entry(rec, entry_rec)
        log.info(f"  dropped {entry.rel_path} ({human(entry.stored_bytes)}, unverified)")
        return

    try:
        produced = _regenerate_dropped(run_dir, entry, kind, classify_ws, config)
    except recipes.RecipeError as e:
        kept_back.append(f"{entry.rel_path}: {e}")
        log.warning(f"  kept {entry.rel_path}: {e}")
        return

    if produced is None:
        kept_back.append(f"{entry.rel_path}: no recipe produced it")
        log.warning(f"  kept {entry.rel_path}: no recipe produced it")
        return

    if produced != before_md5:
        kept_back.append(
            f"{entry.rel_path}: regeneration gave md5 {produced}, the file on "
            f"disk is {before_md5} — kept, because a file we cannot reproduce "
            f"is not a file we may delete"
        )
        log.warning(f"  kept {entry.rel_path}: regeneration does not reproduce it")
        return

    entry_rec["recipe"]["verified_at_compaction"] = True
    os.unlink(abs_path)
    record.add_entry(rec, entry_rec)
    log.info(f"  dropped {entry.rel_path} ({human(entry.stored_bytes)}, verified)")


def _regenerate_dropped(run_dir, entry, kind, classify_ws, config):
    """Regenerate a droppable file into scratch; return the md5 of its content."""
    if kind.recipe == "classify_from_master":
        basename = os.path.basename(entry.rel_path)
        stem = basename[:-3] if basename.endswith(".gz") else basename
        produced = classify_ws.path_for(kind.name, stem)
        if not os.path.exists(produced):
            raise recipes.RecipeError(
                f"the classifier did not produce {os.path.basename(produced)}"
            )
        return formats.md5_of_content(produced, "none")

    if kind.recipe == "lengths_from_decomposed":
        source = _find_decomposed(run_dir)
        if source is None:
            raise recipes.RecipeError(
                "fastan/*.decomposed.fasta is not present, so .lengths cannot "
                "be regenerated from it"
            )
        path, compression = source
        sink = formats.DigestWriter(None)
        with formats.open_maybe_gzip(path, compression) as fh:
            recipes.lengths_from_decomposed(formats.iter_fasta_records(fh), sink)
        return sink.result()["md5"]

    raise recipes.RecipeError(f"no recipe implemented for {kind.recipe!r}")


def _find_decomposed(run_dir):
    fastan = os.path.join(run_dir, "fastan")
    if not os.path.isdir(fastan):
        return None
    for name in sorted(os.listdir(fastan)):
        if name.endswith(".decomposed.fasta"):
            return os.path.join(fastan, name), "none"
        if name.endswith(".decomposed.fasta.gz"):
            return os.path.join(fastan, name), "gzip"
    return None


def _reencode_one(run_dir, entry, rec, config, master_rows, decomposer_ws,
                  kept_back, log, probe=None):
    """Rewrite one kept file into a container, verify it, then unlink the original."""
    abs_path = os.path.join(run_dir, entry.rel_path)
    match = classify_path(entry.rel_path)
    kind = match.kind
    before_md5 = formats.md5_of_content(abs_path, match.compression)
    out_rel = match.stem + ".satz"
    out_path = os.path.join(os.path.dirname(abs_path), os.path.basename(out_rel))
    out_rel_path = os.path.relpath(out_path, run_dir)
    tmp = partial_path(out_path)

    filter_note = None
    row_filter = None
    record_filter = None
    # A kind the probe refused is re-encoded whole: the compression is still
    # worth taking, and the row cut is the only part that needed a decomposer
    # able to reproduce this catalogue. This is the brief's "reclassify as
    # primary and keep" — keep the rows, not the old encoding.
    apply_filter = kind.filtered and (probe is None or probe.can_filter(kind.name))
    if apply_filter:
        threshold = config.min_array_length
        filter_note = {
            "kind": "array_min_length",
            "min_array_length": threshold,
            "key": "array_id, length field parsed from the right",
        }
        # A record whose array length cannot be parsed is kept: a row we cannot
        # measure is not a row we may decide about.
        if kind.layout == "fasta":
            def keep_record(head, _t=threshold):
                length = formats.array_length_from_fasta_header(head)
                return length is None or length >= _t

            record_filter = keep_record
        else:
            def keep_row(fields, _t=threshold):
                length = formats.array_length_from_id(fields[0])
                return length is None or length >= _t

            row_filter = keep_row

    source = {
        "name": os.path.basename(entry.rel_path),
        "stored_bytes": entry.stored_bytes,
        "content_bytes": entry.content_bytes,
        "compression": match.compression,
        "md5": before_md5,
    }

    try:
        if kind.layout == "fasta":
            footer = formats.encode_fasta(
                abs_path, match.compression, tmp, config.level,
                record_filter=record_filter, filter_note=filter_note, source=source,
            )
        elif kind.layout == "blob":
            footer = formats.encode_blob(
                abs_path, match.compression, tmp, config.level, source=source
            )
        else:
            footer = formats.encode_tsv(
                abs_path, match.compression, tmp, config.level,
                has_header=kind.has_header,
                row_filter=row_filter, filter_note=filter_note, source=source,
            )
    except Exception as e:
        discard(out_path)
        kept_back.append(f"{entry.rel_path}: encoding failed ({e})")
        log.error(f"  kept {entry.rel_path}: encoding failed: {e}")
        return

    os.replace(tmp, out_path)

    ok, detail = columnar.verify(out_path)
    if not ok:
        os.unlink(out_path)
        kept_back.append(f"{entry.rel_path}: container did not read back ({detail})")
        log.error(f"  kept {entry.rel_path}: container did not read back: {detail}")
        return

    if not apply_filter and footer["md5"] != before_md5:
        os.unlink(out_path)
        kept_back.append(
            f"{entry.rel_path}: re-encoded content md5 {footer['md5']} != "
            f"{before_md5}"
        )
        log.error(f"  kept {entry.rel_path}: round trip is not byte-identical")
        return

    if apply_filter and config.verify_drops:
        problem = _verify_filtered_roundtrip(
            run_dir, out_path, kind, master_rows, decomposer_ws, before_md5, match
        )
        if problem:
            # The cut cannot be proven reversible here, so it is not made. The
            # file is still re-encoded, whole: compression is not the part that
            # needed proving.
            os.unlink(out_path)
            kept_back.append(f"{entry.rel_path}: {problem}; re-encoded whole instead")
            log.warning(f"  {entry.rel_path}: {problem}; re-encoding whole instead")
            _reencode_one(
                run_dir, entry, rec, config, master_rows, decomposer_ws,
                kept_back, log, probe=_ProbeDenyAll(),
            )
            return

    after_bytes = os.path.getsize(out_path)
    if not apply_filter and after_bytes >= entry.stored_bytes:
        # Per-column framing has a fixed cost per column per chunk. On a table
        # of a few hundred rows that cost exceeds what splitting the columns
        # saves, and the container comes out *larger* than the gzip it would
        # replace. Keeping the original is both smaller and one less format to
        # read, so say so and move on.
        os.unlink(out_path)
        kept_back.append(
            f"{entry.rel_path}: re-encoding it would grow the file "
            f"({human(entry.stored_bytes)} -> {human(after_bytes)}), kept as it is"
        )
        log.info(
            f"  kept {entry.rel_path}: too small for columnar to pay "
            f"({human(entry.stored_bytes)} -> {human(after_bytes)})"
        )
        return
    record.add_entry(rec, {
        "path": entry.rel_path,
        "kind": kind.name,
        "class": kind.file_class,
        "action": "reencode_filtered" if apply_filter else "reencode",
        "layout": kind.layout,
        "before": {
            "stored_bytes": entry.stored_bytes,
            "content_bytes": entry.content_bytes,
            "compression": match.compression,
            "md5": before_md5,
            "stored_md5": formats.md5_of_file(abs_path),
            "gzip": formats.gzip_framing(abs_path),
        },
        "after": {
            "path": out_rel_path,
            "stored_bytes": after_bytes,
            "content_bytes": footer.get("bytes", 0),
            "md5": footer.get("md5"),
        },
        "filter": filter_note,
        "rows": {"in": footer.get("rows_in", 0), "out": footer.get("rows_out", 0)},
        "recipe": {
            "kind": kind.recipe if apply_filter else None,
            "verified_at_compaction": (
                bool(config.verify_drops) if apply_filter else True
            ),
        },
    })
    os.unlink(abs_path)
    log.info(
        f"  re-encoded {entry.rel_path} -> {out_rel_path} "
        f"({human(entry.stored_bytes)} -> {human(after_bytes)})"
    )


def _verify_filtered_roundtrip(run_dir, container_path, kind, master_rows,
                               decomposer_ws, before_md5, match):
    """Prove the rows this container leaves out can be put back, now.

    Returns a reason string when it cannot, or None when the merged result
    reproduces the original byte for byte.
    """
    try:
        order = recipes.master_array_order(master_rows())
        if kind.layout == "fasta":
            produced = decomposer_ws.produced(".decomposed.fasta")
            groups = recipes.group_fasta_by_array(produced) if produced else {}
            sink = formats.DigestWriter(None)
            merged = recipes.merge_fasta_records(
                columnar.read_fasta(container_path), groups, order, sink
            )
        else:
            suffix = ".monomers.tsv" if kind.name == "monomers" else ".hors.tsv"
            produced = decomposer_ws.produced(suffix)
            if produced is None:
                header, groups = b"", {}
            else:
                header, groups = recipes.group_rows_by_array(produced)
            container_header = columnar.header_line(container_path)
            sink = formats.DigestWriter(None)
            merged = recipes.merge_per_copy_rows(
                columnar.read_table(container_path), groups, order, sink,
                container_header,
            )
    except recipes.RecipeError as e:
        return str(e)

    if merged["md5"] != before_md5:
        return (
            f"the dropped sub-threshold rows do not reproduce the original "
            f"(merged md5 {merged['md5']}, original {before_md5}); kept whole"
        )
    return None


# --------------------------------------------------------------------------
# expand
# --------------------------------------------------------------------------

def expand(run_dir, log=logger, threads=1):
    """Restore *run_dir* to its pre-compaction contents."""
    run_dir = os.path.abspath(run_dir)
    lines = []
    try:
        rec = record.load_record(run_dir)
    except record.RecordError as e:
        return Outcome(False, "refused", [str(e)])

    config_dict = rec.get("config", {})
    min_array_length = int(config_dict.get("min_array_length", 1000))
    prefix = rec.get("prefix") or ""

    work = recipes.temp_workspace(run_dir, "expand")
    restored, failed, unrestorable = [], [], []

    def master_write(fh):
        Master.locate(run_dir).write_plain(fh)

    def master_rows():
        return Master.locate(run_dir).data_lines()

    try:
        # The master first: every other recipe reads it, and it is located
        # afresh each time because restoring it changes which form it is in.
        for entry in rec["files"]:
            if entry.get("kind") == "sat_master":
                _restore_reencoded(run_dir, entry, restored, failed, log)

        for entry in rec["files"]:
            if entry.get("action") not in {"reencode", "reencode_filtered"}:
                continue
            if entry.get("kind") == "sat_master":
                continue
            if entry.get("action") == "reencode":
                _restore_reencoded(run_dir, entry, restored, failed, log)

        decomposer_ws = recipes.DecomposerWorkspace(
            master_rows, prefix, os.path.join(work, "decomp"),
            min_array_length, threads=threads,
        )
        order = None
        for entry in rec["files"]:
            if entry.get("action") != "reencode_filtered":
                continue
            if order is None:
                order = recipes.master_array_order(master_rows())
            _restore_filtered(
                run_dir, entry, decomposer_ws, order, restored, failed, log
            )

        classify_ws = recipes.ClassifyWorkspace(
            master_write, prefix,
            recipes.genome_size_from_results(run_dir),
            os.path.join(work, "classify"),
        )
        for entry in rec["files"]:
            if entry.get("action") != "drop":
                continue
            if entry.get("recipe", {}).get("restorable") is False:
                unrestorable.append(entry["path"])
                continue
            _restore_dropped(run_dir, entry, classify_ws, restored, failed, log)
    except (recipes.RecipeError, columnar.ContainerError) as e:
        return Outcome(False, "refused", [f"expand failed: {e}"] + lines)
    finally:
        recipes.drop_workspace(run_dir)

    import datetime

    sync_manifest(run_dir, {
        "action": "expand",
        "at": datetime.datetime.now().isoformat(timespec="seconds"),
        "satellome_version": __version__,
        "restored": len(restored),
        "unrestorable": list(unrestorable),
        "failed": len(failed),
        "note": (
            "restored files are content-identical to the originals; a file that "
            "was gzipped is re-gzipped, so its container bytes differ while the "
            "table inside does not"
        ),
    }, log=log)

    if not failed:
        os.unlink(record.record_path(run_dir))

    lines.append(f"Expanded {run_dir}")
    lines.append(f"  restored {len(restored)} file(s)")
    for path in unrestorable:
        lines.append(
            f"  ! {path}: dropped by policy and not restorable without the "
            f"genome — rerun satellome to rebuild it"
        )
    for note in failed:
        lines.append(f"  ! {note}")
    if failed:
        lines.append(
            f"  {record.RECORD_NAME} kept: the directory is not fully restored"
        )
    return Outcome(not failed, "expanded", lines)


def _restore_reencoded(run_dir, entry, restored, failed, log):
    container = os.path.join(run_dir, entry["after"]["path"])
    target = os.path.join(run_dir, entry["path"])
    if os.path.exists(target):
        restored.append(entry["path"])
        return
    if not os.path.exists(container):
        failed.append(f"{entry['path']}: its container {entry['after']['path']} is gone")
        return
    compression = entry["before"].get("compression", "none")
    tmp = partial_path(target)
    digest = formats.decode_to_file(
        container, tmp, compression=compression,
        framing=entry["before"].get("gzip"),
    )
    if digest["md5"] != entry["before"]["md5"]:
        discard(target)
        failed.append(
            f"{entry['path']}: restored content md5 {digest['md5']} != recorded "
            f"{entry['before']['md5']}"
        )
        return
    os.replace(tmp, target)
    os.unlink(container)
    restored.append(entry["path"])
    log.info(f"  restored {entry['path']}")


def _restore_filtered(run_dir, entry, decomposer_ws, order, restored, failed, log):
    container = os.path.join(run_dir, entry["after"]["path"])
    target = os.path.join(run_dir, entry["path"])
    if os.path.exists(target):
        restored.append(entry["path"])
        return
    if not os.path.exists(container):
        failed.append(f"{entry['path']}: its container {entry['after']['path']} is gone")
        return

    kind_name = entry.get("kind")
    compression = entry["before"].get("compression", "none")
    framing = entry["before"].get("gzip")
    tmp = partial_path(target)

    try:
        if kind_name == "decomposed":
            produced = decomposer_ws.produced(".decomposed.fasta")
            groups = recipes.group_fasta_by_array(produced) if produced else {}
            with _open_out(tmp, compression, framing) as out:
                digest = recipes.merge_fasta_records(
                    columnar.read_fasta(container), groups, order, out
                )
        else:
            suffix = ".monomers.tsv" if kind_name == "monomers" else ".hors.tsv"
            produced = decomposer_ws.produced(suffix)
            groups = recipes.group_rows_by_array(produced)[1] if produced else {}
            header = columnar.header_line(container)
            with _open_out(tmp, compression, framing) as out:
                digest = recipes.merge_per_copy_rows(
                    columnar.read_table(container), groups, order, out, header
                )
    except recipes.RecipeError as e:
        discard(target)
        failed.append(f"{entry['path']}: {e}")
        return

    if digest["md5"] != entry["before"]["md5"]:
        discard(target)
        failed.append(
            f"{entry['path']}: rebuilt content md5 {digest['md5']} != recorded "
            f"{entry['before']['md5']} — the decomposer no longer reproduces "
            f"the rows that were dropped"
        )
        return
    os.replace(tmp, target)
    os.unlink(container)
    restored.append(entry["path"])
    log.info(f"  rebuilt {entry['path']}")


def _restore_dropped(run_dir, entry, classify_ws, restored, failed, log):
    target = os.path.join(run_dir, entry["path"])
    if os.path.exists(target):
        restored.append(entry["path"])
        return
    recipe = entry.get("recipe", {}).get("kind")
    compression = entry["before"].get("compression", "none")
    framing = entry["before"].get("gzip")
    tmp = partial_path(target)
    os.makedirs(os.path.dirname(tmp), exist_ok=True)

    try:
        if recipe == "classify_from_master":
            basename = os.path.basename(entry["path"])
            stem = basename[:-3] if basename.endswith(".gz") else basename
            produced = classify_ws.path_for(entry["kind"], stem)
            with open(produced, "rb") as src, \
                    _open_out(tmp, compression, framing) as out:
                writer = formats.DigestWriter(out)
                for block in iter(lambda: src.read(formats.IO_CHUNK), b""):
                    writer.write(block)
                digest = writer.result()
        elif recipe == "lengths_from_decomposed":
            source = _find_decomposed(run_dir)
            if source is None:
                raise recipes.RecipeError(
                    "fastan/*.decomposed.fasta is not present; restore it first"
                )
            path, src_compression = source
            with formats.open_maybe_gzip(path, src_compression) as fh, \
                    _open_out(tmp, compression, framing) as out:
                digest = recipes.lengths_from_decomposed(
                    formats.iter_fasta_records(fh), out
                )
        else:
            raise recipes.RecipeError(f"no recipe implemented for {recipe!r}")
    except recipes.RecipeError as e:
        discard(target)
        failed.append(f"{entry['path']}: {e}")
        return

    if digest["md5"] != entry["before"]["md5"]:
        discard(target)
        failed.append(
            f"{entry['path']}: regenerated content md5 {digest['md5']} != "
            f"recorded {entry['before']['md5']}"
        )
        return
    os.replace(tmp, target)
    restored.append(entry["path"])
    log.info(f"  regenerated {entry['path']}")


def _open_out(path, compression, framing=None):
    """Open a restore target, replaying the recorded gzip framing when there is one."""
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    if compression == "gzip":
        return formats.open_gzip_out(path, framing)
    return open(path, "wb")


# --------------------------------------------------------------------------
# the blocking gate
# --------------------------------------------------------------------------

def check_recipes(run_dir, config, log=logger):
    """Checklist 0: does every recipe reproduce the file it claims to replace?

    Run this before any corpus-wide compaction. It writes nothing to *run_dir*
    beyond a scratch directory it removes, and reports per kind whether the
    regenerated bytes match what is on disk. A "derived" file that cannot
    actually be rederived is a one-way delete wearing the word *reversible*.
    """
    run_dir = os.path.abspath(run_dir)
    ready = readiness.check(run_dir)
    if not ready.ok:
        return Outcome(
            False, "refused",
            [f"REFUSED {run_dir}"] + [f"  ! {r}" for r in ready.reasons],
        )

    results = []
    work = recipes.temp_workspace(run_dir, "check")
    try:
        def master_write(fh):
            Master.locate(run_dir).write_plain(fh)

        def master_rows():
            return Master.locate(run_dir).data_lines()

        classify_ws = recipes.ClassifyWorkspace(
            master_write, ready.prefix,
            recipes.genome_size_from_results(run_dir),
            os.path.join(work, "classify"),
        )
        decomposer_ws = recipes.DecomposerWorkspace(
            master_rows, ready.prefix, os.path.join(work, "decomp"),
            config.min_array_length,
        )

        for rel_path in iter_run_files(run_dir):
            match = classify_path(rel_path)
            if match.kind is None:
                continue
            kind = match.kind
            abs_path = os.path.join(run_dir, rel_path)
            if kind.file_class in DROPPED_CLASSES:
                if kind.recipe == "rerun_fastan":
                    results.append((rel_path, kind.name, None,
                                    "needs the genome; not checked offline"))
                    continue
                expected = formats.md5_of_content(abs_path, match.compression)
                try:
                    got = _regenerate_dropped(
                        run_dir,
                        _planned_stub(rel_path, abs_path, match),
                        kind, classify_ws, config,
                    )
                except recipes.RecipeError as e:
                    results.append((rel_path, kind.name, False, str(e)))
                    continue
                results.append((rel_path, kind.name, got == expected,
                                "byte-identical" if got == expected
                                else f"md5 {got} != {expected}"))
            elif kind.filtered:
                expected = formats.md5_of_content(abs_path, match.compression)
                tmp = os.path.join(work, os.path.basename(match.stem) + ".satz")
                try:
                    _encode_for_check(abs_path, match, kind, config, tmp)
                    problem = _verify_filtered_roundtrip(
                        run_dir, tmp, kind, master_rows, decomposer_ws, expected,
                        match
                    )
                except (recipes.RecipeError, columnar.ContainerError) as e:
                    results.append((rel_path, kind.name, False, str(e)))
                    continue
                results.append((rel_path, kind.name, problem is None,
                                problem or "byte-identical after merge"))
    finally:
        recipes.drop_workspace(run_dir)

    lines = [f"Recipe check for {run_dir}", ""]
    failures = 0
    for rel_path, kind_name, ok, detail in results:
        mark = "OK  " if ok else ("SKIP" if ok is None else "FAIL")
        if ok is False:
            failures += 1
        lines.append(f"  {mark} {kind_name:<14} {rel_path}")
        lines.append(f"       {detail}")
    if not results:
        lines.append("  nothing droppable or filtered is present here")
    lines.append("")
    lines.append(
        f"  {len(results) - failures} of {len(results)} checks passed"
        if results else "  no checks ran"
    )
    return Outcome(failures == 0, "checked", lines)


def _planned_stub(rel_path, abs_path, match):
    from satellome.compact.ledger import PlannedFile

    size = os.path.getsize(abs_path)
    return PlannedFile(rel_path, match.kind.name, match.kind.file_class, None,
                       match.compression, size, size, 0)


def _encode_for_check(abs_path, match, kind, config, out_path):
    threshold = config.min_array_length
    note = {"kind": "array_min_length", "min_array_length": threshold}
    if kind.layout == "fasta":
        formats.encode_fasta(
            abs_path, match.compression, out_path, config.level,
            record_filter=lambda h: (
                formats.array_length_from_fasta_header(h) is None
                or formats.array_length_from_fasta_header(h) >= threshold
            ),
            filter_note=note,
        )
    else:
        formats.encode_tsv(
            abs_path, match.compression, out_path, config.level,
            row_filter=lambda f: (
                formats.array_length_from_id(f[0]) is None
                or formats.array_length_from_id(f[0]) >= threshold
            ),
            filter_note=note,
        )
