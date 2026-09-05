"""What compaction would do, priced, before anything is touched.

``satellome compact --dry-run`` prints this. It is a real measurement, not a
table of assumed ratios: the freed fraction of the per-copy cut ranges from 19%
to 92% across assemblies — a T2T assembly's arrays are long and intact, so most
of its copy rows belong to arrays above the threshold — and any estimate built
on a corpus average would be wrong by 5x on individual directories.

So the estimate reads each file once: the kept fraction is counted exactly, and
the codec ratio is measured by compressing a bounded sample of the same file
rather than by guessing. What it cannot know is how the codec behaves on the
part of a large file it did not sample, which is why the ledger prints
``estimate`` and the realised figures are reported again after the run.
"""

import os
import shutil
import tempfile
from typing import List, NamedTuple, Optional

from satellome.compact import columnar, formats
from satellome.compact.policy import DROPPED_CLASSES, FileClass, classify_path

#: Content bytes compressed to measure the ratio of one file.
SAMPLE_BYTES = 32 << 20


class PlannedFile(NamedTuple):
    rel_path: str
    kind_name: str
    file_class: str
    layout: Optional[str]
    compression: str
    #: Bytes the file occupies now.
    stored_bytes: int
    #: Bytes of content (== stored_bytes unless gzipped).
    content_bytes: int
    #: Estimated bytes after compaction (0 for a drop).
    after_bytes: int
    #: Rows/records kept out of the total, for filtered kinds.
    rows_in: int = 0
    rows_out: int = 0


class Plan(NamedTuple):
    run_dir: str
    prefix: str
    files: List[PlannedFile]
    unclassified: List[str]
    #: {kind name: why the per-copy cut will not be applied here}. None rather
    #: than {} so the default is never a shared mutable object.
    refused: Optional[dict] = None

    @property
    def before(self):
        return sum(f.stored_bytes for f in self.files)

    @property
    def after(self):
        return sum(f.after_bytes for f in self.files)


def iter_run_files(run_dir):
    """Yield paths relative to *run_dir*, skipping our own scratch directory."""
    for root, dirs, files in os.walk(run_dir):
        dirs[:] = sorted(d for d in dirs if d != ".compact_work")
        for name in sorted(files):
            if name.endswith(".partial"):
                continue
            yield os.path.relpath(os.path.join(root, name), run_dir)


def _content_bytes(path, compression):
    if compression != "gzip":
        return os.path.getsize(path)
    total = 0
    with formats.open_maybe_gzip(path, compression) as fh:
        for block in iter(lambda: fh.read(formats.IO_CHUNK), b""):
            total += len(block)
    return total


def _blob_ratio(path, compression, level):
    """Compressed/uncompressed ratio of a bounded sample, compressed whole."""
    sample = bytearray()
    with formats.open_maybe_gzip(path, compression) as fh:
        while len(sample) < SAMPLE_BYTES:
            block = fh.read(min(formats.IO_CHUNK, SAMPLE_BYTES - len(sample)))
            if not block:
                break
            sample += block
    if not sample:
        return 1.0
    return len(columnar._compress(bytes(sample), level)) / len(sample)


def _columnar_ratio(path, compression, level, layout, has_header, tmp_dir):
    """Ratio measured by actually encoding a bounded prefix of the file.

    Estimating a columnar container by compressing the file whole is wrong in
    the direction that matters: it misses the entire effect being measured, and
    the ledger then under-promises by a fixed margin an operator cannot see. So
    the sample is encoded the way the real thing will be.
    """
    sample_path = os.path.join(tmp_dir, "sample" + (".fasta" if layout == "fasta" else ".tsv"))
    written = 0
    with formats.open_maybe_gzip(path, compression) as src, \
            open(sample_path, "wb") as out:
        for raw in src:
            out.write(raw)
            written += len(raw)
            if written >= SAMPLE_BYTES:
                break
    if written == 0:
        return 1.0
    encoded = os.path.join(tmp_dir, "sample.satz")
    if layout == "fasta":
        formats.encode_fasta(sample_path, "none", encoded, level)
    else:
        formats.encode_tsv(sample_path, "none", encoded, level, has_header=has_header)
    ratio = os.path.getsize(encoded) / written
    os.unlink(sample_path)
    os.unlink(encoded)
    return ratio


def _kept_fraction_tsv(path, compression, min_array_length):
    """Exact ``(rows_in, rows_out, content_bytes, kept_bytes)`` for a per-copy TSV."""
    rows_in = rows_out = 0
    total = kept = 0
    with formats.open_maybe_gzip(path, compression) as fh:
        header = fh.readline()
        total += len(header)
        kept += len(header)
        for raw in fh:
            rows_in += 1
            total += len(raw)
            array_id = raw.split(b"\t", 1)[0]
            length = formats.array_length_from_id(array_id)
            if length is None or length >= min_array_length:
                rows_out += 1
                kept += len(raw)
    return rows_in, rows_out, total, kept


def _kept_fraction_fasta(path, compression, min_array_length):
    rows_in = rows_out = 0
    total = kept = 0
    with formats.open_maybe_gzip(path, compression) as fh:
        for head, body in formats.iter_fasta_records(fh):
            rows_in += 1
            size = len(head) + len(body)
            total += size
            length = formats.array_length_from_fasta_header(head)
            if length is None or length >= min_array_length:
                rows_out += 1
                kept += size
    return rows_in, rows_out, total, kept


def build_plan(run_dir, config, prefix="", measure=True, probe=None):
    """Inventory *run_dir* and price every file against the policy.

    Args:
        measure: when False, skip the codec sampling and the exact row counts.
            The ledger then reports ``after`` as 0 for re-encoded files, which
            is honest about not having measured rather than filling in a guess.
    """
    files = []
    unclassified = []
    scratch = tempfile.mkdtemp(prefix="satellome_ledger_", dir=run_dir)

    try:
        return _build_plan(run_dir, config, prefix, measure, files, unclassified,
                           scratch, probe)
    finally:
        shutil.rmtree(scratch, ignore_errors=True)


def _build_plan(run_dir, config, prefix, measure, files, unclassified, scratch,
                probe):
    for rel_path in iter_run_files(run_dir):
        abs_path = os.path.join(run_dir, rel_path)
        match = classify_path(rel_path)
        stored = os.path.getsize(abs_path)

        if match.kind is None:
            unclassified.append(rel_path)
            files.append(
                PlannedFile(rel_path, "unknown", FileClass.UNKNOWN, None,
                            match.compression, stored, stored, stored)
            )
            continue

        kind = match.kind
        if kind.file_class in DROPPED_CLASSES:
            files.append(
                PlannedFile(rel_path, kind.name, kind.file_class, None,
                            match.compression, stored, stored, 0)
            )
            continue

        if kind.file_class == FileClass.KEEP or kind.layout is None:
            files.append(
                PlannedFile(rel_path, kind.name, kind.file_class, kind.layout,
                            match.compression, stored, stored, stored)
            )
            continue

        if not measure:
            files.append(
                PlannedFile(rel_path, kind.name, kind.file_class, kind.layout,
                            match.compression, stored, stored, 0)
            )
            continue

        filtered_here = kind.filtered and (probe is None or probe.can_filter(kind.name))
        if filtered_here:
            counter = (
                _kept_fraction_fasta if kind.layout == "fasta" else _kept_fraction_tsv
            )
            rows_in, rows_out, content, kept_bytes = counter(
                abs_path, match.compression, config.min_array_length
            )
        else:
            content = _content_bytes(abs_path, match.compression)
            kept_bytes = content
            rows_in = rows_out = 0

        if kind.layout in {"tsv", "fasta"}:
            ratio = _columnar_ratio(
                abs_path, match.compression, config.level, kind.layout,
                kind.has_header, scratch,
            )
        else:
            ratio = _blob_ratio(abs_path, match.compression, config.level)
        estimated = int(round(kept_bytes * ratio))
        if not filtered_here:
            # An unfiltered re-encode is only taken when it is actually smaller.
            # The engine keeps a file whose container would be no smaller — on a
            # table of a few hundred rows the per-column framing costs more than
            # the split saves. Predicting a saving it will refuse to take is how
            # a ledger ends up over-promising by a fixed, invisible margin.
            estimated = min(estimated, stored)
        files.append(
            PlannedFile(
                rel_path, kind.name, kind.file_class, kind.layout,
                match.compression, stored, content,
                estimated, rows_in, rows_out,
            )
        )

    return Plan(os.path.abspath(run_dir), prefix, files, unclassified,
                dict(probe.refused) if probe else {})


def human(n):
    """Bytes as a short human-readable string."""
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if abs(n) < 1024 or unit == "TB":
            return f"{n:.1f} {unit}" if unit != "B" else f"{int(n)} B"
        n /= 1024.0
    return f"{n:.1f} TB"


def format_ledger(plan, config, realised=False):
    """Render the per-kind byte ledger."""
    by_kind = {}
    for entry in plan.files:
        row = by_kind.setdefault(
            entry.kind_name, {"class": entry.file_class, "n": 0, "before": 0,
                              "after": 0, "rows_in": 0, "rows_out": 0}
        )
        row["n"] += 1
        row["before"] += entry.stored_bytes
        row["after"] += entry.after_bytes
        row["rows_in"] += entry.rows_in
        row["rows_out"] += entry.rows_out

    word = "actual" if realised else "estimate"
    lines = [
        f"Compaction ledger for {plan.run_dir}",
        f"  policy: min_array_length={config.min_array_length}, "
        f"codec=zstd level={config.level}",
        "",
        f"  {'kind':<16}{'class':<18}{'n':>4}{'now':>12}{'after':>12}{'freed':>12}   rows",
    ]
    for name in sorted(by_kind, key=lambda k: -by_kind[k]["before"]):
        row = by_kind[name]
        freed = row["before"] - row["after"]
        rows = (
            f"{row['rows_out']}/{row['rows_in']}"
            if row["rows_in"] else ""
        )
        lines.append(
            f"  {name:<16}{row['class']:<18}{row['n']:>4}"
            f"{human(row['before']):>12}{human(row['after']):>12}"
            f"{human(freed):>12}   {rows}"
        )
    before, after = plan.before, plan.after
    freed = before - after
    share = (100.0 * freed / before) if before else 0.0
    lines.append("")
    lines.append(
        f"  total ({word}): {human(before)} -> {human(after)}  "
        f"freed {human(freed)} ({share:.1f}%)"
    )
    if plan.refused:
        lines.append("")
        lines.append(
            "  the per-copy cut will NOT be applied to these, so they are "
            "priced as kept whole:"
        )
        for name, reason in sorted(plan.refused.items()):
            lines.append(f"    {name}: {reason}")
    if plan.unclassified:
        lines.append("")
        lines.append(
            f"  kept, unclassified ({len(plan.unclassified)}) — no policy row "
            f"matches these, so they are left exactly as they are:"
        )
        for rel in plan.unclassified[:20]:
            lines.append(f"    {rel}")
        if len(plan.unclassified) > 20:
            lines.append(f"    ... and {len(plan.unclassified) - 20} more")
    return lines
