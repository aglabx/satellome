"""``satellome compact`` and ``satellome expand``.

These are subcommands rather than flags because they act on a *finished*
directory and take none of the pipeline's options. Dispatch happens before the
main flag parser, so nothing about the existing command line changes.
"""

import argparse
import logging
import os
import time

from satellome.compact import engine, formats
from satellome.compact.ledger import human
from satellome.compact.policy import (
    KINDS,
    REENCODED_CLASSES,
    CompactConfig,
    classify_path,
)

logger = logging.getLogger("satellome")

SUBCOMMANDS = ("compact", "expand")


def build_compact_parser():
    parser = argparse.ArgumentParser(
        prog="satellome compact",
        description=(
            "Reduce finished satellome output directories to their primary "
            "layer. Every dropped file is regenerated and compared before it is "
            "removed, and .compact.json records how to put it back."
        ),
    )
    parser.add_argument("directories", nargs="*", metavar="DIR",
                        help="run directories to compact")
    parser.add_argument("--dry-run", action="store_true",
                        help="print the byte ledger and write nothing")
    parser.add_argument("--min-array-length", type=int, default=1000, metavar="BP",
                        help="keep per-copy rows only for arrays at least this "
                             "long [1000]; this is a scope decision, not a "
                             "storage one — see the policy notes")
    parser.add_argument("--level", type=int, default=15, metavar="N",
                        help="zstd level for the columnar codec [15]; measure "
                             "with --sweep before changing it")
    parser.add_argument("--no-verify-drops", dest="verify_drops",
                        action="store_false", default=True,
                        help="do not regenerate each dropped file to check it "
                             "before unlinking (much faster; only after "
                             "--check-recipes has passed on a sample)")
    parser.add_argument("--check-recipes", action="store_true",
                        help="Checklist 0: report whether every recipe "
                             "reproduces the file it would replace, and change "
                             "nothing")
    parser.add_argument("--sweep", metavar="LEVELS", default=None,
                        help="measure the codec at these levels (e.g. "
                             "12,15,19) on this directory's tables and exit")
    parser.add_argument("--explain", action="store_true",
                        help="print the policy table and exit")
    parser.add_argument("--from-file", metavar="PATH", default=None,
                        help="read run directories from a file, one per line")
    parser.add_argument("--continue-on-error", action="store_true",
                        help="keep going after a directory is refused")
    return parser


def build_expand_parser():
    parser = argparse.ArgumentParser(
        prog="satellome expand",
        description=(
            "Restore compacted directories to their pre-compaction contents, "
            "verifying every restored file against the md5 recorded before "
            "compaction."
        ),
    )
    parser.add_argument("directories", nargs="*", metavar="DIR")
    parser.add_argument("--from-file", metavar="PATH", default=None,
                        help="read run directories from a file, one per line")
    parser.add_argument("--threads", type=int, default=1, metavar="N",
                        help="threads for the decomposer [1]; the content is "
                             "thread-count independent but the row order is "
                             "not, so the default reproduces the original")
    parser.add_argument("--continue-on-error", action="store_true")
    return parser


def _directories(args):
    directories = list(args.directories)
    if args.from_file:
        try:
            with open(args.from_file, "r") as fh:
                directories.extend(
                    line.strip() for line in fh if line.strip() and not line.startswith("#")
                )
        except OSError as e:
            logger.error(f"--from-file {args.from_file}: {e}")
            return None
    return directories


def explain_policy():
    lines = ["Compaction policy", ""]
    width = max(len(k.suffix) for k in KINDS)
    for kind in KINDS:
        lines.append(f"  {kind.suffix:<{width}}  {kind.name:<16}{kind.file_class}")
        if kind.why:
            lines.append(f"  {'':<{width}}  {kind.why}")
    return lines


def sweep_codec(run_dir, levels, config):
    """Measure the codec at several levels on this directory's own tables.

    The level is a parameter to measure, not to inherit: on real tables the
    difference between 12 and 19 is 1.7x in time for well under 10% in size, and
    which end of that trade is right depends on how much of the corpus is left
    to process.
    """
    lines = [f"Codec sweep on {run_dir}", ""]
    lines.append(f"  {'file':<40}{'level':>6}{'bytes':>12}{'vs now':>9}{'sec':>8}")
    for rel_path in sorted(os.listdir(os.path.join(run_dir, "fastan"))
                           if os.path.isdir(os.path.join(run_dir, "fastan")) else []):
        rel = os.path.join("fastan", rel_path)
        _sweep_one(run_dir, rel, levels, config, lines)
    for name in sorted(os.listdir(run_dir)):
        _sweep_one(run_dir, name, levels, config, lines)
    return lines


def _sweep_one(run_dir, rel, levels, config, lines):
    abs_path = os.path.join(run_dir, rel)
    if not os.path.isfile(abs_path):
        return
    match = classify_path(rel)
    if match.kind is None or match.kind.file_class not in REENCODED_CLASSES:
        return
    if match.kind.layout not in {"tsv", "fasta"}:
        return
    now = os.path.getsize(abs_path)
    for level in levels:
        out = abs_path + f".sweep{level}"
        start = time.time()
        try:
            if match.kind.layout == "fasta":
                formats.encode_fasta(abs_path, match.compression, out, level)
            else:
                formats.encode_tsv(abs_path, match.compression, out, level,
                                   has_header=(match.kind.name != "bed"))
            size = os.path.getsize(out)
            lines.append(
                f"  {rel[:40]:<40}{level:>6}{size:>12}"
                f"{now / size:>8.2f}x{time.time() - start:>8.1f}"
            )
        finally:
            if os.path.exists(out):
                os.unlink(out)


def run_compact(argv):
    args = build_compact_parser().parse_args(argv)

    if args.explain:
        for line in explain_policy():
            print(line)
        return 0

    directories = _directories(args)
    if directories is None:
        return 2
    if not directories:
        logger.error("satellome compact needs at least one directory (or --from-file)")
        return 2

    config = CompactConfig(
        min_array_length=args.min_array_length,
        level=args.level,
        verify_drops=args.verify_drops,
    )

    if args.sweep:
        try:
            levels = [int(x) for x in args.sweep.split(",") if x.strip()]
        except ValueError:
            logger.error(f"--sweep expects comma-separated integers, got {args.sweep!r}")
            return 2
        for directory in directories:
            for line in sweep_codec(directory, levels, config):
                logger.info(line)
        return 0

    failures = 0
    total_before = total_after = 0
    for directory in directories:
        if args.check_recipes:
            outcome = engine.check_recipes(directory, config, log=logger)
        else:
            outcome = engine.compact(directory, config, dry_run=args.dry_run, log=logger)
        for line in outcome.lines:
            (logger.info if outcome.ok else logger.error)(line)
        total_before += outcome.before
        total_after += outcome.after
        if not outcome.ok:
            failures += 1
            if not args.continue_on_error:
                return 1

    if len(directories) > 1:
        freed = total_before - total_after
        logger.info(
            f"{len(directories) - failures} of {len(directories)} directories "
            f"processed; {human(total_before)} -> {human(total_after)} "
            f"(freed {human(freed)})"
        )
    return 1 if failures else 0


def run_expand(argv):
    args = build_expand_parser().parse_args(argv)
    directories = _directories(args)
    if directories is None:
        return 2
    if not directories:
        logger.error("satellome expand needs at least one directory (or --from-file)")
        return 2

    failures = 0
    for directory in directories:
        outcome = engine.expand(directory, log=logger, threads=args.threads)
        for line in outcome.lines:
            (logger.info if outcome.ok else logger.error)(line)
        if not outcome.ok:
            failures += 1
            if not args.continue_on_error:
                return 1
    return 1 if failures else 0


def dispatch(argv):
    """Run a compact/expand subcommand. Returns the exit code, or None if not ours."""
    if not argv or argv[0] not in SUBCOMMANDS:
        return None
    if argv[0] == "compact":
        return run_compact(argv[1:])
    return run_expand(argv[1:])
