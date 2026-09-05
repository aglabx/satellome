"""A synthetic satellome run directory, self-consistent by construction.

The per-copy tables are produced by the same stub decomposer the recipes will
call to rebuild them, which is exactly the production situation: a catalogue can
only be restored by the decomposer that made it. The derived views and the gff
are produced by the real classifier, so the fixture exercises the real recipe
rather than a mock of it.
"""

import gzip
import os
import stat
import sys
import textwrap

MASTER_HEADER = (
    "project\ttrf_id\ttrf_head\ttrf_l_ind\ttrf_r_ind\ttrf_period\ttrf_n_copy\t"
    "trf_pmatch\ttrf_pvar\ttrf_entropy\ttrf_consensus\ttrf_array\ttrf_array_gc\t"
    "trf_consensus_gc\ttrf_array_length\ttrf_joined\ttrf_family\ttrf_ref_annotation\n"
)

PREAMBLE = (
    "# FasTAN results converted to TRF format\n"
    "# Source FASTA: {prefix}.fna\n"
    "# Source BED: {prefix}.bed\n"
    "# Note: pmatch, pvar, entropy are calculated from sequence data\n"
)

#: A decomposer that is deterministic and per-array independent, like the real
#: one, but small enough to live in a test.
STUB_DECOMPOSER = textwrap.dedent(
    '''
    import sys

    args = sys.argv[1:]
    opts = dict(zip(args[::2], args[1::2]))
    records = []
    header = None
    body = []
    for line in open(opts["-i"]):
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(body)))
            header = line[1:].strip()
            body = []
        else:
            body.append(line.strip())
    if header is not None:
        records.append((header, "".join(body)))

    prefix = opts["-o"]
    mono = open(prefix + ".monomers.tsv", "w")
    hors = open(prefix + ".hors.tsv", "w")
    dec = open(prefix + ".decomposed.fasta", "w")
    lengths = open(prefix + ".lengths", "w")
    mono.write("array_id\\ttype\\tidx\\tlength\\tsequence\\n")
    hors.write("array_id\\ttype\\tidx\\tlength\\tsequence\\n")
    for array_id, seq in records:
        period = int(array_id.rsplit("_", 1)[1])
        pieces = [seq[i:i + period] for i in range(0, len(seq), period)] or [""]
        for idx, piece in enumerate(pieces):
            mono.write("%s\\tmono\\t%d\\t%d\\t%s\\n" % (array_id, idx, len(piece), piece))
        hors.write("%s\\thor\\t0\\t%d\\t%s\\n" % (array_id, len(seq), pieces[0]))
        head = ">%s n_monomers=%d\\n" % (array_id, len(pieces))
        dec.write(head + " ".join(pieces) + "\\n")
        lengths.write(head + " ".join(str(len(p)) for p in pieces) + "\\n")
    for fh in (mono, hors, dec, lengths):
        fh.close()
    '''
).strip()


def write_stub_decomposer(tmp_path):
    """Write the stub and return a path usable as $SATELLOME_ARRAYSPLITTER."""
    script = tmp_path / "stub_arraysplitter.py"
    script.write_text(STUB_DECOMPOSER + "\n")
    launcher = tmp_path / "stub_arraysplitter"
    launcher.write_text(f'#!/bin/sh\nexec "{sys.executable}" "{script}" "$@"\n')
    launcher.chmod(launcher.stat().st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)
    return str(launcher)


def _array(seed, length):
    unit = "ACGT" * 8
    return (unit * (length // len(unit) + 1))[:length]


def master_rows(prefix, arrays):
    """arrays: list of (chrom, start, length, period)."""
    rows = []
    for i, (chrom, start, length, period) in enumerate(arrays, start=1):
        seq = _array(i, length)
        rows.append(
            "\t".join([
                prefix, str(i), chrom, str(start + 1), str(start + length),
                str(period), f"{length / period:.1f}", "95", "5", "1.9",
                seq[:period], seq, "0.50", "0.50", str(length), "0", "", "",
            ])
        )
    return rows


def build_run_dir(root, prefix="TESTASM_v1", arrays=None, gzip_tables=True,
                  with_manifest=True, decomposer=None):
    """Create a finished-looking run directory and return its path."""
    if arrays is None:
        # Five named arrays the tests reason about by hand, then enough filler
        # that the tables are big enough for columnar encoding to be worth
        # doing. A handful of rows is a case where the container is legitimately
        # larger than the gzip, and testing only that case would test the guard
        # instead of the codec.
        arrays = [
            ("NC_073249.2", 100, 2400, 12),
            ("NC_073249.2", 5000, 480, 8),
            ("ABRL03000001.1", 0, 1600, 16),
            ("CM002639.2", 77, 300, 10),
            ("CM002639.2", 900, 12000, 24),
        ]
        position = 100_000
        for i in range(400):
            length = 200 + (i * 37) % 5000
            arrays.append((f"scaffold_{i % 7}.1", position, length, 6 + i % 30))
            position += length + 500
    run_dir = os.path.join(str(root), "RUN_DIR")
    os.makedirs(os.path.join(run_dir, "fastan"), exist_ok=True)

    rows = master_rows(prefix, arrays)
    master_text = PREAMBLE.format(prefix=prefix) + MASTER_HEADER + "".join(
        row + "\n" for row in rows
    )
    master_plain = os.path.join(run_dir, f"{prefix}.sat")
    with open(master_plain, "w") as fh:
        fh.write(master_text)

    # fastan/*.bed — headerless, five columns.
    with open(os.path.join(run_dir, "fastan", f"{prefix}.bed"), "w") as fh:
        for chrom, start, length, period in arrays:
            fh.write(f"{chrom}\t{start}\t{start + length}\t{period}\t900\n")

    # results.yaml, which carries the genome size every recipe needs.
    genome_size = sum(a[2] for a in arrays) * 40
    with open(os.path.join(run_dir, "results.yaml"), "w") as fh:
        fh.write(
            "pid: project\n"
            "work_files:\n"
            "  assembly_stats:\n"
            "    dataset:\n"
            f"      genome_size: {genome_size}\n"
            "  ref_assembly_name_for_trf: dataset\n"
        )

    # The derived views and gff come from the real classifier.
    from satellome.steps.trf_classify import classify_trf_data

    classify_trf_data(os.path.join(run_dir, prefix), run_dir, genome_size)

    # The per-copy layer comes from the stub decomposer, run over every array.
    import subprocess

    arrays_fasta = os.path.join(run_dir, "arrays.fasta")
    with open(arrays_fasta, "w") as fh:
        for i, (chrom, start, length, period) in enumerate(arrays, start=1):
            fh.write(f">{chrom}_{start}_{start + length}_{length}_{period}\n")
            fh.write(_array(i, length) + "\n")
    subprocess.run(
        [decomposer, "-i", arrays_fasta, "-o", os.path.join(run_dir, "fastan", prefix)],
        check=True,
    )
    os.unlink(arrays_fasta)

    # summary.tsv, one row per array.
    with open(os.path.join(run_dir, "fastan", f"{prefix}.summary.tsv"), "w") as fh:
        fh.write("array_id\tarray_length\torientation\n")
        for chrom, start, length, period in arrays:
            fh.write(f"{chrom}_{start}_{start + length}_{length}_{period}\t{length}\tfwd\n")

    # The alignment intermediate: dropped by policy, not restorable offline.
    with open(os.path.join(run_dir, "fastan", f"{prefix}.1aln"), "wb") as fh:
        fh.write(b"\x00binary alignment intermediate\x00" * 100)

    with open(os.path.join(run_dir, "_satellome.log"), "w") as fh:
        fh.write("2026-09-05 00:00:00 - satellome - INFO - run finished\n" * 40)

    if gzip_tables:
        for rel in [
            f"{prefix}.sat", f"{prefix}.1kb.sat", f"{prefix}.10kb.sat",
            f"{prefix}.100kb.sat", f"{prefix}.1000kb.sat", f"{prefix}.micro.sat",
            f"{prefix}.pmicro.sat", f"{prefix}.tssr.sat", f"{prefix}.complex.sat",
            f"fastan/{prefix}.monomers.tsv", f"fastan/{prefix}.hors.tsv",
            f"fastan/{prefix}.decomposed.fasta", f"fastan/{prefix}.summary.tsv",
        ]:
            path = os.path.join(run_dir, rel)
            if not os.path.exists(path):
                continue
            with open(path, "rb") as src:
                data = src.read()
            with gzip.GzipFile(filename=os.path.basename(rel), mode="wb",
                               fileobj=open(path + ".gz", "wb"), mtime=1700000000) as out:
                out.write(data)
            os.unlink(path)

    # fasta/ is left over from an earlier campaign in 1,015 tree directories and
    # must not be read as "still computing".
    os.makedirs(os.path.join(run_dir, "fasta"), exist_ok=True)

    if with_manifest:
        from satellome.core_functions.tools.run_manifest import (
            build_manifest,
            write_run_manifest,
        )

        manifest = build_manifest(run_dir, project=prefix, version="test",
                                 steps={"fastan": "ok", "classification": "ok"})
        write_run_manifest(run_dir, manifest)

    return run_dir


def content_md5(path):
    import hashlib

    opener = gzip.open if path.endswith(".gz") else open
    digest = hashlib.md5()
    with opener(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def snapshot(run_dir):
    """{relative path: content md5} for every file under *run_dir*."""
    out = {}
    for root, _dirs, files in os.walk(run_dir):
        for name in files:
            path = os.path.join(root, name)
            out[os.path.relpath(path, run_dir)] = content_md5(path)
    return out
