# `satellome compact` / `satellome expand`

Reduce a finished output directory to the layer that cannot be recomputed from
it, and put the rest back on demand.

## Why

The eukaryotic satellite panel produces about **190 MB per assembly**. The roster
holds **26,266** assemblies, so the finished corpus is roughly **4.9 TB**, and
the volume that has to hold it has 1.4 TB free. The runners carry disk floors so
that they refuse rather than write a truncated catalogue — which means the panel
does not slow down as the volume fills, it *stops*, at roughly 77% of the roster.

With the policy below the same 26,266 assemblies occupy about **1.0 TB**.

Two halves, and both are needed:

* `satellome compact` retrofits the directories already on disk.
* Runs without `--extended-output` no longer produce the droppable files at all,
  so the remaining assemblies land compacted.

## Using it

```bash
# What would happen, priced, without touching anything
satellome compact --dry-run /mnt/data/analysis/satellome/vgp/GCA_000001405.29

# Checklist 0: does every recipe reproduce the file it would replace?
satellome compact --check-recipes <dir>

# Do it
satellome compact <dir> [<dir> ...]
satellome compact --from-file dirs.txt --continue-on-error

# Put everything back
satellome expand <dir>

# The policy table itself, with the reason for every row
satellome compact --explain
```

Exit codes: `0` done (or nothing to do), `1` at least one directory refused or
failed, `2` bad arguments.

## What happens to each kind

`satellome compact --explain` prints the live table. In summary:

| kind | class | action |
|---|---|---|
| `<p>.sat` | primary | keep, re-encode columnar |
| `<p>.10kb.sat` | primary (working view) | keep, re-encode |
| `fastan/<p>.summary.tsv` | primary | keep, re-encode |
| `fastan/<p>.bed` | primary | keep, re-encode |
| `fastan/<p>.monomers.tsv` | primary, filtered | keep rows of arrays ≥ 1 kb |
| `fastan/<p>.hors.tsv` | primary, filtered | keep rows of arrays ≥ 1 kb |
| `fastan/<p>.decomposed.fasta` | primary, filtered | keep records of arrays ≥ 1 kb |
| `fastan/<p>.lengths` | projection | drop — regenerate from `decomposed` |
| `<p>.{1kb,100kb,1000kb}.sat` | view | drop — length filter on the master |
| `<p>.{micro,pmicro,tssr,complex}.sat` | view | drop — classifier over the master |
| `gff3/*.gff` | not needed | drop — fSSR is an early-version atavism |
| `fastan/<p>.1aln` | not needed | drop — **the one kind expand cannot rebuild offline** |
| `reports/`, `*.html`, logs | primary | keep, re-encode whole |
| anything else | unknown | **kept, and named in the report** |

Classification is on the **trailing suffix**, never on a prefix: the basename
comes from the input FASTA, so directory `GCA_029289425.3` legitimately holds
`GCF_029289425.2_NHGRI_mPanPan1-v2.0_pri_genomic.*`.

## The guarantees, and their exact edges

**Nothing is removed that has not just been reproduced.** Before a file is
unlinked it is regenerated into scratch and compared by md5 with the file about
to go. Before a re-encoded file replaces its original, the container is decoded
and checked against its own footer. "Reversible" is a fact established per file,
not a property claimed for the design.

**Content-exact, not container-exact.** A restored `.gz` holds exactly the
original bytes, but the gzip *framing* around them is re-encoded: a deflate
stream cannot be reproduced without the compressor that made it, and this corpus
was gzipped by a tool that is not this Python. The embedded original file name
and mtime are preserved; the compressed bytes are not. Compare with `zcat | md5`,
not `md5`.

**A failure narrows the compaction; it never widens the damage.** A kind whose
rows cannot be shown recomputable here is re-encoded whole — no rows dropped —
and the reason is printed and written into `.compact.json`.

**Idempotent.** A directory carrying `.compact.json` is a no-op: zero writes.

**`satellome --verify-run` keeps working.** Compaction and expansion both
refresh `run_manifest.json` and append to its `compaction_history`, so the
driver's completeness check stays meaningful at every point.

## What it refuses, and why the obvious guards are wrong

Checked on the live tree, not assumed:

* **`.rc_complete` is not a completeness marker.** 1,020 finished catalogues do
  not carry it — they predate it. Gating on it would exclude them forever.
* **A `fasta/` subdirectory is not an in-flight signal.** 1,015 tree directories
  carry one left from an earlier campaign.

What actually marks a directory under a live worker is an **uncompressed `.sat`
sitting beside its `.gz`** — compression happens at the end. The positive test
is: a master `.sat`/`.sat.gz` that decompresses to its last byte, a
`fastan/*.bed`, a `results.yaml`, no `.partial` files, and no uncompressed `.sat`
beside a compressed one.

Above all, compaction runs **only on the finished tree, never on `_incoming/`**.
That is a guarantee by construction; a check can be written wrong.

Every refusal is a named line in the report. A silent skip is indistinguishable
from success.

## The 1 kb threshold

`--min-array-length` defaults to 1000. **This is a scope decision, not a storage
one**, and must not be presented as if it were: the fraction it frees ranges from
19% to 92% across assemblies, so no storage argument can fix a value that moves
by 5×. It comes from the tool being for satellites rather than for tandem repeats
in general, and it matches the boundary at which `1kb.sat` was already cut.

What the cut does **not** remove is the sub-kilobase *arrays*. They keep their
master row — coordinates, period, copy number, consensus, entropy, GC, and the
full array sequence. What goes is their per-copy monomer decomposition, which is
recomputable from `trf_array`.

## The decomposer, and the one thing to know before restoring an old catalogue

The per-copy layer can only be rebuilt by the **arraysplitter build that
produced it**. Measured on `GCA_022385595.1`: a locally installed arraysplitter
0.1.0 writes 16 columns instead of 17, emits no `hors.tsv` at all, and decomposes
**3.4%** of the arrays differently — a different period, a different monomer
count. Point at the right one with:

```bash
export SATELLOME_ARRAYSPLITTER=/path/to/the/arraysplitter/that/made/this/corpus
```

`--dry-run` probes the decomposer over 200 sub-threshold arrays before pricing
anything, so a ledger produced on a machine with the wrong build says so instead
of promising a saving it will not deliver.

What *is* established about the decomposer, on real data:

* **deterministic** — three runs at `-t 4` and two at `-t 1` produce byte-identical output;
* **per-array independent** — run over only the 13,233 sub-kilobase arrays it
  reproduces their rows byte-identically, 13,233 of 13,233, against the full run;
* the corpus files are in **master order**, with each array's rows contiguous,
  which is why re-inserting regenerated rows needs no stored positions.

## The columnar codec

A row of `monomers.tsv` interleaves seventeen value types and gzip's 32 KB window
sees them mixed. Split by column and each stream is homogeneous: `array_id`
collapses from 45.6 MB to 0.9 MB. Measured against whole-file compression at the
same level, on six real tables: **1.98× byte-weighted**, and the gain grows with
table size. `--long` windows buy 1%; the split buys 2×.

Against the corpus's existing gzip, at level 15, on `GCA_022385595.1`:

| file | now | after | gain |
|---|---:|---:|---:|
| `monomers.tsv.gz` | 5.0 MB | 2.7 MB | 1.87× |
| `bed` | 402 KB | 90 KB | 4.0× |
| `sat.gz` | 1.3 MB | 1.0 MB | 1.35× |

**The level is a parameter to measure, not to inherit.** `--level` defaults to 15;
`satellome compact --sweep 12,15,19 <dir>` measures it on the tables in front of
you. Level 19 costs 2–3× the time for under 10% more.

A file the container would make *larger* is kept as it is — on a table of a few
hundred rows the per-column framing costs more than the split saves — and the
ledger prices it that way, so it never promises a saving the engine will refuse.

## `.compact.json`

Written into every compacted directory. Per file it records the class, the
action, the md5 and byte length of the content **before**, the same **after**,
the gzip framing, and the recipe. `expand` reads nothing else. Without this
record, compaction is deletion.

## Runs that never produce the files

Every satellome run now prunes the droppable kinds at the end, before writing its
manifest, and says what it removed:

```
Compact output: removed 12 derived file(s), 2.4 MB freed (gff x9, lengths x1, aln x1, sat_micro x1).
Run with --extended-output to keep them, or 'satellome expand <dir>' after a compaction to rebuild them.
```

`--extended-output` keeps everything, as before. `<p>.1kb.sat` and
`<p>.10kb.sat` are kept regardless: the drawing and annotation steps read them,
and so does `--rerun drawing`.

## Retrofitting the corpus

```bash
# 1. Prove the recipes on a stratified sample first. Nothing is deleted
#    corpus-wide until this passes on every tier.
for d in $(shuf -n 100 dirs.txt); do satellome compact --check-recipes "$d"; done

# 2. Price it.
satellome compact --dry-run --from-file dirs.txt | tee ledger.txt

# 3. Do it.
satellome compact --from-file dirs.txt --continue-on-error
```

By default every drop is re-derived and compared before it is unlinked, which
means a full decomposer run per directory. That is the right default for one
directory and too slow for 15,620 of them, so once step 1 has passed on a
stratified sample, add `--no-verify-drops`:

```bash
satellome compact --from-file dirs.txt --continue-on-error --no-verify-drops
```

That skips the per-file re-derivation, **not** the checks: the 200-array probe
still runs on every directory, the readiness test still refuses anything that is
not verifiably finished, containers are still decoded and compared to their own
footers, and `.compact.json` still records every digest. What it gives up is the
proof that *these particular* sub-threshold rows come back — which is exactly
what step 1 established for the corpus as a whole. The record says which files
were dropped unverified, so the distinction is never lost.

A T2T assembly frees only about 19% on the per-copy cut, and that is correct —
its arrays are long and intact, so most copy rows belong to arrays above the
threshold. A T2T assembly reporting a *high* freed fraction means the assembly or
its catalogue is wrong, not that compaction worked well.
