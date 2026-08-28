# ADR 0001 — Rewriting the Python engine in Rust

**Status:** deferred (measured, not acted on)
**Date:** 2026-08-28
**Context version:** satellome 1.14.0

## Question

Can the engine be rewritten from Python to Rust, preserving logic and
compatibility?

## What the "engine" actually is today

| | Lines | What it does |
|---|---:|---|
| Python (`src/satellome`) | 17 491 | orchestration, parsing, classification, drawing, reports, installers |
| Rust (`rust/`) | 1 132 | sat-family, telomere-check, find-gaps, bed-extract, genome-size |
| External native | — | FasTAN, tanbed, trf, arraysplitter |

There are 18 `subprocess` call sites outside the installers. **Every
compute-heavy stage is already native**: the tandem-repeat search, monomer and
array extraction, family clustering, telomere detection, gap finding and genome
size. Python is the glue around them plus the handling of their results.

So the question is not whether to rewrite the engine — it is whether to rewrite
what is left.

## Measurements

On a synthetic `.sat` of 300 000 records (79 MB), Python 3.11, M-series laptop:

```
raw line read            0.09 s    868 MB/s     <- the floor
parse into TRModel       2.97 s     26 MB/s     33x off the floor, 9.9 us/record
parse + build UCSC track 3.26 s
```

`cProfile` of the parse:

```
tottime  calls       function
  2.781  300 000     AbstractModel.set_with_dict
  0.713  300 001     AbstractModel.__init__
  0.637  14 100 029  builtins.setattr          <- 47 fields per record
  0.447  6 000 000   set_with_dict.<genexpr>   <- legacy detection, per record
  0.257  5 400 000   preprocess_pair           <- a no-op hook, called per field
```

`set_with_dict` recomputes the legacy-format field set on **every record**, and
tests membership against `int_attributes` / `float_attributes` /
`list_attributes`, which are lists, not sets.

### The decisive measurement

A prototype of an *optimal* Python parser — column plan resolved once,
`__slots__`, frozensets, no per-record recomputation:

```
optimised python parse   2.12 s     37 MB/s     7.1 us/record  ->  only 1.4x faster
```

**Python cannot be fixed into being fast here.** 7.1 µs across 47 fields is
~150 ns per field: the cost is the interpreter's per-attribute overhead, not
sloppy code. Even an ideal Python parser stays 24× off the I/O floor. A Rust
TSV parser sustains 200–500 MB/s, i.e. 20–40× the current implementation.

For a T2T assembly (~5 M repeats, ~1.3 GB `.sat`) that is ~50 s per pass in
Python versus ~2 s in Rust, and the pipeline makes several passes.

## What is not measured

**Python's share of a real run.** If FasTAN takes two hours on a 3 Gb genome,
four minutes of parsing is noise and none of this is worth doing. If FasTAN
takes ten minutes, Python is a third of the wall clock. Without that number any
decision is a guess.

Getting it is cheap: record per-step wall time into `run_manifest.json` during a
normal run. That is worth having regardless of this decision.

## Decision

**Do not rewrite wholesale. Do not start now.** Revisit when a real run shows
Python to be a meaningful share of it.

A full rewrite is the wrong shape for three reasons:

1. **Drawing and reports (~2–3k lines)** are plotly/matplotlib. Rust has no
   equivalent; this would trade functionality for speed that these stages do
   not need.
2. **Orchestration (`main.py`, ~2k lines)** — installers, the directory lock,
   the run manifest, `--rerun`, `--doctor`. No performance to gain, maximum
   behavioural risk.
3. **The file formats are the compatibility contract** and are already stable
   and tested: `.sat` TSV, GFF3, BED, `run_manifest.json`, CLI flags, exit
   codes. A rewrite must not touch them, which means it buys nothing on that
   axis either.

## If and when it is done: the shape

Extract the **data plane** — `.sat`/TRF parsing, classification, annotation
intersection, statistics (~4–5k lines) — into a `satellome-core` crate exposed
through PyO3. CLI, orchestration, drawing and reports stay in Python.

This is the pattern the project already uses successfully: `genome-size`,
`find-gaps` and `bed-extract` are native, invoked from Python, with a Python
fallback when the binary is absent.

Compatibility is verifiable rather than hoped for:

* the 568-test suite as the behavioural net;
* differential testing — run both implementations over the same input and
  compare the output files byte for byte;
* the existing `--verify-run` manifest check, which already asserts that a run
  produced exactly the files it claims, at the sizes it claims.

Order of extraction, hottest and most isolated first:

1. `.sat` reader/writer and the record model
2. classification
3. annotation intersection
4. statistics

Each step keeps its Python implementation as a fallback until the Rust one has
run clean on real genomes.

## Revisit when

* per-step timings show Python above ~20% of a real run, or
* a genome large enough that the parse passes cost more than the search, or
* a new stage is needed that is intrinsically per-record and hot.
