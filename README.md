# Satellome

[![Tests](https://github.com/aglabx/satellome/workflows/Tests/badge.svg)](https://github.com/aglabx/satellome/actions?query=workflow%3ATests)
[![codecov](https://codecov.io/gh/aglabx/satellome/branch/main/graph/badge.svg)](https://codecov.io/gh/aglabx/satellome)
[![Python Version](https://img.shields.io/badge/python-3.9%20%7C%203.10%20%7C%203.11-blue)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![PyPI version](https://badge.fury.io/py/satellome.svg)](https://badge.fury.io/py/satellome)

A comprehensive bioinformatics tool for analyzing satellite DNA (tandem repeats) in telomere-to-telomere (T2T) genome assemblies.

## Overview

Satellome uses **FasTAN** (Fast Tandem Repeat Finder) as its default tandem repeat detection engine, providing fast and accurate identification of repetitive DNA sequences. The tool classifies and visualizes tandem repeats with a focus on centromeric and telomeric regions.

The tool is designed to work with various genome assembly projects including:
- T2T (Telomere-to-Telomere) Consortium assemblies
- DNA Zoo chromosome-length assemblies
- VGP (Vertebrate Genome Project) assemblies
- NCBI RefSeq and GenBank assemblies

## Features

- **Fast Tandem Repeat Detection**: Uses FasTAN by default for rapid, accurate detection
- **Smart Classification**: Categorizes repeats into microsatellites, complex repeats, and other types
- **Rich Visualizations**: Generates karyotype plots and chromosome-level visualizations
- **Annotation Integration**: Supports GFF3 and RepeatMasker annotations
- **Parallel Processing**: Efficient handling of large genomes
- **Smart Pipeline**: Automatically skips completed steps (override with `--force`)
- **Compressed File Support**: Direct processing of .gz compressed FASTA files
- **Optional TRF Support**: Traditional TRF analysis available with `--run-trf` flag

## Quick Start

```bash
# Install from PyPI
pip install satellome

# Install required binaries (FasTAN, tanbed)
satellome --install-all

# Run on a genome
satellome -i genome.fasta -o output_dir -p project_name -t 8
```

## Installation

### From PyPI (Recommended)

```bash
pip install satellome

# Install external tools
satellome --install-all
```

### Checking the installation

```bash
satellome --doctor
```

`--doctor` reports where this install lives, whether your shell can actually
see the launcher, and where each external tool resolves from. It exits 0 when
everything is healthy and 1 when it finds a problem, so it can be used as a
setup gate in a driver script.

It is the answer to the most common post-install surprise:

```
WARNING: The script satellome is installed in '/home/user/.local/bin'
which is not on PATH.
```

That means `pip` used a user install and your shell cannot see the launcher —
`satellome` will report `command not found` even though the package is fine.
Any one of these fixes it:

```bash
export PATH="$HOME/.local/bin:$PATH"                       # this shell only
echo 'export PATH="$HOME/.local/bin:$PATH"' >> ~/.bashrc   # permanent
python -m satellome ...                                    # no PATH change needed
```

`python -m satellome` is the most robust form in batch scripts: it always
follows the active environment, whereas a console script keeps the interpreter
its shebang was written with. Satellome resolves pip-installed companion tools
(such as `arraysplitter`) through the same interpreter-aware lookup, so a
directory missing from `PATH` degrades nothing silently — it is reported and
the tool is still used.

Set `SATELLOME_NO_ENV_CHECK=1` to silence the startup warning on machines where
this layout is deliberate; `--doctor` still reports it.

### From Source

```bash
git clone https://github.com/aglabx/satellome.git
cd satellome
pip install -e .
satellome --install-all
```

### External Tools

Satellome requires **FasTAN** and **tanbed** for default operation. Install them automatically:

```bash
# Install all tools (FasTAN, tanbed, modified TRF)
satellome --install-all

# Or install individually
satellome --install-fastan
satellome --install-tanbed
satellome --install-trf-large  # For genomes with chromosomes >2GB
```

**Build requirements:** git, make, C compiler (gcc/clang)

```bash
# Ubuntu/Debian
sudo apt-get install build-essential git

# macOS
xcode-select --install
```

## Usage

### Basic Command

```bash
satellome -i genome.fasta -o output_dir -p project_name -t 8
```

### Common Options

```bash
# With GFF3 annotations
satellome -i genome.fasta -o output_dir -p project_name -t 8 --gff annotations.gff3

# With RepeatMasker annotations
satellome -i genome.fasta -o output_dir -p project_name -t 8 --rm repeatmasker.out

# Force rerun all steps
satellome -i genome.fasta -o output_dir -p project_name -t 8 --force

# Also run traditional TRF analysis
satellome -i genome.fasta -o output_dir -p project_name -t 8 --run-trf
```

### Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `-i, --input` | Input FASTA file (.fa, .fasta, .gz) | Required |
| `-o, --output` | Output directory | Required |
| `-p, --project` | Project name | Required |
| `-t, --threads` | Number of threads | 1 |
| `--gff` | GFF3 annotation file | None |
| `--rm` | RepeatMasker output file | None |
| `--run-trf` | Also run TRF analysis | False |
| `--force` | Force rerun all steps | False |
| `--taxid` | NCBI taxonomy ID | None |

## Output Structure

```
output_dir/
├── genome.sat                    # Main SAT output (all arrays)
├── genome.1kb.sat                # Arrays >1kb
├── genome.3kb.sat                # Arrays >3kb
├── genome.10kb.sat               # Arrays >10kb
├── genome.micro.sat              # Microsatellites (1-9 bp monomers)
├── genome.complex.sat            # Complex repeats (>9 bp monomers)
├── genome.pmicro.sat             # Potential microsatellites
├── genome.tssr.sat               # Tandem simple sequence repeats
├── genome.gaps.bed               # Gaps annotation
├── results.yaml                  # Analysis statistics
├── run_manifest.json             # What the run produced (written last) + step statuses
├── fastan/                       # FasTAN intermediate files
│   ├── genome.1aln               # FasTAN alignment output
│   └── genome.bed                # FasTAN BED format
├── fasta/                        # FASTA sequences
│   └── genome.arrays.fasta       # All array sequences
├── gff3/                         # GFF3 annotations
│   ├── genome.1kb.gff
│   ├── genome.complex.gff
│   └── ...
├── images/                       # Visualizations
│   └── *.png
└── reports/                      # HTML reports
    └── satellome_report.html
```

### Taxon Names in File Names

Karyotype charts are named after the taxon (`--taxon`, or the name resolved from
`--taxid`). Organism names are free text and NCBI strain designations often
contain slashes — `Leishmania braziliensis MHOM/BR/75/M2904`, `Chlorella
vulgaris CCAP 1055/1`. Such a name is reduced to a single safe file-name
component (`..._MHOM_BR_75_M2904.karyo.*`) and the substitution is logged; the
plot titles keep the original name. Without this the slash was read as a path
separator and the drawing step died with `FileNotFoundError` after the whole
pipeline had already written its data.

## Verifying a Run

Do not decide that an output directory is complete by checking that some files
exist. A file that was read or copied while satellome was still writing it exists
just as hard as a complete one, and a `gzip` of such a partial read is a valid
archive that `gzip -t` accepts — so truncated data can enter downstream analysis
unnoticed.

Every run writes `run_manifest.json` **last**, recording each file it produced
with its byte size plus the status of every step. Verify against it:

```bash
satellome --verify-run output_dir
```

* exit `0` — the directory matches its manifest and no step failed
* exit `1` — not a verifiably complete run: no/corrupt manifest, a failed step, a
  missing file, a leftover `*.partial`, or a file whose size no longer matches
  what the run wrote (the truncated-copy case)
* exit `2` — the argument is not a directory

Files already compressed by your own pipeline are still checked: for a missing
`X` with an `X.gz` next to it, the gzip ISIZE trailer (the uncompressed length
the compressor actually consumed) is compared to the recorded size, which catches
a `.gz` made from an incomplete read. Above 4 GiB that comparison is modulo
4 GiB and the report says so.

Two other guarantees back this up:

* **Atomic outputs** — files are written as `<path>.partial` and renamed into
  place, so a final name never refers to a half-written file. If you compress or
  copy an output directory concurrently, you either get the complete file or no
  file, never a truncated one.
* **Output-directory lock** — a second satellome run into the same `-o` is
  refused, naming the pid and host that holds it, instead of overwriting the
  first run's files mid-write. Override with `--ignore-lock` only if you are sure.

A run whose drawing step fails still writes a manifest — with `drawing: failed` —
and exits non-zero. The data files are complete and usable; the failed step is
recorded in the run's own artifact rather than only in an exit code.

## SAT File Format

The SAT format is a tab-delimited file with the following columns:

| Column | Description |
|--------|-------------|
| project | Project name |
| trf_id | Unique array ID |
| trf_head | Chromosome/scaffold name |
| trf_l_ind | Left coordinate (1-based) |
| trf_r_ind | Right coordinate |
| trf_period | Monomer period length |
| trf_n_copy | Number of copies |
| trf_pmatch | Percent match |
| trf_pvar | Percent variation |
| trf_entropy | Shannon entropy |
| trf_consensus | Consensus monomer sequence |
| trf_array | Full array sequence |
| trf_array_gc | Array GC content |
| trf_consensus_gc | Consensus GC content |
| trf_array_length | Array length in bp |
| trf_joined | Join status |
| trf_family | Repeat family |
| trf_ref_annotation | Reference annotation |

## Classification System

Satellome classifies tandem repeats into four categories:

| Category | Description | Criteria |
|----------|-------------|----------|
| **micro** | Microsatellites | Monomer 1-9 bp |
| **complex** | Complex repeats | Monomer >9 bp, entropy >1.82 |
| **pmicro** | Potential microsatellites | Intermediate characteristics |
| **tssr** | Tandem simple sequence repeats | Simple patterns |

## Example Results

Analysis of CHM13 v2.0 human genome (3.1 GB):

| Category | Arrays | % Genome |
|----------|--------|----------|
| Total | 614,616 | - |
| Complex | 20,373 | 5.27% |
| Microsatellites | 319,489 | 1.96% |
| TSSR | 296,475 | 0.47% |
| >1kb | 14,438 | 7.69% |
| >10kb | 1,223 | 6.67% |

## Utility Scripts

### Format Conversion
```bash
python scripts/trf_to_fasta.py -i repeats.sat -o repeats.fasta
python scripts/trf_to_gff3.py -i repeats.sat -o repeats.gff3
```

### Analysis Tools
```bash
python scripts/trf_get_large.py -i repeats.sat -m 1000 -o large_repeats.sat
python scripts/trf_get_micro_stat.py -i repeats.sat -o micro_stats.txt
python scripts/check_telomeres.py -i genome.fasta -t repeats.sat
```

## Testing

```bash
pytest tests/unit/ -v
```

## Contributing

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## Citation

If you use Satellome in your research, please cite:

```
Komissarov A. et al. (2026). Satellome: A comprehensive tool for satellite DNA
analysis in T2T genome assemblies. [Publication details]
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Support

- **Issues**: [GitHub Issues](https://github.com/aglabx/satellome/issues)
- **Documentation**: [Wiki](https://github.com/aglabx/satellome/wiki)
- **Email**: ad3002@gmail.com

## Acknowledgments

- [FasTAN](https://github.com/thegenemyers/FASTAN) by Gene Myers
- [Tandem Repeat Finder](https://github.com/Benson-Genomics-Lab/TRF) by Gary Benson
- [T2T Consortium](https://www.genome.gov/about-nhgri/telomere-to-telomere)
- [DNA Zoo](https://www.dnazoo.org/)
- [Vertebrate Genome Project](https://vertebrategenomesproject.org/)
