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
Komissarov A. et al. (2024). Satellome: A comprehensive tool for satellite DNA
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
