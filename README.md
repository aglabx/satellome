# Satellome

[![Tests](https://github.com/aglabx/satellome/workflows/Tests/badge.svg)](https://github.com/aglabx/satellome/actions?query=workflow%3ATests)
[![codecov](https://codecov.io/gh/aglabx/satellome/branch/main/graph/badge.svg)](https://codecov.io/gh/aglabx/satellome)
[![Python Version](https://img.shields.io/badge/python-3.9%20%7C%203.10%20%7C%203.11-blue)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![PyPI version](https://badge.fury.io/py/satellome.svg)](https://badge.fury.io/py/satellome)

A comprehensive bioinformatics tool for analyzing satellite DNA (tandem repeats) in telomere-to-telomere (T2T) genome assemblies.

## Overview

Satellome integrates Tandem Repeat Finder (TRF) to identify, classify, and visualize repetitive DNA sequences, with a particular focus on centromeric and telomeric regions. It provides a complete pipeline from raw genome sequences to detailed visualizations and reports of tandem repeat patterns.

The tool is designed to work with various genome assembly projects including:
- T2T (Telomere-to-Telomere) Consortium assemblies
- DNA Zoo chromosome-length assemblies
- VGP (Vertebrate Genome Project) assemblies
- NCBI RefSeq and GenBank assemblies

## Features

- **Tandem Repeat Detection**: Automated detection using TRF with optimized parameters
- **Smart Classification**: Categorizes repeats into microsatellites, complex repeats, and other types
- **Rich Visualizations**: Generates karyotype plots, 3D visualizations, and distance matrices
- **Annotation Integration**: Supports GFF3 and RepeatMasker annotations
- **Parallel Processing**: Efficient handling of multiple genomes
- **Smart Pipeline**: Automatically skips completed steps (override with `--force`)
- **Compressed File Support**: Direct processing of .gz compressed FASTA files
- **K-mer Based Filtering**: Optional k-mer profiling to focus on repeat-rich regions and skip repeat-poor areas

## Installation

### Prerequisites

- Python 3.9 or higher
- Conda (recommended) or pip
- TRF (Tandem Repeat Finder) binary

### Quick Setup

1. **Clone the repository**
```bash
git clone https://github.com/aglabx/satellome.git
cd satellome
```

2. **Create conda environment**
```bash
conda create -n satellome python=3.9
conda activate satellome
```

3. **Install dependencies**
```bash
pip install -r requirements.txt
```

4. **Install satellome**
```bash
pip install -e .  # Development mode
# or
pip install .     # Production mode
```

**Note**: During installation, Satellome will **automatically** attempt to install external tools (FasTAN, tanbed, modified TRF). This process:
- Compiles tools from source (requires: git, make, gcc/clang)
- Installs binaries to `<site-packages>/satellome/bin/` (or `~/.satellome/bin/` if no write permissions)
- Takes 2-5 minutes depending on your system
- **Can be skipped**: `SATELLOME_SKIP_AUTO_INSTALL=1 pip install satellome`
- If compilation fails, Satellome will still install successfully
- Failed tools can be installed later with `satellome --install-all`

5. **Download TRF binary**
```bash
# Linux
wget https://github.com/Benson-Genomics-Lab/TRF/releases/download/v4.09.1/trf409.linux64
chmod +x trf409.linux64
mv trf409.linux64 trf

# macOS
wget https://github.com/Benson-Genomics-Lab/TRF/releases/download/v4.09.1/trf409.macosx
chmod +x trf409.macosx
mv trf409.macosx trf
```

### TRF for Large Genomes (Chromosomes >1-2 GB)

**Important**: The standard TRF binary has limitations with very large chromosomes (>1-2 GB) and may crash during analysis. For large genome assemblies (e.g., some plant genomes, salamander genomes), use our modified TRF version.

#### Automatic Installation (Linux)

```bash
# Install modified TRF automatically (Linux only)
satellome --install-trf-large
```

Binary will be installed to `<site-packages>/satellome/bin/trf-large` (or `~/.satellome/bin/trf-large` as fallback).

**Note**: Automatic installation works best on Linux. macOS users may encounter compilation issues and should use manual installation or pre-compiled binaries.

#### Manual Installation

```bash
# Clone and build modified TRF
git clone https://github.com/aglabx/trf.git
cd trf
mkdir build && cd build
../configure
make

# Copy to system or Satellome directory
cp src/trf ~/.satellome/bin/trf-large
```

For pre-compiled binaries, visit: https://github.com/aglabx/trf/releases

**When to use the modified TRF:**
- Working with genomes containing chromosomes larger than 1-2 GB
- Experiencing crashes or "Segmentation fault" errors with standard TRF
- Processing large plant or amphibian genomes

The modified TRF includes memory optimizations and can handle chromosomes up to several gigabases in size. Specify the path using: `--trf ~/.satellome/bin/trf-large`

### FasTAN and tanbed (Optional)

Satellome supports **FasTAN** as an alternative tandem repeat finder. FasTAN and its companion tool **tanbed** can be automatically installed:

#### Automatic Installation

```bash
# Install FasTAN only
satellome --install-fastan

# Install tanbed only
satellome --install-tanbed

# Install both FasTAN and tanbed
satellome --install-all
```

**Note**: These tools are automatically installed during `pip install satellome`. Manual installation is only needed if automatic installation failed or was skipped.

Binaries will be installed to `<site-packages>/satellome/bin/` (or `~/.satellome/bin/` as fallback).

#### Requirements for Installation

The automatic installer requires:
- **git**: For cloning repositories
- **make**: For building
- **C compiler**: gcc, clang, or cc

**On macOS:**
```bash
xcode-select --install
```

**On Ubuntu/Debian:**
```bash
sudo apt-get install build-essential git
```

**On CentOS/RHEL:**
```bash
sudo yum groupinstall 'Development Tools'
sudo yum install git
```

#### Manual Installation

If you prefer manual installation or encounter issues:

**FasTAN:**
```bash
git clone https://github.com/thegenemyers/FASTAN.git
cd FASTAN
make
cp FasTAN ~/.satellome/bin/fastan
```

**tanbed:**
```bash
git clone https://github.com/richarddurbin/alntools.git
cd alntools
make
cp tanbed ~/.satellome/bin/tanbed
```

## Usage

### Basic Command

```bash
# Note: Output directory must be an absolute path
satellome -i genome.fasta -o /absolute/path/to/output_dir -p project_name -t 8
```

### Advanced Options

```bash
# With GFF3 annotations
satellome -i genome.fasta -o output_dir -p project_name -t 8 --gff annotations.gff3

# With RepeatMasker annotations
satellome -i genome.fasta -o output_dir -p project_name -t 8 --rm repeatmasker.out

# Force rerun all steps
satellome -i genome.fasta -o output_dir -p project_name -t 8 --force

# Smart recompute: only process chromosomes that failed TRF analysis
satellome -i genome.fasta -o output_dir -p project_name -t 8 --recompute-failed

# Custom TRF binary path (if not in PATH)
satellome -i genome.fasta -o /absolute/path/to/output_dir -p project_name -t 8 --trf /path/to/trf409.macosx

# Parallel processing of multiple genomes
python scripts/run_satellome_parallel.py -i genomes_list.txt -o results_dir -t 32

# With k-mer filtering to skip repeat-poor regions
satellome -i genome.fasta -o output_dir -p project_name -t 8 --use_kmer_filter

# Use pre-computed k-mer profile
varprofiler genome.fasta genome.varprofile.bed 17 100000 25000 20
satellome -i genome.fasta -o output_dir -p project_name -t 8 --kmer_bed genome.varprofile.bed

# Adjust k-mer threshold (default 90000)
satellome -i genome.fasta -o output_dir -p project_name -t 8 --use_kmer_filter --kmer_threshold 70000

# Continue with partial results if some TRF runs fail
satellome -i genome.fasta -o output_dir -p project_name -t 8 --continue-on-error

# Skip FasTAN analysis (run TRF only)
satellome -i genome.fasta -o output_dir -p project_name -t 8 --nofastan

# Skip TRF analysis (run FasTAN only)
satellome -i genome.fasta -o output_dir -p project_name -t 8 --notrf
```

### Parameters

- `-i, --input`: Input FASTA file (supports .fa, .fasta, .fa.gz, .fasta.gz)
- `-o, --output`: Output directory (required, must be an absolute path)
- `-p, --project`: Project name (required)
- `-t, --threads`: Number of threads (default: 1)
- `--gff`: GFF3 annotation file (optional)
- `--rm`: RepeatMasker output file (optional)
- `--trf`: Path to TRF binary (default: "trf")
- `--force`: Force rerun all steps
- `--recompute-failed`: Smart recompute - only process chromosomes/contigs that failed TRF analysis (missing from results)
- `--nofastan`: Skip FasTAN analysis (TRF runs by default)
- `--notrf`: Skip TRF analysis (FasTAN runs by default)
- `--use_kmer_filter`: Enable k-mer based filtering of repeat-poor regions
- `--kmer_threshold`: Threshold for unique k-mers (default: 90000)
- `--kmer_bed`: Pre-computed k-mer profile BED file from varprofiler
- `--continue-on-error`: Continue pipeline even if some TRF runs fail (results may be incomplete)

**Note:** By default, both TRF and FasTAN run on every analysis. Use `--nofastan` or `--notrf` to skip either tool. At least one tool must run.

## Output Structure

```
output_dir/
├── genome_name.trf                   # Main TRF output file
├── genome_name.gaps.bed              # Gaps annotation in BED format
├── genome_name.1kb.trf               # Repeats >1kb
├── genome_name.3kb.trf               # Repeats >3kb
├── genome_name.10kb.trf              # Repeats >10kb
├── genome_name.micro.trf             # Microsatellites (1-9 bp monomers)
├── genome_name.complex.trf           # Complex repeats (>9 bp monomers)
├── genome_name.pmicro.trf            # Potential microsatellites
├── genome_name.tssr.trf              # Tandem simple sequence repeats
├── genome_name.*.gff3                # GFF3 format files for each category
├── genome_name.*.fa                  # FASTA files with repeat sequences
├── distances.tsv.*                   # Distance matrices with various extensions
├── fastan/
│   ├── project_name.1aln             # FasTAN alignment output
│   └── project_name.bed              # FasTAN results in BED format
├── images/
│   ├── *.png                         # Karyotype and other visualizations
│   └── *.svg                         # Vector graphics versions
└── reports/
    ├── satellome_report.html         # Comprehensive HTML report
    └── annotation_report.txt         # Annotation intersection report (if GFF provided)
```

**Note:** The `fastan/` directory and gap annotation file are generated by default. Use `--nofastan` to skip FasTAN analysis.

## Classification System

Satellome classifies tandem repeats into four categories:

1. **micro**: Microsatellites (monomer length 1-9 bp)
2. **complex**: Complex repeats (monomer length >9 bp)
3. **pmicro**: Potential microsatellites
4. **tssr**: Tandem simple sequence repeats

## Utility Scripts

### Format Conversion
```bash
# Convert TRF to FASTA
python scripts/trf_to_fasta.py -i repeats.trf -o repeats.fasta

# Convert TRF to GFF3
python scripts/trf_to_gff3.py -i repeats.trf -o repeats.gff3

# Extract coordinates
python scripts/trf_to_coordinates.py -i repeats.trf -o coordinates.txt
```

### Analysis Tools
```bash
# Check TRF consistency - verify all large scaffolds have results
python scripts/check_trf_consistency.py -f genome.fasta -t output_dir/genome.trf
python scripts/check_trf_consistency.py -f genome.fasta -t output_dir/genome.trf -s 500000 -o report.txt

# Extract large tandem repeats
python scripts/trf_get_large.py -i repeats.trf -m 1000 -o large_repeats.trf

# Get microsatellite statistics
python scripts/trf_get_micro_stat.py -i repeats.trf -o micro_stats.txt

# Check telomeric repeats
python scripts/check_telomeres.py -i genome.fasta -t repeats.trf

# Check TRF results consistency
python scripts/check_trf_consistency.py -f genome.fna -t genome.trf

# Batch check TRF consistency for multiple genomes
python scripts/batch_check_trf_consistency.py reptiles mammals birds
```

### Quality Control Scripts

#### check_trf_consistency.py
Verifies that TRF analysis completed successfully for all contigs/scaffolds above a certain size threshold.

```bash
# Basic usage
python scripts/check_trf_consistency.py -f genome.fna -t genome.trf

# With custom minimum scaffold size (default: 1Mb)
python scripts/check_trf_consistency.py -f genome.fna -t genome.trf -s 500000

# With debug information for troubleshooting
python scripts/check_trf_consistency.py -f genome.fna -t genome.trf --debug

# Save detailed report
python scripts/check_trf_consistency.py -f genome.fna -t genome.trf -o report.txt
```

#### batch_check_trf_consistency.py
Batch process multiple genome assemblies to check TRF consistency.

```bash
# Check multiple directories
python scripts/batch_check_trf_consistency.py reptiles mammals birds

# Auto-skip failed assemblies
python scripts/batch_check_trf_consistency.py reptiles --auto-skip

# Show assemblies that need TRF analysis
python scripts/batch_check_trf_consistency.py reptiles --check-missing

# With progress tracking and debug info
python scripts/batch_check_trf_consistency.py reptiles --debug --verbose

# Save summary report
python scripts/batch_check_trf_consistency.py reptiles -o consistency_report.txt
```

**Interactive mode options:**
- `[s]` Skip - continue to next assembly
- `[d]` Delete - remove TRF directory and re-run TRF
- `[v]` View - show TRF directory contents
- `[q]` Quit - exit the script

### Smart Recompute Mode

If TRF analysis fails for some chromosomes (e.g., due to memory issues or signal errors), you can use the `--recompute-failed` flag to reprocess only the failed chromosomes without redoing the entire analysis.

**How it works:**
1. Checks which chromosomes/contigs are missing from existing TRF results
2. Extracts only those chromosomes to a temporary FASTA file
3. Runs TRF only on the missing chromosomes
4. Merges results back into the existing TRF file
5. Continues with the rest of the pipeline

**Usage example:**
```bash
# First, check which chromosomes failed
python scripts/check_trf_consistency.py -f genome.fna -t output_dir/project.trf

# Then recompute only the failed ones
satellome -i genome.fasta -o output_dir -p project_name -t 8 --recompute-failed
```

**When to use:**
- TRF failed for specific chromosomes (visible in error messages like "TRF failed for 94.fa")
- `check_trf_consistency.py` reports missing chromosomes
- You want to save time by not reprocessing successful chromosomes

**Benefits:**
- Much faster than `--force` (only processes failed chromosomes)
- Preserves successful results
- Creates automatic backup before merging (`.before_recompute` suffix)
- More informative error messages with actual TRF output

## Example Workflow

### 1. Download Test Dataset
```bash
# Download S. cerevisiae genome
curl -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_000146045.2/download?include_annotation_type=GENOME_FASTA,GENOME_GFF&filename=GCF_000146045.2.zip" -H "Accept: application/zip"
unzip GCF_000146045.2.zip
```

### 2. Run Analysis
```bash
# Run satellome pipeline
satellome -i ncbi_dataset/data/GCF_000146045.2/GCF_000146045.2_R64_genomic.fna \
          -o results \
          -p scerevisiae \
          -t 8 \
          --gff ncbi_dataset/data/GCF_000146045.2/genomic.gff

# View results
open results/scerevisiae_report.html
```

### 3. Analyzing DNA Zoo Assemblies
```bash
# Download a DNA Zoo assembly (example: Cheetah)
wget https://dnazoo.s3.wasabisys.com/Acinonyx_jubatus/aciJub1_HiC.fasta.gz

# Run satellome directly on compressed file (no need to decompress!)
satellome -i aciJub1_HiC.fasta.gz \
          -o dnazoo_results \
          -p cheetah \
          -t 8
```

## Configuration

The pipeline uses `settings.yaml` for tool parameters. Key settings include:

- TRF parameters (match/mismatch scores, indel penalties)
- Minimum/maximum repeat lengths
- Classification thresholds
- Visualization parameters

## Testing

Run the test suite:
```bash
python tests/test_overlapping.py
python test_standalone.py
python test_chromosome_sorting.py
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

- [Tandem Repeat Finder](https://github.com/Benson-Genomics-Lab/TRF) by Gary Benson
- [T2T Consortium](https://www.genome.gov/about-nhgri/telomere-to-telomere) for inspiring this work
- [DNA Zoo](https://www.dnazoo.org/) for providing chromosome-length assemblies
- [Vertebrate Genome Project](https://vertebrategenomesproject.org/) for high-quality reference genomes