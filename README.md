# Satellome

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

## Usage

### Basic Command

```bash
satellome -i genome.fasta -o output_dir -p project_name -t 8
```

### Advanced Options

```bash
# With GFF3 annotations
satellome -i genome.fasta -o output_dir -p project_name -t 8 --gff annotations.gff3

# With RepeatMasker annotations
satellome -i genome.fasta -o output_dir -p project_name -t 8 --rm repeatmasker.out

# Force rerun all steps
satellome -i genome.fasta -o output_dir -p project_name -t 8 --force

# Custom TRF binary path
satellome -i genome.fasta -o output_dir -p project_name -t 8 --trf /path/to/trf

# Parallel processing of multiple genomes
python scripts/run_satellome_parallel.py -i genomes_list.txt -o results_dir -t 32
```

### Parameters

- `-i, --input`: Input FASTA file (required)
- `-o, --output`: Output directory (required)
- `-p, --project`: Project name (required)
- `-t, --threads`: Number of threads (default: 1)
- `--gff`: GFF3 annotation file (optional)
- `--rm`: RepeatMasker output file (optional)
- `--trf`: Path to TRF binary (default: "trf")
- `--force`: Force rerun all steps

## Output Structure

```
output_dir/
├── genome_name.trf                   # Main TRF output file
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
├── images/
│   ├── *.png                         # Karyotype and other visualizations
│   └── *.svg                         # Vector graphics versions
└── reports/
    ├── satellome_report.html         # Comprehensive HTML report
    └── annotation_report.txt         # Annotation intersection report (if GFF provided)
```

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
# Extract large tandem repeats
python scripts/trf_get_large.py -i repeats.trf -m 1000 -o large_repeats.trf

# Get microsatellite statistics
python scripts/trf_get_micro_stat.py -i repeats.trf -o micro_stats.txt

# Check telomeric repeats
python scripts/check_telomeres.py -i genome.fasta -t repeats.trf
```

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
# Download a DNA Zoo assembly (example: Aardvark)
wgethttps://dnazoo.s3.wasabisys.com/Acinonyx_jubatus/aciJub1_HiC.fasta.gz
gzip -d aciJub1_HiC.fasta.gz

# Run satellome on DNA Zoo assembly
satellome -i aciJub1_HiC.fasta \
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