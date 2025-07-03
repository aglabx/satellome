#!/bin/bash

# Run Satellome analysis on multiple genomes in parallel
# Usage: ./run_satellome_batch.sh <input_dir> <output_dir> [jobs] [threads_per_job]

if [ $# -lt 2 ]; then
    echo "Usage: $0 <input_dir> <output_dir> [jobs] [threads_per_job]"
    echo "Example: $0 /path/to/genomes /path/to/output 10 10"
    echo ""
    echo "Arguments:"
    echo "  input_dir       - Directory containing .fa.gz genome files"
    echo "  output_dir      - Directory for analysis results"
    echo "  jobs            - Number of parallel jobs (default: 10)"
    echo "  threads_per_job - Threads per Satellome job (default: 10)"
    exit 1
fi

INPUT_DIR="$1"
OUTPUT_DIR="$2"
JOBS="${3:-10}"
THREADS="${4:-10}"

# Calculate total CPU usage
TOTAL_CPUS=$((JOBS * THREADS))

echo "=========================================="
echo "Satellome Batch Analysis"
echo "=========================================="
echo "Input directory: $INPUT_DIR"
echo "Output directory: $OUTPUT_DIR"
echo "Parallel jobs: $JOBS"
echo "Threads per job: $THREADS"
echo "Total CPU cores: $TOTAL_CPUS"
echo "=========================================="

# Check if input directory exists
if [ ! -d "$INPUT_DIR" ]; then
    echo "Error: Input directory does not exist: $INPUT_DIR"
    exit 1
fi

# Count genome files
GENOME_COUNT=$(find "$INPUT_DIR" \( -name "*.fa.gz" -o -name "*.fasta.gz" -o -name "*.fna.gz" \) -type f | wc -l)
echo "Found $GENOME_COUNT genome files (.fa.gz, .fasta.gz, .fna.gz)"

if [ "$GENOME_COUNT" -eq 0 ]; then
    echo "Error: No FASTA files found in $INPUT_DIR"
    echo "Looking for: *.fa.gz, *.fasta.gz, *.fna.gz"
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Create log directory
LOG_DIR="$OUTPUT_DIR/logs"
mkdir -p "$LOG_DIR"

# Run the analysis
SCRIPT_DIR="$(dirname "$0")"
LOG_FILE="$LOG_DIR/satellome_batch_$(date +%Y%m%d_%H%M%S).log"

echo ""
echo "Starting analysis..."
echo "Log file: $LOG_FILE"
echo ""

# Run with nohup for long-running jobs
nohup python3 "$SCRIPT_DIR/run_satellome_parallel.py" \
    -i "$INPUT_DIR" \
    -o "$OUTPUT_DIR" \
    -j "$JOBS" \
    -t "$THREADS" \
    > "$LOG_FILE" 2>&1 &

PID=$!
echo "Process started with PID: $PID"
echo ""
echo "To monitor progress:"
echo "  tail -f $LOG_FILE"
echo ""
echo "To check if still running:"
echo "  ps -p $PID"