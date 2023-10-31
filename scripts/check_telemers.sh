#!/bin/bash

# Check if the argument is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <BASE_DIR>"
    exit 1
fi

# Base directory containing the fasta files
BASE_DIR="$1"

# Loop through each .fna file in the directory tree
find "$BASE_DIR" -type f \( -name "*.fna" -o -name "*.fasta" \) ! -path "*/trf/*" | while read -r fasta_file; do


    # Get the base name without the extension
    base_name="${fasta_file%.*}"
    echo "Processing $fasta_file"
    # Run the check_telomeres.py script
    python ~/Dropbox/workspace/PyBioSnippets/satelome/scripts/check_telomeres.py -i "$fasta_file" -o "${base_name}.telomeres"
done
