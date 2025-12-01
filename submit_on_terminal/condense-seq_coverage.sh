#!/bin/bash
# ----------------------------------------------------------
# Script Name: run_coverage_by_chr.sh
# Description:
#   Loop through each chromosome and compute genome coverage
#   using coverage_ver3.py for a given experiment.
#
# Usage:
#   ./run_coverage_by_chr.sh <input_bam_dir> <experiment_prefix> <output_dir> <reference_fasta>
#
# Example:
#   ./run_coverage_by_chr.sh alignment Exp3 analysis ../genomes/mm10/Bowtie2Index/genome.fa
#
# Notes:
#   - All BAM files matching <experiment_prefix>*.bam inside <input_bam_dir>
#     will be processed.
#   - The script will run one job per chromosome.
# ----------------------------------------------------------

# This is a secure step but it might cause error when running in virtual environment. 
#   So, be careful when you are using it.
# set -euo pipefail

# --------------------------
# Parse arguments
# --------------------------
inDir=${1:-}
expPrefix=${2:-}
outDir=${3:-}
refFile=${4:-}

# Print help message if arguments are missing or -h/--help provided
if [[ $# -lt 4 || "$1" == "-h" || "$1" == "--help" ]]; then
    echo "Usage:"
    echo "  $0 <input_bam_dir> <experiment_prefix> <output_dir> <reference_fasta>"
    echo
    echo "Example:"
    echo "  $0 alignment Exp3 analysis ../genomes/mm10/Bowtie2Index/genome.fa"
    echo
    echo "Arguments:"
    echo "  <input_bam_dir>      Directory containing BAM files (e.g. alignment)"
    echo "  <experiment_prefix>  Common prefix for BAM files (e.g. Exp3)"
    echo "  <output_dir>         Directory to save coverage output files"
    echo "  <reference_fasta>    Reference genome FASTA file (used to detect chromosomes)"
    exit 0
fi

# --------------------------
# Validate input
# --------------------------
if [[ ! -d "$inDir" ]]; then
    echo "Error: Input directory '$inDir' does not exist." >&2
    exit 1
fi

if [[ ! -f "$refFile" ]]; then
    echo "Error: Reference FASTA file '$refFile' not found." >&2
    exit 1
fi

mkdir -p "$outDir"

# --------------------------
# Find BAM files
# --------------------------
files=($(find "$inDir" -type f -name "${expPrefix}*.bam"))
if [[ ${#files[@]} -eq 0 ]]; then
    echo "No BAM files found in $inDir with prefix '${expPrefix}'." >&2
    exit 1
fi

echo "Found ${#files[@]} BAM files for experiment prefix '${expPrefix}'."

# --------------------------
# Detect chromosome names from reference
# --------------------------
chr_list=($(grep '^>' "$refFile" | sed 's/>//' | awk '{print $1}'))
echo "Detected ${#chr_list[@]} chromosomes from $refFile."

# --------------------------
# Run coverage script per chromosome
# --------------------------
echo "[$(date)] Step I in condense-seq pipeline: calculating the coverage"
localScriptFolder=/Volumes/BackupHaLab/condense-seq
for chr_name in "${chr_list[@]}"; do
    echo "Processing chromosome: ${chr_name}"

    python3 ${localScriptFolder}/prepro_scripts/coverage_ver3.py \
        -f "${files[@]}" \
        -x "$refFile" \
        --chr "$chr_name" \
        --skip \
        -o "${outDir}/${expPrefix}_${chr_name}"
done

echo "[$(date)] ✅ All chromosomes processed successfully."
