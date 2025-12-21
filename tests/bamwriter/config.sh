#!/bin/bash
# Configuration file for BAM/CRAM writer test suite
# Source this file to set up all test paths and parameters

# Binaries
CHROMAP_BIN="${CHROMAP_BIN:-/mnt/pikachu/chromap/chromap}"
BASELINE_BIN="${BASELINE_BIN:-/usr/local/bin/chromap}"

# Reference files (defaults from runChromap.sh)
REF_FA="${REF_FASTA:-/storage/scRNAseq_output/indices-110-44/chromap/genome.fa}"
REF_INDEX="${INDEX_FILE:-/storage/scRNAseq_output/indices-110-44/chromap/index}"
CHROM_SIZES="${CHROM_SIZES:-/storage/scRNAseq_output/indices-110-44/chromap/chrNameLength.txt}"

# Test datasets
DATA_SMALL1=/storage/ATAC-Seq-10000
DATA_SMALL2=/storage/ATAC-short-test
DATA_LARGE=/mnt/pikachu/NW-5-21/ATAC-Seq

# Common chromap flags (morphic pipeline defaults)
CHROMAP_ARGS_BASE="--trim-adapters --min-frag-length 30 -l 2000 --remove-pcr-duplicates"

# Output root directory
OUT_ROOT="${OUT_ROOT:-/mnt/pikachu/chromap/tests/out/bamwriter}"

# Test parameters
NUM_THREADS="${NUM_THREADS:-1}"  # Use 1 thread for deterministic output
HTS_THREADS="${HTS_THREADS:-4}"  # Compression threads for BAM/CRAM

# Y-filter stream flags
Y_STREAMS="--emit-noY-bam --emit-Y-bam"

# Export all variables
export CHROMAP_BIN BASELINE_BIN REF_FA REF_INDEX CHROM_SIZES
export DATA_SMALL1 DATA_SMALL2 DATA_LARGE
export CHROMAP_ARGS_BASE OUT_ROOT NUM_THREADS HTS_THREADS Y_STREAMS

