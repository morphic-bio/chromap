#!/bin/bash
# Index validation script: Check BAM/CRAM index files
# Usage: check_index.sh BAM_OR_CRAM_FILE [REF_FASTA]

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"
source "$SCRIPT_DIR/lib/utils.sh"

BAM_FILE="${1:-}"
REF_FASTA="${2:-$REF_FA}"

if [ -z "$BAM_FILE" ]; then
    log_error "Usage: $0 BAM_OR_CRAM_FILE [REF_FASTA]"
    exit 1
fi

if [ ! -f "$BAM_FILE" ]; then
    log_error "File not found: $BAM_FILE"
    exit 1
fi

# Determine file type and expected index extension
BAM_EXT="${BAM_FILE##*.}"
if [ "$BAM_EXT" = "bam" ]; then
    INDEX_EXT="csi"
    INDEX_FILE="${BAM_FILE}.csi"
elif [ "$BAM_EXT" = "cram" ]; then
    INDEX_EXT="crai"
    INDEX_FILE="${BAM_FILE}.crai"
else
    log_error "Unknown file type: $BAM_EXT (expected bam or cram)"
    exit 1
fi

log_info "Checking index for $(basename $BAM_FILE)"

# Check index file exists
if [ ! -f "$INDEX_FILE" ]; then
    log_error "Index file not found: $INDEX_FILE"
    exit 1
fi

# Check index file is non-empty
if [ ! -s "$INDEX_FILE" ]; then
    log_error "Index file is empty: $INDEX_FILE"
    exit 1
fi

log_info "Index file exists: $(basename $INDEX_FILE) ($(file_size $INDEX_FILE))"

# Run samtools quickcheck (index is auto-detected, don't pass explicitly)
# For CRAM, set REF_PATH environment variable instead of -r flag
log_info "Running samtools quickcheck..."
if [ "$BAM_EXT" = "cram" ]; then
    REF_PATH="$REF_FASTA" samtools quickcheck -v "$BAM_FILE" 2>&1 || {
        log_error "samtools quickcheck failed for CRAM"
        exit 1
    }
else
    samtools quickcheck -v "$BAM_FILE" 2>&1 || {
        log_error "samtools quickcheck failed for BAM"
        exit 1
    }
fi

# Verify index is usable with idxstats
log_info "Verifying index usability with samtools idxstats..."
if [ "$BAM_EXT" = "cram" ]; then
    REF_PATH="$REF_FASTA" samtools idxstats "$BAM_FILE" > /dev/null 2>&1 || {
        log_error "samtools idxstats failed for CRAM"
        exit 1
    }
else
    samtools idxstats "$BAM_FILE" > /dev/null 2>&1 || {
        log_error "samtools idxstats failed for BAM"
        exit 1
    }
fi

log_success "Index validation passed: $INDEX_FILE"
exit 0

