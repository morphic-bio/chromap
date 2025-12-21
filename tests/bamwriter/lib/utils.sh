#!/bin/bash
# Shared utility functions for BAM/CRAM test suite

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Logging functions
log_info() {
    echo -e "${BLUE}[INFO]${NC} $*"
}

log_success() {
    echo -e "${GREEN}[PASS]${NC} $*"
}

log_warning() {
    echo -e "${YELLOW}[WARN]${NC} $*"
}

log_error() {
    echo -e "${RED}[FAIL]${NC} $*" >&2
}

# Check if command exists
check_command() {
    if ! command -v "$1" &> /dev/null; then
        log_error "$1 not found in PATH"
        return 1
    fi
    return 0
}

# Check prerequisites
check_prereqs() {
    local errors=0
    
    if [ ! -f "$CHROMAP_BIN" ]; then
        log_error "CHROMAP_BIN not found: $CHROMAP_BIN"
        errors=$((errors + 1))
    fi
    
    if [ ! -f "$BASELINE_BIN" ]; then
        log_warning "BASELINE_BIN not found: $BASELINE_BIN (memory comparison will be skipped)"
    fi
    
    if ! check_command samtools; then
        errors=$((errors + 1))
    fi
    
    if [ ! -f "$REF_FA" ]; then
        log_error "Reference FASTA not found: $REF_FA"
        errors=$((errors + 1))
    fi
    
    if [ ! -f "$REF_FA.fai" ]; then
        log_warning "Reference FASTA index not found: $REF_FA.fai (creating...)"
        samtools faidx "$REF_FA" || {
            log_error "Failed to create FASTA index"
            errors=$((errors + 1))
        }
    fi
    
    if [ ! -f "$REF_INDEX" ]; then
        log_error "Chromap index not found: $REF_INDEX"
        errors=$((errors + 1))
    fi
    
    return $errors
}

# Extract sample_id from R1 filename (same logic as runChromap.sh derive_sample_id)
# Handles: ..._S<digits>_R1_001.fastq.gz, ..._R1_001.fastq.gz, ..._R1.fastq.gz
derive_sample_id() {
    local r1base="$1"
    local sid="${r1base}"
    sid="${sid%.fastq.gz}"
    sid=$(echo "${sid}" | sed -E 's/_S[0-9]+_R1_001$|_R1_001$|_R1$//')
    echo "${sid}"
}

# Resolve read files for a dataset directory using runChromap.sh pairing logic
# Usage: resolve_reads DATASET_DIR READ1_VAR READ2_VAR
resolve_reads() {
    local dataset_dir="$1"
    local read1_var="$2"
    local read2_var="$3"
    
    # Find first R1 file matching runChromap.sh pattern: .*_R1(_001)?\.fastq\.gz$
    local read1=$(find "$dataset_dir" -maxdepth 1 -type f -regextype posix-extended \
        -regex '.*_R1(_001)?\.fastq\.gz$' | sort | head -1)
    
    if [ -z "$read1" ]; then
        log_error "Could not find R1 file in $dataset_dir"
        return 1
    fi
    
    local base=$(basename "$read1")
    local sample_id=$(derive_sample_id "$base")
    
    # Try direct replacement (runChromap.sh logic)
    local read2="${read1/_R1_/_R2_}"
    if [ "$read2" = "$read1" ]; then
        read2="${read1/_R1./_R2.}"
    fi
    
    # If still missing, try sample_id-based patterns
    if [ ! -f "$read2" ]; then
        local candidates=(
            "${dataset_dir}/${sample_id}_S"*"_R2_001.fastq.gz"
            "${dataset_dir}/${sample_id}_R2_001.fastq.gz"
            "${dataset_dir}/${sample_id}_R2.fastq.gz"
        )
        local found=""
        for pattern in "${candidates[@]}"; do
            for f in $pattern; do
                if [ -f "$f" ]; then
                    found="$f"
                    break
                fi
            done
            [ -n "$found" ] && break
        done
        if [ -n "$found" ]; then
            read2="$found"
        fi
    fi
    
    if [ ! -f "$read2" ]; then
        log_error "Could not find matching R2 file for $read1"
        log_error "Sample ID: $sample_id"
        return 1
    fi
    
    eval "$read1_var=\"$read1\""
    eval "$read2_var=\"$read2\""
    
    log_info "Paired reads: $(basename $read1) / $(basename $read2)"
    
    return 0
}

# Create output directory structure
setup_output_dir() {
    local dataset_name="$1"
    local output_dir="$OUT_ROOT/$dataset_name"
    
    mkdir -p "$output_dir"/{sam,bam,cram,sam_lowmem,bam_lowmem,cram_lowmem}
    mkdir -p "$OUT_ROOT/memory"
    mkdir -p "$OUT_ROOT/reports"
    
    echo "$output_dir"
}

# Get file size in human-readable format
file_size() {
    local file="$1"
    if [ -f "$file" ]; then
        du -h "$file" | cut -f1
    else
        echo "N/A"
    fi
}

# Count reads in SAM/BAM/CRAM file
count_reads() {
    local file="$1"
    local ref="${2:-}"
    
    if [ -z "$file" ] || [ ! -f "$file" ]; then
        echo "0"
        return
    fi
    
    local ext="${file##*.}"
    local cmd="samtools view -c"
    
    if [ "$ext" = "cram" ] && [ -n "$ref" ]; then
        cmd="$cmd -T $ref"
    fi
    
    $cmd "$file" 2>/dev/null || echo "0"
}

# Check if file is empty
is_empty() {
    local file="$1"
    [ ! -f "$file" ] || [ ! -s "$file" ]
}

