#!/bin/bash
# Parity check script: Compare SAM vs BAM/CRAM outputs
# Usage: check_parity.sh SAM_FILE BAM_OR_CRAM_FILE [REF_FASTA]

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"
source "$SCRIPT_DIR/lib/utils.sh"

SAM_FILE="${1:-}"
BAM_FILE="${2:-}"
REF_FASTA="${3:-$REF_FA}"
SORT_OUTPUTS="${4:-false}"  # If true, sort outputs before comparison (for lowmem mode)

if [ -z "$SAM_FILE" ] || [ -z "$BAM_FILE" ]; then
    log_error "Usage: $0 SAM_FILE BAM_OR_CRAM_FILE [REF_FASTA]"
    exit 1
fi

# Check files exist
if [ ! -f "$SAM_FILE" ]; then
    log_error "SAM file not found: $SAM_FILE"
    exit 1
fi

if [ ! -f "$BAM_FILE" ]; then
    log_error "BAM/CRAM file not found: $BAM_FILE"
    exit 1
fi

# Determine file type
BAM_EXT="${BAM_FILE##*.}"
TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

SAM_BASE=$(basename "$SAM_FILE" .sam)
BAM_BASE=$(basename "$BAM_FILE" .${BAM_EXT})

log_info "Comparing $SAM_BASE vs $BAM_BASE"

# Step 1: Convert to canonical SAM
log_info "Converting to canonical SAM..."

# SAM input (re-emit through samtools for canonical tag ordering)
samtools view -h -O SAM "$SAM_FILE" > "$TMPDIR/sam.canon.sam" 2>/dev/null || {
    log_error "Failed to convert SAM file"
    exit 1
}

# BAM/CRAM input
if [ "$BAM_EXT" = "cram" ]; then
    samtools view -h -T "$REF_FASTA" -O SAM "$BAM_FILE" > "$TMPDIR/bam.canon.sam" 2>/dev/null || {
        log_error "Failed to convert CRAM file (check reference: $REF_FASTA)"
        exit 1
    }
else
    samtools view -h -O SAM "$BAM_FILE" > "$TMPDIR/bam.canon.sam" 2>/dev/null || {
        log_error "Failed to convert BAM file"
        exit 1
    }
fi

# Step 2: Split header and body
if [ "$SORT_OUTPUTS" = "true" ]; then
    log_info "Splitting headers and bodies (will sort outputs for comparison)..."
else
    log_info "Splitting headers and bodies..."
fi

grep '^@' "$TMPDIR/sam.canon.sam" | LC_ALL=C sort > "$TMPDIR/sam.header.txt"
# Sort body by entire line for stable comparison (as per runbook)
# If SORT_OUTPUTS=true, this ensures both files are sorted regardless of native sort order
grep -v '^@' "$TMPDIR/sam.canon.sam" | LC_ALL=C sort > "$TMPDIR/sam.body.txt"

grep '^@' "$TMPDIR/bam.canon.sam" | LC_ALL=C sort > "$TMPDIR/bam.header.txt"
# Sort body by entire line for stable comparison (as per runbook)
grep -v '^@' "$TMPDIR/bam.canon.sam" | LC_ALL=C sort > "$TMPDIR/bam.body.txt"

# Step 3: Compare headers (only @SQ lines, as per runbook)
# BAM/CRAM includes @HD with SO:coordinate which SAM may not have
# CRAM also includes M5 and UR fields in @SQ which SAM doesn't
# Normalize to only compare SN and LN fields
log_info "Comparing @SQ headers (SN and LN only)..."
grep '^@SQ' "$TMPDIR/sam.header.txt" | awk -F'\t' '{for(i=2;i<=NF;i++){if($i~/^SN:/){sn=$i} if($i~/^LN:/){ln=$i}} print sn"\t"ln}' | LC_ALL=C sort > "$TMPDIR/sam.header.sq.txt" || true
grep '^@SQ' "$TMPDIR/bam.header.txt" | awk -F'\t' '{for(i=2;i<=NF;i++){if($i~/^SN:/){sn=$i} if($i~/^LN:/){ln=$i}} print sn"\t"ln}' | LC_ALL=C sort > "$TMPDIR/bam.header.sq.txt" || true

HEADER_DIFF=$(diff -u "$TMPDIR/sam.header.sq.txt" "$TMPDIR/bam.header.sq.txt" || true)
if [ -n "$HEADER_DIFF" ]; then
    log_error "@SQ header mismatch (SN/LN fields)"
    echo "$HEADER_DIFF"
    exit 1
fi

# Step 4: Compare bodies
# BAM must match SAM exactly (including QUAL and all tags)
# CRAM may encode quality strings differently, so exclude QUAL (field 11) but include all tags
if [ "$BAM_EXT" = "bam" ]; then
    log_info "Comparing alignment bodies (BAM: all fields including QUAL and tags)..."
    # BAM: Compare full lines (all fields)
    LC_ALL=C sort "$TMPDIR/sam.body.txt" > "$TMPDIR/sam.full.txt"
    LC_ALL=C sort "$TMPDIR/bam.body.txt" > "$TMPDIR/bam.full.txt"
    
    BODY_DIFF=$(diff -u "$TMPDIR/sam.full.txt" "$TMPDIR/bam.full.txt" || true)
    if [ -n "$BODY_DIFF" ]; then
        log_error "Body mismatch (BAM should match SAM exactly, including QUAL and tags)"
        echo "$BODY_DIFF" | head -50
        exit 1
    fi
else
    # CRAM: Compare all fields except QUAL (field 11)
    log_info "Comparing alignment bodies (CRAM: all fields except QUAL, tags included)..."
    # Extract fields 1-10, then fields 12+ (skip field 11 which is QUAL)
    # Sort tags within each record for stable comparison (CRAM may reorder tags)
    awk 'BEGIN{OFS="\t"} {
        # Fields 1-10 (core fields)
        core = ""
        for(i=1;i<=10;i++) {
            if(i==1) core = $i
            else core = core "\t" $i
        }
        # Collect fields 12+ (tags) and sort them
        ntags = 0
        for(i=12;i<=NF;i++) {
            tags[ntags++] = $i
        }
        # Sort tags alphabetically for stable comparison
        asort(tags)
        tagstr = ""
        for(i=1;i<=ntags;i++) {
            if(i==1) tagstr = tags[i]
            else tagstr = tagstr "\t" tags[i]
        }
        print core "\t" tagstr
    }' "$TMPDIR/sam.body.txt" | LC_ALL=C sort > "$TMPDIR/sam.noqual.txt"
    
    awk 'BEGIN{OFS="\t"} {
        # Fields 1-10 (core fields)
        core = ""
        for(i=1;i<=10;i++) {
            if(i==1) core = $i
            else core = core "\t" $i
        }
        # Collect fields 12+ (tags) and sort them
        ntags = 0
        for(i=12;i<=NF;i++) {
            tags[ntags++] = $i
        }
        # Sort tags alphabetically for stable comparison
        asort(tags)
        tagstr = ""
        for(i=1;i<=ntags;i++) {
            if(i==1) tagstr = tags[i]
            else tagstr = tagstr "\t" tags[i]
        }
        print core "\t" tagstr
    }' "$TMPDIR/bam.body.txt" | LC_ALL=C sort > "$TMPDIR/bam.noqual.txt"
    
    BODY_DIFF=$(diff -u "$TMPDIR/sam.noqual.txt" "$TMPDIR/bam.noqual.txt" || true)
    if [ -n "$BODY_DIFF" ]; then
        log_error "Body mismatch (CRAM should match SAM except QUAL; tags must match)"
        echo "$BODY_DIFF" | head -50
        exit 1
    fi
fi

# Step 5: Verify read counts
log_info "Verifying read counts..."
SAM_COUNT=$(count_reads "$SAM_FILE")
BAM_COUNT=$(count_reads "$BAM_FILE" "$REF_FASTA")

if [ "$SAM_COUNT" != "$BAM_COUNT" ]; then
    log_error "Read count mismatch: SAM=$SAM_COUNT, ${BAM_EXT^^}=$BAM_COUNT"
    exit 1
fi

log_success "Parity check passed: $SAM_COUNT reads match"
exit 0

