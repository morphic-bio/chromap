#!/bin/bash
# Memory comparison script: Compare Max RSS between current and baseline binaries
# Usage: compare_memory.sh [OUTPUT_DIR]

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"
source "$SCRIPT_DIR/lib/utils.sh"

OUTPUT_DIR="${1:-$OUT_ROOT/memory}"

if [ ! -f "$BASELINE_BIN" ]; then
    log_warning "Baseline binary not found: $BASELINE_BIN"
    log_warning "Skipping memory comparison"
    exit 0
fi

mkdir -p "$OUTPUT_DIR"

# Check for time command
if ! command -v /usr/bin/time &> /dev/null; then
    log_error "/usr/bin/time not found (required for memory comparison)"
    exit 1
fi

# Resolve read files (use small dataset for memory test)
resolve_reads "$DATA_SMALL2" READ1 READ2 || {
    log_error "Failed to resolve read files"
    exit 1
}

log_info "Running memory comparison test"
log_info "Dataset: $(basename $DATA_SMALL2)"
log_info "Reads: $(basename $READ1) / $(basename $READ2)"

# Build common args
COMMON_ARGS="$CHROMAP_ARGS_BASE -x $REF_INDEX -r $REF_FA -1 $READ1 -2 $READ2 -t $NUM_THREADS"

# Run current binary
log_info "Running current binary: $CHROMAP_BIN"
(/usr/bin/time -v "$CHROMAP_BIN" $COMMON_ARGS --SAM -o "$OUTPUT_DIR/new.sam" 2>&1) > "$OUTPUT_DIR/new.time" || {
    log_error "Current binary failed"
    exit 1
}

# Run baseline binary
log_info "Running baseline binary: $BASELINE_BIN"
(/usr/bin/time -v "$BASELINE_BIN" $COMMON_ARGS --SAM -o "$OUTPUT_DIR/old.sam" 2>&1) > "$OUTPUT_DIR/old.time" || {
    log_error "Baseline binary failed"
    exit 1
}

# Parse Max RSS from time output
NEW_RSS=$(grep "Maximum resident set size" "$OUTPUT_DIR/new.time" | awk '{print $6}' || echo "0")
OLD_RSS=$(grep "Maximum resident set size" "$OUTPUT_DIR/old.time" | awk '{print $6}' || echo "0")

if [ "$NEW_RSS" = "0" ] || [ "$OLD_RSS" = "0" ]; then
    log_error "Failed to parse Max RSS from time output"
    log_info "New binary output:"
    grep -A 2 "Maximum resident set size" "$OUTPUT_DIR/new.time" || true
    log_info "Baseline binary output:"
    grep -A 2 "Maximum resident set size" "$OUTPUT_DIR/old.time" || true
    exit 1
fi

# Calculate difference
DIFF=$((NEW_RSS - OLD_RSS))
PERCENT=$(awk "BEGIN {printf \"%.1f\", ($NEW_RSS / $OLD_RSS) * 100}")

# Convert to MB for readability
NEW_MB=$((NEW_RSS / 1024))
OLD_MB=$((OLD_RSS / 1024))
DIFF_MB=$((DIFF / 1024))

# Report results
log_info "Memory comparison results:"
echo "  Current binary:  ${NEW_MB} MB (${NEW_RSS} KB)"
echo "  Baseline binary: ${OLD_MB} MB (${OLD_RSS} KB)"
echo "  Difference:      ${DIFF_MB} MB (${DIFF} KB)"
echo "  Percentage:     ${PERCENT}%"

# Check threshold (110% = 1.1)
THRESHOLD=$(awk "BEGIN {printf \"%.0f\", $OLD_RSS * 1.1}")
if [ "$NEW_RSS" -gt "$THRESHOLD" ]; then
    log_warning "Current binary uses more than 110% of baseline memory"
    exit 1
fi

log_success "Memory comparison passed (within 110% threshold)"
exit 0

