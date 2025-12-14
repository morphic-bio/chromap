#!/bin/bash
# Regression test for Y/NoY split feature
# Compares current chromap against baseline to ensure no unintended changes

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

# Get script directory and repo root
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

# Configuration
BASELINE_COMMIT="${BASELINE_COMMIT:-5bd17e1}"
TEST_DIR="${TEST_DIR:-$(mktemp -d)}"
DATA_DIR="$TEST_DIR/data"
WORKTREE_DIR="$TEST_DIR/baseline_worktree"

# Binary paths
CURRENT_BIN="${CURRENT_BIN:-$REPO_ROOT/chromap}"
BASELINE_BIN="${BASELINE_BIN:-}"

echo "=============================================="
echo "Y/NoY Split Regression Test"
echo "=============================================="
echo "Test directory: $TEST_DIR"
echo "Baseline commit: $BASELINE_COMMIT"
echo ""

# Cleanup function
cleanup() {
    if [ -n "$WORKTREE_DIR" ] && [ -d "$WORKTREE_DIR" ]; then
        echo "Cleaning up worktree..."
        cd "$REPO_ROOT"
        git worktree remove --force "$WORKTREE_DIR" 2>/dev/null || true
    fi
    if [ "${KEEP_TEST_DIR:-0}" != "1" ] && [ -n "$TEST_DIR" ] && [ -d "$TEST_DIR" ]; then
        echo "Cleaning up test directory: $TEST_DIR"
        rm -rf "$TEST_DIR"
    fi
}

# Only cleanup on exit if not keeping test dir
if [ "${KEEP_TEST_DIR:-0}" != "1" ]; then
    trap cleanup EXIT
fi

# Check for samtools
if ! command -v samtools &> /dev/null; then
    echo -e "${RED}ERROR: samtools is required but not found${NC}"
    exit 1
fi

# Check for python3
if ! command -v python3 &> /dev/null; then
    echo -e "${RED}ERROR: python3 is required but not found${NC}"
    exit 1
fi

# Check current binary exists
if [ ! -f "$CURRENT_BIN" ]; then
    echo -e "${YELLOW}Current binary not found at $CURRENT_BIN, building...${NC}"
    cd "$REPO_ROOT" && make
fi

if [ ! -f "$CURRENT_BIN" ]; then
    echo -e "${RED}ERROR: Failed to build current chromap binary${NC}"
    exit 1
fi
echo "Current binary: $CURRENT_BIN"

# Build or find baseline binary
if [ -z "$BASELINE_BIN" ] || [ ! -f "$BASELINE_BIN" ]; then
    echo ""
    echo "Building baseline binary from commit $BASELINE_COMMIT..."
    
    # Create worktree for baseline
    cd "$REPO_ROOT"
    git worktree add "$WORKTREE_DIR" "$BASELINE_COMMIT" 2>/dev/null || {
        echo -e "${RED}ERROR: Failed to create worktree for baseline commit${NC}"
        exit 1
    }
    
    # Build baseline
    cd "$WORKTREE_DIR"
    make clean && make || {
        echo -e "${RED}ERROR: Failed to build baseline chromap${NC}"
        exit 1
    }
    
    BASELINE_BIN="$WORKTREE_DIR/chromap"
fi

if [ ! -f "$BASELINE_BIN" ]; then
    echo -e "${RED}ERROR: Baseline binary not found at $BASELINE_BIN${NC}"
    exit 1
fi
echo "Baseline binary: $BASELINE_BIN"

# Generate test data
echo ""
echo "Generating test data..."
mkdir -p "$DATA_DIR"
python3 "$SCRIPT_DIR/data/generate_test_data.py" -o "$DATA_DIR"

# Build index (use current binary - index format should be compatible)
echo ""
echo "Building index..."
"$CURRENT_BIN" --build-index \
    -r "$DATA_DIR/test_ref.fa" \
    -o "$DATA_DIR/test_ref.idx" \
    > /dev/null 2>&1

# ============================================
# BASELINE RUNS
# ============================================
echo ""
echo "=============================================="
echo "Running baseline mappings..."
echo "=============================================="

# Single-end baseline
echo "  Single-end..."
"$BASELINE_BIN" --SAM \
    -r "$DATA_DIR/test_ref.fa" \
    -x "$DATA_DIR/test_ref.idx" \
    -1 "$DATA_DIR/test_se.fq" \
    -o "$TEST_DIR/baseline_se.sam" \
    --min-num-seeds 1 --error-threshold 10 \
    > /dev/null 2>&1

# Paired-end baseline
echo "  Paired-end..."
"$BASELINE_BIN" --SAM \
    -r "$DATA_DIR/test_ref.fa" \
    -x "$DATA_DIR/test_ref.idx" \
    -1 "$DATA_DIR/test_pe_R1.fq" \
    -2 "$DATA_DIR/test_pe_R2.fq" \
    -o "$TEST_DIR/baseline_pe.sam" \
    --min-num-seeds 1 --error-threshold 10 \
    > /dev/null 2>&1

# ============================================
# CURRENT RUNS
# ============================================
echo ""
echo "=============================================="
echo "Running current mappings..."
echo "=============================================="

# Single-end current (primary only)
echo "  Single-end (primary)..."
"$CURRENT_BIN" --SAM \
    -r "$DATA_DIR/test_ref.fa" \
    -x "$DATA_DIR/test_ref.idx" \
    -1 "$DATA_DIR/test_se.fq" \
    -o "$TEST_DIR/current_se.sam" \
    --min-num-seeds 1 --error-threshold 10 \
    > /dev/null 2>&1

# Single-end current (with split)
echo "  Single-end (with Y/NoY split)..."
"$CURRENT_BIN" --SAM --emit-noY-bam --emit-Y-bam \
    -r "$DATA_DIR/test_ref.fa" \
    -x "$DATA_DIR/test_ref.idx" \
    -1 "$DATA_DIR/test_se.fq" \
    -o "$TEST_DIR/current_se_split.sam" \
    --min-num-seeds 1 --error-threshold 10 \
    > /dev/null 2>&1

# Paired-end current (primary only)
echo "  Paired-end (primary)..."
"$CURRENT_BIN" --SAM \
    -r "$DATA_DIR/test_ref.fa" \
    -x "$DATA_DIR/test_ref.idx" \
    -1 "$DATA_DIR/test_pe_R1.fq" \
    -2 "$DATA_DIR/test_pe_R2.fq" \
    -o "$TEST_DIR/current_pe.sam" \
    --min-num-seeds 1 --error-threshold 10 \
    > /dev/null 2>&1

# Paired-end current (with split)
echo "  Paired-end (with Y/NoY split)..."
"$CURRENT_BIN" --SAM --emit-noY-bam --emit-Y-bam \
    -r "$DATA_DIR/test_ref.fa" \
    -x "$DATA_DIR/test_ref.idx" \
    -1 "$DATA_DIR/test_pe_R1.fq" \
    -2 "$DATA_DIR/test_pe_R2.fq" \
    -o "$TEST_DIR/current_pe_split.sam" \
    --min-num-seeds 1 --error-threshold 10 \
    > /dev/null 2>&1

# ============================================
# COMPARISONS
# ============================================
echo ""
echo "=============================================="
echo "Running comparisons..."
echo "=============================================="

FAILED=0

# Helper function to extract SAM body (no headers) and sort
extract_and_sort() {
    samtools view "$1" | sort
}

# Test 1: Primary SE output vs baseline
echo ""
echo "Test 1: Single-end primary output vs baseline"
extract_and_sort "$TEST_DIR/baseline_se.sam" > "$TEST_DIR/baseline_se.sorted.body"
extract_and_sort "$TEST_DIR/current_se.sam" > "$TEST_DIR/current_se.sorted.body"

if diff -q "$TEST_DIR/baseline_se.sorted.body" "$TEST_DIR/current_se.sorted.body" > /dev/null 2>&1; then
    echo -e "  ${GREEN}PASS${NC}: SE primary output matches baseline"
else
    echo -e "  ${RED}FAIL${NC}: SE primary output differs from baseline"
    diff "$TEST_DIR/baseline_se.sorted.body" "$TEST_DIR/current_se.sorted.body" | head -20
    FAILED=1
fi

# Test 2: Primary PE output vs baseline
echo ""
echo "Test 2: Paired-end primary output vs baseline"
extract_and_sort "$TEST_DIR/baseline_pe.sam" > "$TEST_DIR/baseline_pe.sorted.body"
extract_and_sort "$TEST_DIR/current_pe.sam" > "$TEST_DIR/current_pe.sorted.body"

if diff -q "$TEST_DIR/baseline_pe.sorted.body" "$TEST_DIR/current_pe.sorted.body" > /dev/null 2>&1; then
    echo -e "  ${GREEN}PASS${NC}: PE primary output matches baseline"
else
    echo -e "  ${RED}FAIL${NC}: PE primary output differs from baseline"
    diff "$TEST_DIR/baseline_pe.sorted.body" "$TEST_DIR/current_pe.sorted.body" | head -20
    FAILED=1
fi

# Test 3: SE split integrity - Y file has only chrY
echo ""
echo "Test 3: SE split integrity - Y file exclusivity"
Y_NON_CHRY=$(samtools view "$TEST_DIR/current_se_split.Y.sam" | awk '$3 != "chrY" && $3 != "*" {count++} END {print count+0}')
if [ "$Y_NON_CHRY" -eq 0 ]; then
    echo -e "  ${GREEN}PASS${NC}: SE Y file contains only chrY alignments"
else
    echo -e "  ${RED}FAIL${NC}: SE Y file contains $Y_NON_CHRY non-chrY alignments"
    FAILED=1
fi

# Test 4: SE split integrity - noY file has no chrY
echo ""
echo "Test 4: SE split integrity - noY file exclusivity"
NOY_CHRY=$(samtools view "$TEST_DIR/current_se_split.noY.sam" | awk '$3 == "chrY" {count++} END {print count+0}')
if [ "$NOY_CHRY" -eq 0 ]; then
    echo -e "  ${GREEN}PASS${NC}: SE noY file contains no chrY alignments"
else
    echo -e "  ${RED}FAIL${NC}: SE noY file contains $NOY_CHRY chrY alignments"
    FAILED=1
fi

# Test 5: SE count conservation
echo ""
echo "Test 5: SE count conservation (Y + noY == primary)"
SE_TOTAL=$(samtools view -c "$TEST_DIR/current_se_split.sam")
SE_Y=$(samtools view -c "$TEST_DIR/current_se_split.Y.sam")
SE_NOY=$(samtools view -c "$TEST_DIR/current_se_split.noY.sam")
SE_SUM=$((SE_Y + SE_NOY))

if [ "$SE_TOTAL" -eq "$SE_SUM" ]; then
    echo -e "  ${GREEN}PASS${NC}: SE count conserved: $SE_TOTAL = $SE_Y + $SE_NOY"
else
    echo -e "  ${RED}FAIL${NC}: SE count mismatch: $SE_TOTAL != $SE_Y + $SE_NOY"
    FAILED=1
fi

# Test 6: PE split integrity - Y file has only chrY
echo ""
echo "Test 6: PE split integrity - Y file exclusivity"
PE_Y_NON_CHRY=$(samtools view "$TEST_DIR/current_pe_split.Y.sam" | awk '$3 != "chrY" && $3 != "*" {count++} END {print count+0}')
if [ "$PE_Y_NON_CHRY" -eq 0 ]; then
    echo -e "  ${GREEN}PASS${NC}: PE Y file contains only chrY alignments"
else
    echo -e "  ${RED}FAIL${NC}: PE Y file contains $PE_Y_NON_CHRY non-chrY alignments"
    FAILED=1
fi

# Test 7: PE split integrity - noY file has no chrY
echo ""
echo "Test 7: PE split integrity - noY file exclusivity"
PE_NOY_CHRY=$(samtools view "$TEST_DIR/current_pe_split.noY.sam" | awk '$3 == "chrY" {count++} END {print count+0}')
if [ "$PE_NOY_CHRY" -eq 0 ]; then
    echo -e "  ${GREEN}PASS${NC}: PE noY file contains no chrY alignments"
else
    echo -e "  ${RED}FAIL${NC}: PE noY file contains $PE_NOY_CHRY chrY alignments"
    FAILED=1
fi

# Test 8: PE count conservation
echo ""
echo "Test 8: PE count conservation (Y + noY == primary)"
PE_TOTAL=$(samtools view -c "$TEST_DIR/current_pe_split.sam")
PE_Y=$(samtools view -c "$TEST_DIR/current_pe_split.Y.sam")
PE_NOY=$(samtools view -c "$TEST_DIR/current_pe_split.noY.sam")
PE_SUM=$((PE_Y + PE_NOY))

if [ "$PE_TOTAL" -eq "$PE_SUM" ]; then
    echo -e "  ${GREEN}PASS${NC}: PE count conserved: $PE_TOTAL = $PE_Y + $PE_NOY"
else
    echo -e "  ${RED}FAIL${NC}: PE count mismatch: $PE_TOTAL != $PE_Y + $PE_NOY"
    FAILED=1
fi

# Test 9: Determinism check
echo ""
echo "Test 9: Determinism check (run twice, compare md5)"

# Run SE split again
"$CURRENT_BIN" --SAM --emit-noY-bam --emit-Y-bam \
    -r "$DATA_DIR/test_ref.fa" \
    -x "$DATA_DIR/test_ref.idx" \
    -1 "$DATA_DIR/test_se.fq" \
    -o "$TEST_DIR/current_se_split_run2.sam" \
    --min-num-seeds 1 --error-threshold 10 \
    > /dev/null 2>&1

# Compare md5sums
MD5_RUN1=$(cat "$TEST_DIR/current_se_split.sam" "$TEST_DIR/current_se_split.Y.sam" "$TEST_DIR/current_se_split.noY.sam" | md5sum | cut -d' ' -f1)
MD5_RUN2=$(cat "$TEST_DIR/current_se_split_run2.sam" "$TEST_DIR/current_se_split_run2.Y.sam" "$TEST_DIR/current_se_split_run2.noY.sam" | md5sum | cut -d' ' -f1)

if [ "$MD5_RUN1" = "$MD5_RUN2" ]; then
    echo -e "  ${GREEN}PASS${NC}: Outputs are deterministic (md5: $MD5_RUN1)"
else
    echo -e "  ${RED}FAIL${NC}: Outputs are non-deterministic"
    echo "    Run 1 md5: $MD5_RUN1"
    echo "    Run 2 md5: $MD5_RUN2"
    FAILED=1
fi

# Test 10: Header consistency
echo ""
echo "Test 10: Header consistency across split files"
HEADER_MAIN=$(samtools view -H "$TEST_DIR/current_se_split.sam" | grep '^@SQ' | sort)
HEADER_Y=$(samtools view -H "$TEST_DIR/current_se_split.Y.sam" | grep '^@SQ' | sort)
HEADER_NOY=$(samtools view -H "$TEST_DIR/current_se_split.noY.sam" | grep '^@SQ' | sort)

if [ "$HEADER_MAIN" = "$HEADER_Y" ] && [ "$HEADER_MAIN" = "$HEADER_NOY" ]; then
    echo -e "  ${GREEN}PASS${NC}: Headers are consistent across all split files"
else
    echo -e "  ${RED}FAIL${NC}: Headers differ between split files"
    FAILED=1
fi

# ============================================
# SUMMARY
# ============================================
echo ""
echo "=============================================="
echo "Summary"
echo "=============================================="

echo ""
echo "Mapping statistics:"
echo "  Baseline SE: $(samtools view -c "$TEST_DIR/baseline_se.sam") alignments"
echo "  Current SE:  $(samtools view -c "$TEST_DIR/current_se.sam") alignments"
echo "  SE Y:        $SE_Y alignments"
echo "  SE noY:      $SE_NOY alignments"
echo ""
echo "  Baseline PE: $(samtools view -c "$TEST_DIR/baseline_pe.sam") alignments"
echo "  Current PE:  $(samtools view -c "$TEST_DIR/current_pe.sam") alignments"
echo "  PE Y:        $PE_Y alignments"
echo "  PE noY:      $PE_NOY alignments"

echo ""
if [ "$FAILED" -eq 0 ]; then
    echo -e "${GREEN}========================================${NC}"
    echo -e "${GREEN}All regression tests PASSED!${NC}"
    echo -e "${GREEN}========================================${NC}"
    exit 0
else
    echo -e "${RED}========================================${NC}"
    echo -e "${RED}Some regression tests FAILED!${NC}"
    echo -e "${RED}========================================${NC}"
    echo ""
    echo "To inspect outputs: KEEP_TEST_DIR=1 $0"
    exit 1
fi

