#!/bin/bash
# Main test runner for BAM/CRAM writer test suite
# Usage: ./run_tests.sh [options]

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"
source "$SCRIPT_DIR/lib/utils.sh"

# Parse options
SMALL_ONLY=false
LARGE_ONLY=false
SKIP_LOWMEM=false
SKIP_MEMORY=false
DATASET=""
FORMAT=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --small-only)
            SMALL_ONLY=true
            shift
            ;;
        --large-only)
            LARGE_ONLY=true
            shift
            ;;
        --skip-lowmem)
            SKIP_LOWMEM=true
            shift
            ;;
        --skip-memory)
            SKIP_MEMORY=true
            shift
            ;;
        --dataset)
            DATASET="$2"
            shift 2
            ;;
        --format)
            FORMAT="$2"
            shift 2
            ;;
        -h|--help)
            echo "Usage: $0 [options]"
            echo "Options:"
            echo "  --small-only     Run only small fixture tests"
            echo "  --large-only     Run only large dataset tests"
            echo "  --skip-lowmem    Skip low-mem mode tests"
            echo "  --skip-memory    Skip memory comparison"
            echo "  --dataset NAME   Run single dataset (small1, small2, large)"
            echo "  --format FMT     Run single format (sam, bam, cram)"
            exit 0
            ;;
        *)
            log_error "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Check prerequisites
log_info "Checking prerequisites..."
if ! check_prereqs; then
    log_error "Prerequisites check failed"
    exit 1
fi

# Setup output directories
mkdir -p "$OUT_ROOT"/{reports,memory}

# Track test results
TOTAL_TESTS=0
PASSED_TESTS=0
FAILED_TESTS=0

# Function to run a single test
run_test() {
    local dataset_name="$1"
    local format="$2"
    local mode="$3"  # normal or lowmem
    local dataset_dir="$4"
    local read1="$5"
    local read2="$6"
    
    TOTAL_TESTS=$((TOTAL_TESTS + 1))
    
    local output_dir="$OUT_ROOT/$dataset_name/${format}"
    if [ "$mode" = "lowmem" ]; then
        output_dir="$OUT_ROOT/$dataset_name/${format}_lowmem"
    fi
    
    mkdir -p "$output_dir"
    
    local output_file="$output_dir/output.${format}"
    local format_flag="--${format^^}"
    
    # Build command
    local cmd="$CHROMAP_BIN $CHROMAP_ARGS_BASE -x $REF_INDEX -r $REF_FA"
    cmd="$cmd -1 $read1 -2 $read2 -t $NUM_THREADS"
    cmd="$cmd $format_flag $Y_STREAMS"
    
    if [ "$mode" = "lowmem" ]; then
        cmd="$cmd --low-mem"
    else
        if [ "$format" != "sam" ]; then
            # --write-index now requires --sort-bam
            cmd="$cmd --sort-bam --write-index"
            if [ "$format" = "cram" ]; then
                cmd="$cmd --hts-threads $HTS_THREADS"
            fi
        fi
    fi
    
    cmd="$cmd -o $output_file"
    
    log_info "Running: $dataset_name / $format / $mode"
    
    # Run chromap
    if eval "$cmd" > "$output_dir/run.log" 2>&1; then
        log_success "chromap completed: $output_file"
    else
        log_error "chromap failed: $output_file"
        FAILED_TESTS=$((FAILED_TESTS + 1))
        return 1
    fi
    
    # Verify output files exist
    if [ ! -f "$output_file" ]; then
        log_error "Output file not created: $output_file"
        FAILED_TESTS=$((FAILED_TESTS + 1))
        return 1
    fi
    
    # Check Y/noY streams
    # Y/noY files use pattern: <base>.noY.<ext> and <base>.Y.<ext>
    local output_base="${output_file%.${format}}"
    local noy_file="${output_base}.noY.${format}"
    local y_file="${output_base}.Y.${format}"
    
    if [ ! -f "$noy_file" ] || [ ! -f "$y_file" ]; then
        log_error "Y/noY streams not created"
        FAILED_TESTS=$((FAILED_TESTS + 1))
        return 1
    fi
    
    # Verify counts
    local all_count=$(count_reads "$output_file" "$REF_FA")
    local noy_count=$(count_reads "$noy_file" "$REF_FA")
    local y_count=$(count_reads "$y_file" "$REF_FA")
    local expected_count=$((noy_count + y_count))
    
    if [ "$all_count" != "$expected_count" ]; then
        log_error "Count mismatch: all=$all_count, noY=$noy_count, Y=$y_count, expected=$expected_count"
        FAILED_TESTS=$((FAILED_TESTS + 1))
        return 1
    fi
    
    # Run parity check (compare against SAM baseline)
    # For lowmem mode, sort outputs before comparison to validate k-way merge correctness
    if [ "$format" != "sam" ] && [ -f "$OUT_ROOT/$dataset_name/sam/output.sam" ]; then
        log_info "Running parity check..."
        local sort_flag="false"
        if [ "$mode" = "lowmem" ]; then
            sort_flag="true"
            log_info "Lowmem mode: will sort outputs before comparison"
        fi
        if "$SCRIPT_DIR/check_parity.sh" "$OUT_ROOT/$dataset_name/sam/output.sam" "$output_file" "$REF_FA" "$sort_flag" >> "$output_dir/parity.log" 2>&1; then
            log_success "Parity check passed"
        else
            log_error "Parity check failed"
            FAILED_TESTS=$((FAILED_TESTS + 1))
            return 1
        fi
    fi
    
    # Verify lowmem spill was triggered (for lowmem mode)
    if [ "$mode" = "lowmem" ]; then
        log_info "Checking if lowmem overflow was triggered..."
        if grep -q -E "(overflow files|temp files|Processing.*overflow|k-way merge)" "$output_dir/run.log" 2>/dev/null; then
            log_success "Lowmem overflow detected in logs"
        else
            log_warning "No overflow evidence found - spill may not have been triggered (dataset may be too small)"
            # This is not a failure, just informational
        fi
    fi
    
    # Run index check (for normal mode with indexing)
    if [ "$mode" = "normal" ] && [ "$format" != "sam" ]; then
        log_info "Running index check..."
        if "$SCRIPT_DIR/check_index.sh" "$output_file" "$REF_FA" >> "$output_dir/index.log" 2>&1; then
            log_success "Index check passed"
        else
            log_error "Index check failed"
            FAILED_TESTS=$((FAILED_TESTS + 1))
            return 1
        fi
    fi
    
    PASSED_TESTS=$((PASSED_TESTS + 1))
    return 0
}

# Function to run tests for a dataset
run_dataset() {
    local dataset_name="$1"
    local dataset_dir="$2"
    
    log_info "Processing dataset: $dataset_name"
    
    # Resolve read files
    resolve_reads "$dataset_dir" READ1 READ2 || {
        log_error "Failed to resolve reads for $dataset_name"
        return 1
    }
    
    # Determine if this is a small dataset (for low-mem testing)
    local is_small=false
    if [ "$dataset_name" = "small1" ] || [ "$dataset_name" = "small2" ]; then
        is_small=true
    fi
    
    # Run formats
    for format in sam bam cram; do
        if [ -n "$FORMAT" ] && [ "$FORMAT" != "$format" ]; then
            continue
        fi
        
        # Normal mode
        run_test "$dataset_name" "$format" "normal" "$dataset_dir" "$READ1" "$READ2"
        
        # Low-mem mode (only for small datasets)
        if [ "$is_small" = true ] && [ "$SKIP_LOWMEM" = false ]; then
            run_test "$dataset_name" "$format" "lowmem" "$dataset_dir" "$READ1" "$READ2"
        fi
    done
}

# Main test execution
log_info "Starting BAM/CRAM writer test suite"
log_info "Output root: $OUT_ROOT"

# Run small datasets
if [ "$LARGE_ONLY" = false ]; then
    if [ -z "$DATASET" ] || [ "$DATASET" = "small1" ]; then
        run_dataset "small1" "$DATA_SMALL1"
    fi
    
    if [ -z "$DATASET" ] || [ "$DATASET" = "small2" ]; then
        run_dataset "small2" "$DATA_SMALL2"
    fi
fi

# Run large dataset
if [ "$SMALL_ONLY" = false ]; then
    if [ -z "$DATASET" ] || [ "$DATASET" = "large" ]; then
        run_dataset "large" "$DATA_LARGE"
    fi
fi

# Memory comparison
if [ "$SKIP_MEMORY" = false ]; then
    log_info "Running memory comparison..."
    if "$SCRIPT_DIR/compare_memory.sh" >> "$OUT_ROOT/reports/memory.log" 2>&1; then
        log_success "Memory comparison passed"
        PASSED_TESTS=$((PASSED_TESTS + 1))
    else
        log_warning "Memory comparison failed or skipped"
        FAILED_TESTS=$((FAILED_TESTS + 1))
    fi
    TOTAL_TESTS=$((TOTAL_TESTS + 1))
fi

# Sorting tests (only for small datasets to keep runtime reasonable)
if ([ "$SMALL_ONLY" = false ] && [ -z "$DATASET" ]) || [ "$DATASET" = "small1" ]; then
    log_info "Running sorting tests..."
    
    local sort_test_dir="$OUT_ROOT/small1/sort_tests"
    mkdir -p "$sort_test_dir"
    
    # Test 1: Basic sorting
    TOTAL_TESTS=$((TOTAL_TESTS + 1))
    log_info "Test: Basic sorting (--sort-bam)"
    local sorted_bam="$sort_test_dir/sorted.bam"
    if $CHROMAP_BIN $CHROMAP_ARGS_BASE -x $REF_INDEX -r $REF_FA \
        -1 "$DATA_SMALL1/read1.fq" -2 "$DATA_SMALL1/read2.fq" \
        -t 1 --hts-threads 1 --BAM --sort-bam -o "$sorted_bam" \
        > "$sort_test_dir/basic_sort.log" 2>&1; then
        if python3 "$SCRIPT_DIR/check_chromap_sorted.py" "$sorted_bam" >> "$sort_test_dir/basic_sort.log" 2>&1; then
            log_success "Basic sorting test passed"
            PASSED_TESTS=$((PASSED_TESTS + 1))
        else
            log_error "Basic sorting test failed (not sorted)"
            FAILED_TESTS=$((FAILED_TESTS + 1))
        fi
    else
        log_error "Basic sorting test failed (chromap error)"
        FAILED_TESTS=$((FAILED_TESTS + 1))
    fi
    
    # Test 2: Sorting with index
    TOTAL_TESTS=$((TOTAL_TESTS + 1))
    log_info "Test: Sorting with index (--sort-bam --write-index)"
    local sorted_indexed_bam="$sort_test_dir/sorted_indexed.bam"
    if $CHROMAP_BIN $CHROMAP_ARGS_BASE -x $REF_INDEX -r $REF_FA \
        -1 "$DATA_SMALL1/read1.fq" -2 "$DATA_SMALL1/read2.fq" \
        -t 1 --hts-threads 1 --BAM --sort-bam --write-index -o "$sorted_indexed_bam" \
        > "$sort_test_dir/sort_index.log" 2>&1; then
        if [ -f "${sorted_indexed_bam}.csi" ] && \
           python3 "$SCRIPT_DIR/check_chromap_sorted.py" "$sorted_indexed_bam" >> "$sort_test_dir/sort_index.log" 2>&1 && \
           "$SCRIPT_DIR/check_index.sh" "$sorted_indexed_bam" "$REF_FA" >> "$sort_test_dir/sort_index.log" 2>&1; then
            log_success "Sorting with index test passed"
            PASSED_TESTS=$((PASSED_TESTS + 1))
        else
            log_error "Sorting with index test failed"
            FAILED_TESTS=$((FAILED_TESTS + 1))
        fi
    else
        log_error "Sorting with index test failed (chromap error)"
        FAILED_TESTS=$((FAILED_TESTS + 1))
    fi
    
    # Test 3: Sorting with low RAM (spill test)
    TOTAL_TESTS=$((TOTAL_TESTS + 1))
    log_info "Test: Sorting with low RAM (--sort-bam --sort-bam-ram 10M)"
    local sorted_lowram_bam="$sort_test_dir/sorted_lowram.bam"
    if $CHROMAP_BIN $CHROMAP_ARGS_BASE -x $REF_INDEX -r $REF_FA \
        -1 "$DATA_SMALL1/read1.fq" -2 "$DATA_SMALL1/read2.fq" \
        -t 1 --hts-threads 1 --BAM --sort-bam --sort-bam-ram 10M -o "$sorted_lowram_bam" \
        > "$sort_test_dir/sort_lowram.log" 2>&1; then
        if python3 "$SCRIPT_DIR/check_chromap_sorted.py" "$sorted_lowram_bam" >> "$sort_test_dir/sort_lowram.log" 2>&1; then
            # Check if spill was triggered (look for spill file mentions)
            if grep -q -E "(spill|chromap_sort_spill)" "$sort_test_dir/sort_lowram.log" 2>/dev/null; then
                log_success "Sorting with low RAM test passed (spill detected)"
            else
                log_success "Sorting with low RAM test passed (spill may not have been needed)"
            fi
            PASSED_TESTS=$((PASSED_TESTS + 1))
        else
            log_error "Sorting with low RAM test failed (not sorted)"
            FAILED_TESTS=$((FAILED_TESTS + 1))
        fi
    else
        log_error "Sorting with low RAM test failed (chromap error)"
        FAILED_TESTS=$((FAILED_TESTS + 1))
    fi
    
    # Test 4: Determinism test
    TOTAL_TESTS=$((TOTAL_TESTS + 1))
    log_info "Test: Deterministic output (--sort-bam with fixed threads)"
    local sorted1="$sort_test_dir/deterministic1.bam"
    local sorted2="$sort_test_dir/deterministic2.bam"
    if $CHROMAP_BIN $CHROMAP_ARGS_BASE -x $REF_INDEX -r $REF_FA \
        -1 "$DATA_SMALL1/read1.fq" -2 "$DATA_SMALL1/read2.fq" \
        -t 1 --hts-threads 1 --BAM --sort-bam -o "$sorted1" \
        > "$sort_test_dir/deterministic.log" 2>&1 && \
       $CHROMAP_BIN $CHROMAP_ARGS_BASE -x $REF_INDEX -r $REF_FA \
        -1 "$DATA_SMALL1/read1.fq" -2 "$DATA_SMALL1/read2.fq" \
        -t 1 --hts-threads 1 --BAM --sort-bam -o "$sorted2" \
        >> "$sort_test_dir/deterministic.log" 2>&1; then
        if cmp -s "$sorted1" "$sorted2"; then
            log_success "Determinism test passed (byte-identical output)"
            PASSED_TESTS=$((PASSED_TESTS + 1))
        else
            log_error "Determinism test failed (outputs differ)"
            FAILED_TESTS=$((FAILED_TESTS + 1))
        fi
    else
        log_error "Determinism test failed (chromap error)"
        FAILED_TESTS=$((FAILED_TESTS + 1))
    fi
    
    # Test 5: Y/noY routing with sorting
    TOTAL_TESTS=$((TOTAL_TESTS + 1))
    log_info "Test: Y/noY routing with sorting"
    local sorted_y_bam="$sort_test_dir/sorted_y.bam"
    if $CHROMAP_BIN $CHROMAP_ARGS_BASE -x $REF_INDEX -r $REF_FA \
        -1 "$DATA_SMALL1/read1.fq" -2 "$DATA_SMALL1/read2.fq" \
        -t 1 --hts-threads 1 --BAM --sort-bam --emit-noY-bam --emit-Y-bam -o "$sorted_y_bam" \
        > "$sort_test_dir/sort_y.log" 2>&1; then
        local sorted_y_noy="${sorted_y_bam%.bam}.noY.bam"
        local sorted_y_y="${sorted_y_bam%.bam}.Y.bam"
        if [ -f "$sorted_y_noy" ] && [ -f "$sorted_y_y" ] && \
           python3 "$SCRIPT_DIR/check_chromap_sorted.py" "$sorted_y_bam" >> "$sort_test_dir/sort_y.log" 2>&1 && \
           python3 "$SCRIPT_DIR/check_chromap_sorted.py" "$sorted_y_noy" >> "$sort_test_dir/sort_y.log" 2>&1 && \
           python3 "$SCRIPT_DIR/check_chromap_sorted.py" "$sorted_y_y" >> "$sort_test_dir/sort_y.log" 2>&1; then
            log_success "Y/noY routing with sorting test passed"
            PASSED_TESTS=$((PASSED_TESTS + 1))
        else
            log_error "Y/noY routing with sorting test failed"
            FAILED_TESTS=$((FAILED_TESTS + 1))
        fi
    else
        log_error "Y/noY routing with sorting test failed (chromap error)"
        FAILED_TESTS=$((FAILED_TESTS + 1))
    fi
fi

# Summary
log_info "Test suite completed"
echo ""
echo "=========================================="
echo "Test Summary"
echo "=========================================="
echo "Total tests:  $TOTAL_TESTS"
echo "Passed:       $PASSED_TESTS"
echo "Failed:       $FAILED_TESTS"
echo "=========================================="

if [ $FAILED_TESTS -eq 0 ]; then
    log_success "All tests passed!"
    exit 0
else
    log_error "$FAILED_TESTS test(s) failed"
    exit 1
fi

