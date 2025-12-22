# Implementation Plan: Three-Stream SAM Output (all / Y-only / noY)

## Overview

Add CLI flags to emit up to three SAM output streams when running in `--SAM` mode:
1. **Primary (all)**: Unchanged behavior - all mappings
2. **Y-only**: Only reads that have any alignment (primary or secondary) to Y contigs
3. **noY**: Reads that have NO alignments to Y contigs

Filtering is **read-ID based**: if any alignment of a read touches Y, ALL its alignments go to Y-only (and are excluded from noY).

For paired-end reads: if either mate maps to Y, the entire pair is filtered.

---

## Phase 1: CLI and Parameter Plumbing

### 1.1 Add CLI Flags

**File**: `src/chromap_driver.cc`

Add to `AddOutputOptions()`:

```cpp
("emit-noY-bam", "Emit additional SAM stream excluding Y-mapped reads (requires --SAM)")(
"noY-output", "Explicit path for noY output (default: <output>.noY.sam)",
    cxxopts::value<std::string>(), "FILE")(
"emit-Y-bam", "Emit additional SAM stream with only Y-mapped reads (requires --SAM)")(
"Y-output", "Explicit path for Y-only output (default: <output>.Y.sam)",
    cxxopts::value<std::string>(), "FILE");
```

Add validation after parsing:
```cpp
if (result.count("emit-noY-bam") || result.count("emit-Y-bam")) {
  if (mapping_parameters.mapping_output_format != MAPPINGFORMAT_SAM) {
    chromap::ExitWithMessage("--emit-noY-bam and --emit-Y-bam require --SAM mode");
  }
}
```

Add path derivation logic:
```cpp
if (result.count("emit-noY-bam")) {
  mapping_parameters.emit_noY_stream = true;
  if (result.count("noY-output")) {
    mapping_parameters.noY_output_path = result["noY-output"].as<std::string>();
  } else {
    // Auto-derive: output.sam -> output.noY.sam
    mapping_parameters.noY_output_path = DeriveSecondaryPath(
        mapping_parameters.mapping_output_file_path, ".noY.sam");
  }
}
// Similar for emit-Y-bam / Y-output
```

### 1.2 Extend MappingParameters

**File**: `src/mapping_parameters.h`

Add new fields:

```cpp
// Y-chromosome filtering
bool emit_noY_stream = false;
bool emit_Y_stream = false;
std::string noY_output_path;
std::string Y_output_path;
```

### 1.3 Path Derivation Helper

**File**: `src/chromap_driver.cc` (in anonymous namespace)

```cpp
std::string DeriveSecondaryPath(const std::string &primary_path, 
                                 const std::string &suffix) {
  // Handle /dev/stdout, /dev/stderr
  if (primary_path == "/dev/stdout" || primary_path == "/dev/stderr") {
    return "chromap_output" + suffix;
  }
  
  // Find extension and insert suffix before it
  size_t dot_pos = primary_path.rfind('.');
  if (dot_pos != std::string::npos && dot_pos > primary_path.rfind('/')) {
    return primary_path.substr(0, dot_pos) + suffix;
  }
  return primary_path + suffix;
}
```

---

## Phase 2: Y Contig Detection

### 2.1 Create Y RID Mask

**File**: `src/sequence_batch.h` (or new `src/y_contig_detector.h`)

Add a utility function to build a set of reference IDs (RIDs) that represent Y contigs:

```cpp
#include <unordered_set>
#include <cctype>
#include <algorithm>

inline std::unordered_set<uint32_t> BuildYContigRidMask(
    uint32_t num_reference_sequences,
    const SequenceBatch &reference) {
  
  std::unordered_set<uint32_t> y_rids;
  
  for (uint32_t rid = 0; rid < num_reference_sequences; ++rid) {
    std::string name(reference.GetSequenceNameAt(rid));
    
    // Convert to lowercase for case-insensitive matching
    std::transform(name.begin(), name.end(), name.begin(), ::tolower);
    
    // Strip optional "chr" prefix
    if (name.substr(0, 3) == "chr") {
      name = name.substr(3);
    }
    
    // Match "y" exactly (not "y_random", etc. unless we want those too)
    if (name == "y") {
      y_rids.insert(rid);
    }
  }
  
  if (y_rids.empty()) {
    std::cerr << "WARNING: No Y contigs found in reference. "
              << "Y-filtering will have no effect.\n";
  } else {
    std::cerr << "Found " << y_rids.size() << " Y contig(s) for filtering.\n";
  }
  
  return y_rids;
}
```

### 2.2 Integrate into Chromap Class

**File**: `src/chromap.h`

Add member variable:
```cpp
std::unordered_set<uint32_t> y_contig_rids_;
```

Initialize after loading reference (in `MapSingleEndReads` and `MapPairedEndReads`):
```cpp
if (mapping_parameters_.emit_noY_stream || mapping_parameters_.emit_Y_stream) {
  y_contig_rids_ = BuildYContigRidMask(num_reference_sequences, reference);
}
```

---

## Phase 3: Y-Hit Tracking During Mapping

### 3.1 Thread-Local Y-Hit Collection

Reads that have ANY alignment to Y need to be tracked. This happens in the mapping generator when SAMMapping records are created.

**Option A**: Track in mapping loop, collect read IDs with Y hits

Add thread-local vectors in the parallel mapping section:

**File**: `src/chromap.h` (in `MapSingleEndReads` / `MapPairedEndReads`)

```cpp
// Thread-local storage for Y-hit read IDs
std::vector<std::vector<uint32_t>> thread_y_hit_read_ids(mapping_parameters_.num_threads);
```

### 3.2 Detect Y Hits When Generating SAMMapping

**File**: `src/mapping_generator.h`

In `EmplaceBackSingleEndMappingRecord` and `EmplaceBackPairedEndMappingRecord` specializations for SAMMapping:

When a SAMMapping is created, check if `rid` (or `mrid` for paired-end) is in the Y mask:

```cpp
// After creating the SAMMapping, check for Y hits
if (!y_contig_rids_.empty()) {
  bool hits_y = y_contig_rids_.count(mapping_in_memory.rid) > 0;
  // For paired-end, also check mrid
  if (hits_y) {
    // Record this read_id as a Y-hitter
    // (will be merged into a set after parallel region)
  }
}
```

**Better approach**: Pass a reference to a thread-local vector and the Y mask into the generator, and record hits there.

### 3.3 Merge Thread-Local Y Hits

After the parallel mapping region, merge all thread-local Y-hit vectors into a single `std::unordered_set<uint32_t>`:

```cpp
std::unordered_set<uint32_t> reads_with_y_hit;
for (const auto &thread_vec : thread_y_hit_read_ids) {
  reads_with_y_hit.insert(thread_vec.begin(), thread_vec.end());
}
```

---

## Phase 4: Extend MappingWriter for Triple Output

### 4.1 Add Secondary File Handles

**File**: `src/mapping_writer.h`

Add new members to the `MappingWriter` class:

```cpp
protected:
  FILE *noY_output_file_ = nullptr;
  FILE *Y_output_file_ = nullptr;
  const std::unordered_set<uint32_t> *reads_with_y_hit_ = nullptr;
```

### 4.2 Modify Constructor

**File**: `src/mapping_writer.h` or `src/mapping_writer.cc`

```cpp
MappingWriter(const MappingParameters &mapping_parameters,
              const uint32_t cell_barcode_length,
              const std::vector<int> &pairs_custom_rid_rank,
              const std::unordered_set<uint32_t> *reads_with_y_hit = nullptr)
    : /* existing initializers */,
      reads_with_y_hit_(reads_with_y_hit) {
  
  // Existing primary file open...
  
  // Open secondary streams if configured
  if (mapping_parameters_.emit_noY_stream) {
    noY_output_file_ = fopen(mapping_parameters_.noY_output_path.c_str(), "w");
    assert(noY_output_file_ != nullptr);
  }
  if (mapping_parameters_.emit_Y_stream) {
    Y_output_file_ = fopen(mapping_parameters_.Y_output_path.c_str(), "w");
    assert(Y_output_file_ != nullptr);
  }
}
```

### 4.3 Modify Destructor

```cpp
~MappingWriter() {
  fclose(mapping_output_file_);
  if (noY_output_file_) fclose(noY_output_file_);
  if (Y_output_file_) fclose(Y_output_file_);
}
```

### 4.4 Mirror Headers to All Streams

**File**: `src/mapping_writer.cc`

Modify `MappingWriter<SAMMapping>::OutputHeader`:

```cpp
template <>
void MappingWriter<SAMMapping>::OutputHeader(uint32_t num_reference_sequences,
                                             const SequenceBatch &reference) {
  for (uint32_t rid = 0; rid < num_reference_sequences; ++rid) {
    const char *reference_sequence_name = reference.GetSequenceNameAt(rid);
    uint32_t reference_sequence_length = reference.GetSequenceLengthAt(rid);
    std::string header_line = "@SQ\tSN:" + std::string(reference_sequence_name) +
                              "\tLN:" + std::to_string(reference_sequence_length) + "\n";
    
    // Write to all streams
    this->AppendMappingOutput(header_line);
    if (noY_output_file_) {
      fwrite(header_line.data(), 1, header_line.size(), noY_output_file_);
    }
    if (Y_output_file_) {
      fwrite(header_line.data(), 1, header_line.size(), Y_output_file_);
    }
  }
}
```

### 4.5 Route Records to Appropriate Streams

**File**: `src/mapping_writer.cc`

Modify `MappingWriter<SAMMapping>::AppendMapping`:

```cpp
template <>
void MappingWriter<SAMMapping>::AppendMapping(uint32_t rid,
                                              const SequenceBatch &reference,
                                              const SAMMapping &mapping) {
  // Build the SAM record string (existing code)
  std::string out;
  // ... existing record building code ...
  
  // Always write to primary (all) stream
  this->AppendMappingOutput(out);
  
  // Route to secondary streams based on Y-hit status
  if (reads_with_y_hit_ && (noY_output_file_ || Y_output_file_)) {
    bool is_y_hit = reads_with_y_hit_->count(mapping.read_id_) > 0;
    
    if (Y_output_file_ && is_y_hit) {
      fwrite(out.data(), 1, out.size(), Y_output_file_);
    }
    if (noY_output_file_ && !is_y_hit) {
      fwrite(out.data(), 1, out.size(), noY_output_file_);
    }
  }
}
```

---

## Phase 5: Integration

### 5.1 Two-Pass Architecture

The challenge: we need to know ALL Y-hit read IDs before we can correctly route output. But mappings are generated in parallel and may be output in low-memory mode before we've seen all reads.

**Solution**: 

For non-low-memory mode:
1. Generate all mappings, collecting Y-hit read IDs
2. Build the final `reads_with_y_hit` set
3. Pass this set to the MappingWriter before output

For low-memory mode:
1. First pass: collect Y-hit read IDs (can be done during normal mapping)
2. When outputting from temp files, use the complete Y-hit set for filtering

### 5.2 Modify Chromap::MapPairedEndReads (and MapSingleEndReads)

**File**: `src/chromap.h`

After the parallel mapping region, before output:

```cpp
// Build Y-hit set from thread-local collections
std::unordered_set<uint32_t> reads_with_y_hit;
if (mapping_parameters_.emit_noY_stream || mapping_parameters_.emit_Y_stream) {
  for (const auto &thread_vec : thread_y_hit_read_ids) {
    reads_with_y_hit.insert(thread_vec.begin(), thread_vec.end());
  }
  std::cerr << "Found " << reads_with_y_hit.size() 
            << " reads with Y-chromosome alignments.\n";
}

// Pass to MappingWriter (may need to reconstruct or update writer)
```

---

## Phase 6: Edge Cases

### 6.1 /dev/stdout and /dev/stderr

Already handled in `DeriveSecondaryPath()` - derives to `chromap_output.noY.sam`.

### 6.2 Low-Memory Mode

The Y-hit set must be complete before final output. Since Y-hits are detected during mapping (not during temp file processing), the set should be complete when we reach the output phase.

### 6.3 Empty Y Set

If no Y contigs are found, warn the user but continue. The noY stream will be identical to the all stream, and the Y stream will be empty.

### 6.4 Secondary Alignments

Already handled by read-ID based filtering - if read X has a secondary alignment to Y, all of X's alignments (including primary on chr1) are filtered from noY.

---

## Phase 7: Testing

### 7.1 Unit Test: Y Mask Detection

**File**: `tests/test_y_detection.cc` (new)

```cpp
// Test cases:
// - "chrY", "Y", "CHR_Y" all detected
// - "chr1", "chrX", "Y_random" NOT detected  
// - Empty reference -> warning
// - Case insensitivity
```

### 7.2 Unit Test: MappingWriter Routing

**File**: `tests/test_mapping_writer_y_filter.cc` (new)

```cpp
// Mock SequenceBatch with contigs: chr1, chrY, chrM
// Create SAMMapping records:
//   - read1 on chr1 -> present in all, present in noY, absent from Y
//   - read2 on chrY -> present in all, absent from noY, present in Y
//   - read3 primary on chr1, secondary on chrY -> present in all, absent from noY, present in Y
//   - read4 paired: mate1 on chr1, mate2 on chrY -> both records in all and Y, neither in noY
// Verify headers identical in all outputs
```

### 7.3 Integration Test: Mapping Generator Y-Hit Collection

```cpp
// Drive MappingGenerator with synthetic MappingInMemory
// Verify thread-local Y-hit vectors collect correct read IDs
```

### 7.4 End-to-End Smoke Test

**File**: `tests/e2e/test_noy_output.sh` (new)

```bash
#!/bin/bash
# Create tiny reference with chr1 and chrY
# Create synthetic reads: some to chr1, some to chrY, some pairs split
# Run chromap --SAM --emit-noY-bam --emit-Y-bam
# Verify:
#   - samtools view parses all three outputs
#   - Y output contains only Y-mapped reads
#   - noY output contains only non-Y reads
#   - all = Y ∪ noY (by read count)
```

### 7.5 Low-Memory Mode Test

```bash
# Same as above but with --low-mem and threads > 1
# Verify filtering still works correctly after temp file merge
```

### 7.6 Path Derivation Test

```bash
# Test with -o /dev/stdout -> verify noY file created as chromap_output.noY.sam
# Test with explicit --noY-output path
```

---

## Implementation Order

1. **Phase 1**: CLI and parameter plumbing (can be tested independently)
2. **Phase 2**: Y contig detection (can be tested independently)  
3. **Phase 4**: MappingWriter triple output (can be stubbed with empty Y-hit set)
4. **Phase 3**: Y-hit tracking in mapping generator
5. **Phase 5**: Integration wiring
6. **Phase 6**: Edge case handling
7. **Phase 7**: Testing throughout

---

## Files Modified Summary

| File | Changes |
|------|---------|
| `src/chromap_driver.cc` | Add CLI flags, validation, path derivation |
| `src/mapping_parameters.h` | Add new fields for Y-filtering config |
| `src/mapping_writer.h` | Add secondary file handles, Y-hit set pointer |
| `src/mapping_writer.cc` | Triple header output, filtered record routing |
| `src/mapping_generator.h` | Detect Y hits during SAMMapping creation |
| `src/chromap.h` | Y contig mask init, thread-local Y-hit collection, integration |
| `src/sequence_batch.h` | (optional) Add Y mask builder utility |

---

## Risks and Mitigations

| Risk | Mitigation |
|------|------------|
| Memory overhead from Y-hit set | Set uses uint32_t read IDs; even 100M reads = ~400MB |
| Thread safety of Y-hit collection | Use thread-local vectors, merge after parallel region |
| Low-memory mode complexity | Y-hits collected during mapping phase, available for output phase |
| Performance impact | Filtering is O(1) hash lookup per record; negligible |

