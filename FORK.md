# Fork Notes: Chromap

This document records detailed changes in this fork, the rationale behind them, and notes for future work.

## Why this fork

- Address rare but impactful output corruption in SAM mode under high concurrency.
- Stabilize and harden the low-memory spill/merge path used for large datasets.

## Major fixes and how they were implemented

### 1) SAM writer: prevent line tearing and unsafe string writes
- Problem: SAM records were emitted via multiple `fprintf()` calls using `%s` with `std::string::data()`. Under concurrency, writes could interleave and `%s` relied on a terminating NUL (not guaranteed pre-C++17), risking over-read.
- Changes:
  - `src/mapping_writer.h`: `AppendMappingOutput` now uses `fwrite(line.data(), 1, line.size(), ...)` (length-safe).
  - `src/mapping_writer.cc`: `MappingWriter<SAMMapping>::AppendMapping` builds the entire record into a single `std::string` and writes it once (atomic at the `FILE*` level).
- Impact: Eliminates interleaved/partial lines and spurious bytes.

### 2) Low-memory spill path: serialize flush and protect counters
- Problem: Output tasks updated `num_mappings_in_mem`, sorted/cleared shared containers, and appended to `temp_mapping_file_handles` concurrently.
- Changes in `src/chromap.h` (both single-end and paired-end paths):
  - Use `#pragma omp atomic` for `num_mappings_in_mem += added_mappings`.
  - Wrap the flush sequence (threshold check → sort → spill → bookkeeping) in `#pragma omp critical(output_flush)`.
  - Temp file naming remains based on vector size but is now inside the critical region, ensuring uniqueness.
- Impact: Prevents data races, container corruption, and temp file name collisions.

### 3) Temp mapping reload: fix size_t → uint32_t read mismatch
- Problem: Spill files store counts as `size_t`. Reload path read directly into a `uint32_t` field, writing 8 bytes into a 4-byte slot on LP64 systems, corrupting adjacent state and causing huge sizes (e.g., 141012992) and truncations during merge.
- Changes in `src/temp_mapping.h`:
  - Read into a local `size_t`, then cast-assign to the `uint32_t` field.
- Impact: Correct in-memory state, preventing downstream corruption.
- This still has problems so we rewrote this and simplified the serialization and spillover

### 4) New overflow system: thread-local writers with coordinated cleanup
- Problem: The original temp file system had subtle race conditions during thread-local storage cleanup, causing some temp files to remain unclosed/unflushed when the main thread tried to read them, leading to remaining malformed SAM lines.
- Changes (enabled with `NEW_OVERFLOW=1` compile flag):
  - `src/overflow_writer.h/cc`: Thread-safe overflow writer using length-prefixed binary format
  - `src/mapping_writer.h/cc`: Each thread owns its own `OverflowWriter` instance during spill operations
  - `src/chromap.h`: Added explicit thread cleanup phase where each thread closes its writer and contributes file paths to shared collection
  - All mapping record types: Added precise `SerializedSize()` methods and single-write `WriteToFile()` implementations
- Impact: Eliminates remaining thread-safety issues, ensures all temp files are properly closed before reading, removes complex seek operations for simpler and more robust I/O.

## Tests and validation

- Scripts (run from repo root):
  - `./scripts/validate_sam_fix.sh`: sanity checks for SAM emission and basic execution.
  - `./scripts/validate_low_mem_fix.sh`: summarizes low-mem concurrency hardening.
  - `./scripts/test_overflow_basic.sh`: validates new overflow system integration.
- Recommended runtime checks:
  - SAM integrity: `awk 'BEGIN{FS="\t"} !/^@/ && NF<11{bad++} END{exit bad!=0}' output.sam`
  - SAM parsing: `samtools view -S output.sam > /dev/null`
- Stress testing: enable low-memory mode with a small threshold and high thread count.

## User-facing behavior and flags

- New flag: `--temp-dir DIR` to specify directory for temporary files (useful for Docker environments and custom temp locations).
- New compile flag: `NEW_OVERFLOW=1` to enable the improved overflow system (recommended for production use).
- Existing presets and options work as before.
- Performance: neutral or slightly improved due to fewer stdio calls per record and simplified I/O patterns.

## Compatibility

- Output formats unchanged. SAM/BED/PAF/PAIRS semantics preserved.
- Binary compatible with previous indices and configs.

## Future work

- Make `NEW_OVERFLOW=1` the default build configuration after validation.
- Consider a dedicated output thread to fully decouple compute from flush operations.
- Add CI checks for SAM field counts and `samtools view` parsing on sample outputs.
- Optional: add `setvbuf` to increase IO buffer size for throughput on very large runs.
- Docker integration improvements for better container-friendly temp directory handling.
