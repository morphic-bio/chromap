# Changelog

All notable changes to this project will be documented in this file.

## [Unreleased]

### Fixed
- SAM writer: emit each record atomically and use length-safe `fwrite` to prevent interleaving and stray bytes.
- Low-memory spill: serialize flush section with `omp critical` and protect memory counter with `omp atomic` to avoid races and temp file name collisions.
- Temp mapping reload: correctly read `size_t` counts into `uint32_t` fields via local `size_t` then cast, preventing in-memory corruption during merges.
- New overflow system: eliminated remaining thread-safety issues in temp file handling by using thread-local writers with coordinated cleanup, preventing malformed SAM lines from unclosed/unflushed temp files.
- Serialization size calculation: fixed `SerializedSize()` methods to account for all variable-length fields, preventing malloc corruption during buffer operations.

### Added
- New CLI flag: `--temp-dir DIR` to specify directory for temporary files.
- New overflow system (compile with `NEW_OVERFLOW=1`): thread-safe temp file handling with simplified I/O patterns.
- Enhanced mapping record serialization with precise size calculation and single-write operations.
- Validation scripts under `scripts/`:
  - `validate_sam_fix.sh` for SAM path sanity checks.
  - `validate_low_mem_fix.sh` for low-memory path summary validation.
  - `test_overflow_basic.sh` for new overflow system integration testing.
- Docker support with multi-stage builds and temp directory improvements.

### Changed
- Temp file operations now use single atomic writes instead of multiple small writes with seeks.
- Thread cleanup is now explicit and coordinated rather than relying on automatic thread-local storage destruction.

### Notes
- Existing CLI flags and behavior remain backward compatible.
- Recommended to compile with `NEW_OVERFLOW=1` for improved robustness.
