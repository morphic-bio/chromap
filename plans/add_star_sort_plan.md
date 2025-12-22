# Plan

Add an internal coordinate‑sorting path to chromap that mirrors the STAR‑Flex
sorter behavior and sorting spec (coordinate keys + deterministic QNAME
tie‑break), with Y/noY split support and index generation, while keeping an
optional samtools‑pipe fallback.

## Requirements
- Implement internal coordinate sorting with spill‑to‑disk and k‑way merge.
- Deterministic tie‑break (QNAME) when numeric keys are equal.
- Compatible with Y/noY streams and CRAM/BAM indexing.
- Document sort spec and expected behavior (not bit‑for‑bit samtools parity).
- Keep pipe‑to‑samtools as optional user path.

## Scope
- In: new sorter implementation, integration into mapping output path, tests.
- Out: full samtools parity guarantee, UI/CLI redesign beyond adding a sorting
  mode flag.

## Files and entry points
- `src/mapping_writer.*` (sorter integration + output flow)
- `src/chromap_driver.cc` (new flag / mode selection)
- New sorter files (e.g., `src/internal_sorter.*`)
- `tests/bamwriter/` (tests for deterministic sorted output and Y/noY split)
- `README.md` / `plans/bamwriter_plan.md` (behavior notes)

## Data model / API changes
- Add a sorting mode flag (e.g., `--bam-sort-method internal|samtools`).
- Record and expose sort spec in documentation.

## Action items
[ ] Extract STAR‑Flex sorter spec and adapt key structure (tid, pos, flag, mtid,
    mpos, isize, QNAME).
[ ] Implement internal sorter with spill‑to‑disk + k‑way merge, deterministic
    ordering.
[ ] Integrate sorter into chromap output path for BAM/CRAM.
[ ] Ensure Y/noY routing preserves the sort order spec.
[ ] Wire indexing for sorted outputs and guard against low‑mem incompatibilities.
[ ] Add tests for coordinate sortedness, deterministic ordering, Y/noY split
    counts, and index validity.
[ ] Update docs to describe ordering spec and optional samtools pipe.

## Testing and validation
- Small fixture tests: sortedness checks, deterministic ordering across runs,
  Y/noY counts.
- Regression checks vs samtools: allow ordering differences, but confirm
  identical alignment content.
- Large dataset run for scale and spill paths.

## Risks and edge cases
- Spill sorting stability and tie‑break correctness.
- Unmapped reads ordering (tid/pos sentinel handling).
- Indexing with CRAM reference requirements.
- Memory accounting vs large datasets.

## Notes

- use --sort-bam as flag
- byte-wise comparison of QNAME
