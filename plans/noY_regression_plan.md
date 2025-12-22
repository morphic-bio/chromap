# Regression Plan: Y/NoY Split vs Baseline

## Goal
Ensure the Y/NoY split changes introduce no unintended regressions by comparing current outputs against a baseline produced with an older chromap commit.

## Baseline Setup
1) Identify a previous known-good commit (pre Y/NoY changes). Check it out in a clean workspace or build the binary from that commit.
2) Build baseline chromap binary (no Y/NoY split changes).
3) Prepare a small, deterministic test dataset (e.g., existing synthetic or downsampled set used in prior tests) with known Y reads:
   - Reference: small ref with chr1 and chrY.
   - Reads: mix of Y and non-Y (single-end and paired-end if applicable).

## Baseline Run
- Command (example for SAM):
  ```
  ./chromap_baseline --ref ref.fa --index index.idx --read1 R1.fq [--read2 R2.fq] \
    --SAM -o baseline.sam
  ```
- If BAM: pipe through samtools or use existing BAM output option.
- Generate sorted SAM/BAM if needed: `samtools sort -O sam -o baseline.sorted.sam baseline.sam`.
- Record summaries: `samtools idxstats`, `samtools flagstat`.

## Current Run (with Y/NoY)
- Build current chromap.
- Run on the same dataset twice:
  1) Primary (no split): same command as baseline to produce `current.sam` (or BAM).
  2) Split enabled: add `--emit-noY-bam --emit-Y-bam` and run both unsorted and sorted modes as applicable. Capture `_Y` and `_noY` outputs.
- Generate sorted SAMs if needed.
- Record summaries: `samtools idxstats`, `samtools flagstat` for primary, Y, and noY.

## Comparisons
1) **Primary output vs baseline**
   - Diff sorted SAMs (`diff -u baseline.sorted.sam current.sorted.sam`) or compare idxstats/flagstat summaries. They should match (allowing for header metadata differences).
2) **Split integrity**
   - Check Y/noY exclusivity: no chrY in noY; only chrY in Y.
   - Count check: reads(Y) + reads(noY) == reads(primary/current) == reads(baseline).
   - Optional: derive Y-only/noY SAMs from the baseline primary (by filtering on chrY) and compare sorted SAMs to the split outputs.
3) **Determinism**
   - Run the split twice and ensure outputs are identical (md5sum) to catch nondeterministic behavior.

## Automation
- Write a simple regression script (e.g., `tests/run_noY_regression.sh`) that:
  - Builds/runs baseline binary (or uses a prebuilt baseline path).
  - Builds/runs current binary.
  - Performs the diffs/summaries and reports pass/fail.
  - Uses a small synthetic/downsampled dataset to keep runtime low and deterministic.

## Success Criteria
- Primary current output matches baseline (no unintended changes).
- Split outputs satisfy exclusivity and count conservation and match filtered baseline partitions.
- Runs are deterministic (repeatable md5sum).
