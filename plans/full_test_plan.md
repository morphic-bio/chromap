# Full BAM/CRAM Writer Test Runbook

This runbook validates BAM/CRAM output parity against SAM (including Y/noY streams),
exercises low-mem spill paths on small fixtures, and compares memory usage against
the pre-fork `/usr/local/bin/chromap`. The same script flow is reused for large
datasets with normal limits.

## Goals
- Verify BAM/CRAM outputs (all, Y-only, noY) match SAM at the read level.
- Validate BAM/CRAM indexes and basic integrity checks.
- Exercise low-mem spill paths on small fixtures.
- Compare max RSS against `/usr/local/bin/chromap` (1 thread, no low-mem).
- Reuse the same scripts on large datasets for end-to-end validation.

## Prereqs
- Build current chromap: `make -j` in repo root (binary at `./chromap`).
- `samtools` available in PATH.
- Reference FASTA for each dataset (and `.fai` index). For CRAM, the reference
  must be identical to the one used for alignment output.
- Chromap index built for the reference (`chromap -i -r ref.fa -o ref.index`).

## Fixture sources (morphic-atac-seq)
The reference paths and base parameters are defined in
`/home/lhhung/morphic-atac-seq` (see `README.md` and `runChromap.sh`). Use these
as fixture defaults:
- `REF_FASTA`: `/storage/scRNAseq_output/indices-110-44/chromap/genome.fa`
- `INDEX_FILE`: `/storage/scRNAseq_output/indices-110-44/chromap/index`
- `CHROM_SIZES`: `/storage/scRNAseq_output/indices-110-44/chromap/chrNameLength.txt`
- Small FASTQs: `/storage/ATAC-Seq-10000` and `/storage/ATAC-short-test`
- Large FASTQs: `/mnt/pikachu/NW-5-21/ATAC-Seq` (confirm case; some hosts use `/mnt/pikachu/NW-5-21/ATAC-seq`)
- Pipeline chromap flags: `--trim-adapters --min-frag-length 30 -l 2000 --remove-pcr-duplicates`

## Paths and variables
Set these before running the tests:
```sh
export CHROMAP_BIN=/mnt/pikachu/chromap/chromap
export BASELINE_BIN=/usr/local/bin/chromap
export FIXTURE_REPO=/home/lhhung/morphic-atac-seq
export DATA_SMALL1=/storage/ATAC-Seq-10000
export DATA_SMALL2=/storage/ATAC-short-test
export DATA_LARGE=/mnt/pikachu/NW-5-21/ATAC-Seq
export REF_FA=/storage/scRNAseq_output/indices-110-44/chromap/genome.fa
export REF_INDEX=/storage/scRNAseq_output/indices-110-44/chromap/index
export CHROM_SIZES=/storage/scRNAseq_output/indices-110-44/chromap/chrNameLength.txt
export CHROMAP_ARGS_BASE="--trim-adapters --min-frag-length 30 -l 2000 --remove-pcr-duplicates"
export OUT_ROOT=/mnt/pikachu/chromap/tests/out/bamwriter
```
Resolve read files explicitly per dataset:
```sh
export READ1=/path/to/read1.fastq.gz
export READ2=/path/to/read2.fastq.gz
```
The pairing and sample-name conventions used in the pipeline are implemented in
`/home/lhhung/morphic-atac-seq/runChromap.sh`; use that logic when selecting
`READ1`/`READ2` for a single-sample test.

## Forcing low-mem spill on small fixtures
Low-mem spill happens when the in-memory mapping buffer exceeds an internal
threshold. With tiny fixtures, it may not trigger. If needed, temporarily lower
the threshold in `src/chromap.h` (both single-end and paired-end blocks):
- Reduce `max_num_mappings_in_mem` from `1 << 30` / `1 << 29` to a smaller value
  (e.g., `1 << 24`) for a test build.
- Rebuild (`make -j`) and revert after testing.

## Test matrix (small fixtures)
For each small dataset (`$DATA_SMALL1`, `$DATA_SMALL2`), run:
- Output formats: SAM, BAM, CRAM
- Streams: all, noY, Y
- Modes: normal and low-mem

Use consistent flags across formats. Example (morphic pipeline flags, adjust as needed):
```sh
COMMON_ARGS="$CHROMAP_ARGS_BASE -x $REF_INDEX -r $REF_FA -1 $READ1 -2 $READ2 -t 1"
Y_STREAMS="--emit-noY-bam --emit-Y-bam"

# Normal mode
$CHROMAP_BIN $COMMON_ARGS --SAM  $Y_STREAMS -o $OUT_ROOT/small1/sam/output.sam
$CHROMAP_BIN $COMMON_ARGS --BAM  $Y_STREAMS -o $OUT_ROOT/small1/bam/output.bam --write-index
$CHROMAP_BIN $COMMON_ARGS --CRAM $Y_STREAMS -o $OUT_ROOT/small1/cram/output.cram --write-index

# Low-mem mode (note: --write-index is incompatible with --low-mem)
$CHROMAP_BIN $COMMON_ARGS --SAM  --low-mem $Y_STREAMS -o $OUT_ROOT/small1/sam_lowmem/output.sam
$CHROMAP_BIN $COMMON_ARGS --BAM  --low-mem $Y_STREAMS -o $OUT_ROOT/small1/bam_lowmem/output.bam
$CHROMAP_BIN $COMMON_ARGS --CRAM --low-mem $Y_STREAMS -o $OUT_ROOT/small1/cram_lowmem/output.cram
```
Notes:
- Use `-t 1` for deterministic output when doing strict parity diffs.
- For CRAM, `-r $REF_FA` is required and must match the output reference.
- Y/noY files default to `<output>.noY.<ext>` and `<output>.Y.<ext>` unless
  overridden with `--noY-output` / `--Y-output`.

## Parity checks (SAM vs BAM/CRAM)
Canonicalize all outputs to SAM, sort alignments, and diff.

1) Convert outputs to canonical SAM:
```sh
# SAM input (re-emit through samtools for canonical tag ordering)
samtools view -h -O SAM output.sam > output.canon.sam

# BAM input
samtools view -h -O SAM output.bam > output.canon.sam

# CRAM input (requires reference)
samtools view -h -T $REF_FA -O SAM output.cram > output.canon.sam
```

2) Split header and body, sort body for stable compare:
```sh
grep '^@SQ'  output.canon.sam | LC_ALL=C sort > output.header.sq.txt
grep -v '^@' output.canon.sam | LC_ALL=C sort > output.body.txt
```

3) Compare against SAM baseline:
```sh
diff -u sam.output.header.sq.txt bam.output.header.sq.txt
diff -u sam.output.body.txt   bam.output.body.txt
diff -u sam.output.header.sq.txt cram.output.header.sq.txt
diff -u sam.output.body.txt   cram.output.body.txt
```
Notes:
- If `--read-group` is enabled, BAM/CRAM will contain `@RG` and `RG:Z` tags that
  SAM output does not; either disable `--read-group` for parity tests or compare
  BAM vs CRAM directly for RG coverage.
- Repeat steps 1-3 for `output.noY.*` and `output.Y.*` streams, comparing each to
  the corresponding SAM stream (`output.noY.sam`, `output.Y.sam`).

4) Verify counts and Y/noY partitioning:
```sh
samtools view -c output.sam
samtools view -c output.noY.sam
samtools view -c output.Y.sam
# Expect: all == noY + Y
```

## Index checks
For normal mode (no low-mem) where `--write-index` is enabled:
```sh
# BAM uses CSI index (.csi), not BAI
samtools quickcheck -v output.bam output.bam.csi
# CRAM uses CRAI index (.crai)
samtools quickcheck -v output.cram output.cram.crai
```
For low-mem outputs, index manually after sorting:
```sh
samtools sort -o output.sorted.bam output.bam
samtools index output.sorted.bam
```
For CRAM:
```sh
samtools sort --reference $REF_FA -O CRAM -o output.sorted.cram output.cram
samtools index output.sorted.cram
```

## Memory comparison vs pre-fork binary
Run the same dataset with current chromap and `/usr/local/bin/chromap` (1 thread,
no low-mem) and compare Max RSS.

```sh
COMMON_ARGS="$CHROMAP_ARGS_BASE -x $REF_INDEX -r $REF_FA -1 $READ1 -2 $READ2 -t 1"

(/usr/bin/time -v $CHROMAP_BIN $COMMON_ARGS --SAM -o $OUT_ROOT/mem/new.sam) 2> $OUT_ROOT/mem/new.time
(/usr/bin/time -v $BASELINE_BIN $COMMON_ARGS --SAM -o $OUT_ROOT/mem/old.sam) 2> $OUT_ROOT/mem/old.time

rg "Maximum resident set size" $OUT_ROOT/mem/new.time $OUT_ROOT/mem/old.time
```
If `/usr/bin/time -v` is unavailable, use `-f "MAXRSS=%M"` and parse the value.

## Large dataset runs
After small fixtures pass:
- Reuse the same script with normal limits (no forced low-mem).
- Use `$DATA_LARGE` inputs and a dedicated output directory.
- Run parity checks and index checks; expect longer runtimes.

## Expected outputs
- Per run: `output.{sam,bam,cram}` plus `output.noY.*` and `output.Y.*`.
- Index files in normal mode: `.csi` for BAM (CSI format) and `.crai` for CRAM.
- Canonicalized SAM and diff logs for parity checks.
- Memory logs for current vs baseline binaries.

## Test Pass/Fail Criteria

| Test | Pass Condition |
|------|----------------|
| Parity | `@SQ` header diff empty, body diff empty (per stream) |
| Y/noY counts | `count(all) == count(noY) + count(Y)` |
| Index (BAM) | `.csi` exists, `samtools quickcheck` passes |
| Index (CRAM) | `.crai` exists, `samtools quickcheck` passes |
| Memory | New binary RSS within 110% of baseline |

## Expected Warnings

- CIGAR/sequence mismatch warnings: Some reads may have CIGAR operations that don't match sequence length. These are written as unmapped records to avoid indexing issues. This is expected behavior and not a failure.

## Troubleshooting
- CRAM errors about missing contigs: ensure `REF_FA` matches the output header
  and `samtools faidx` is up to date.
- Y/noY outputs empty: valid if dataset lacks chrY; counts should still satisfy
  all = noY + Y.
- Low-mem does not spill: lower `max_num_mappings_in_mem` temporarily and rebuild.
