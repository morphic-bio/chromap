# Runbook

Integrate direct BAM/CRAM output via htslib with expected metadata/streaming features, while preserving existing SAM/BED/PAF/pairs behavior.

## Goals
- Add direct BAM and CRAM output using htslib, without altering existing SAM output logic or mapping behavior.
- Support expected htslib features: @HD/@SQ/@PG headers, optional @RG, stdout/pipe output, multi-threaded compression, and optional indexing.
- Keep low-memory and overflow paths working with the new writer.

## Scope
- In: htslib-backed BAM/CRAM output and indexing; CLI flags and docs.
- Out: Sorting in htslib/samtools, CRAM reference caching policies, new alignment tags beyond current SAM fields.

## Inputs and prerequisites
- htslib dev headers and library available on build host (system `-lhts`), or vendored build if required.
- Reference/index for testing:
  - `REF_FASTA=/storage/scRNAseq_output/indices-110-44/chromap/genome.fa`
  - `INDEX=/storage/scRNAseq_output/indices-110-44/chromap/index`
- Small FASTQs for quick tests:
  - `/storage/ATAC-short-test/ATAC-CHD4-0h-2_S2_R1_001.fastq.gz`
  - `/storage/ATAC-short-test/ATAC-CHD4-0h-2_S2_R2_001.fastq.gz`
- Larger FASTQs to force spill (if needed):
  - `/mnt/pikachu/NW-5-21/ATAC-Seq/ATAC-HDAC3-24h-3_S20_R1_001.fastq.gz`
  - `/mnt/pikachu/NW-5-21/ATAC-Seq/ATAC-HDAC3-24h-3_S20_R2_001.fastq.gz`

## Implementation steps
1) **Add output formats and CLI flags**
   - Extend `MappingOutputFormat` with `MAPPINGFORMAT_BAM` and `MAPPINGFORMAT_CRAM` in `src/mapping_parameters.h`.
   - Add CLI flags in `src/chromap_driver.cc`: `--BAM`, `--CRAM`, and `--write-index` (avoid `-i`).
   - Map BAM/CRAM to `Map*<SAMMapping>` in `src/chromap_driver.cc` and `src/chromap.h`.
   - Optionally auto-select based on output extension (`.bam`, `.cram`) if desired.

2) **Integrate htslib writer (SAM specialization)**
   - In `src/mapping_writer.cc`, add a BAM/CRAM output path inside `MappingWriter<SAMMapping>`:
     - Open `samFile*` with `sam_open(path, "wb")` or `sam_open(path, "wc")`.
     - Build `bam_hdr_t` with `@HD`, `@SQ` (reuse current header generation), and `@PG`.
     - Optional `@RG` and per-record `RG:Z` tag if a CLI flag provides read-group ID.
     - For CRAM, set reference with `hts_set_fai_filename` or `hts_set_opt(HTS_OPT_REFERENCE, ref_path)`.
     - Set compression threads via `hts_set_threads` if a CLI option is added (or reuse `num_threads`).
   - Convert `SAMMapping` to `bam1_t`:
     - `core.tid/pos/mtid/mpos` (0-based), `core.flag`, `core.qual`, `core.isize`.
     - Encode CIGAR, sequence, and qualities; append `NM:i`, `MD:Z`, `CB:Z` if present.
   - Write records with `sam_write1` and cleanly close `samFile` + `bam_hdr_t`.

3) **Indexing (optional)**
   - If `--write-index` is set, call `sam_idx_init` + `sam_idx_save` after writing.
   - Ensure output is coordinate-sorted before indexing; warn or refuse when low-mem overflow path does not guarantee global ordering.

4) **Docs and build updates**
   - Update `Makefile` to link with `-lhts` (and include paths if needed).
   - Update `README.md` and `chromap.1` for new flags and BAM/CRAM usage examples.
   - Note CRAM requires a reference FASTA (and `.fai`).

## Testing and validation
- Build:
  - `make clean && make`
- Quick BAM test:
  - `./chromap -x $INDEX -r $REF_FASTA -1 <R1> -2 <R2> --BAM -o /tmp/out.bam`
  - `samtools quickcheck -v /tmp/out.bam`
  - `samtools view -h /tmp/out.bam | head -n 5`
- SAM parity check:
  - `./chromap -x $INDEX -r $REF_FASTA -1 <R1> -2 <R2> --SAM -o /tmp/out.sam`
  - `samtools view -h /tmp/out.bam | diff -u /tmp/out.sam -`
- CRAM test (requires reference):
  - `./chromap -x $INDEX -r $REF_FASTA -1 <R1> -2 <R2> --CRAM -o /tmp/out.cram`
  - `samtools view -T $REF_FASTA /tmp/out.cram | head -n 5`
- Indexing test:
  - `./chromap ... --BAM --write-index -o /tmp/out.bam`
  - `samtools idxstats /tmp/out.bam`

## Troubleshooting
- If htslib is missing, the build will fail at link time; install `libhts-dev` or vendor htslib.
- For CRAM, missing `.fai` will cause htslib to fail; run `samtools faidx $REF_FASTA`.
- If indexing fails, confirm output is coordinate-sorted and not from an unordered overflow merge path.

## Rollback
- Revert to SAM output via `--SAM` or remove BAM/CRAM flags and rebuild without htslib.
