# Goal
Emit an additional SAM/BAM stream that excludes any read with primary or secondary alignments to Y contigs, gated by `--emit-noY-bam`, while keeping the normal output unchanged.

# Plan
- CLI plumbing: add `--emit-noY-bam` flag in `chromap_driver.cc`, guard so it only works with `--SAM`, derive a default secondary output path (e.g., `<output>.noY.sam/bam`), and stash it in `MappingParameters`.
- Y contig detection: once per run, scan reference contig names (case-insensitive, strip optional `chr` prefix) and mark RIDs that represent Y. Warn if none found.
- Collect Y-mapped reads during mapping:
  - Pass the Y RID mask and thread-local vectors into `MappingGenerator<SAMMapping>`.
  - When emitting SAM records (single-end or paired), record read IDs whose rid/mrid hit Y; include secondary alignments.
  - Merge per-thread vectors into a `std::unordered_set<uint32_t>` of read IDs to drop.
- Dual output in writer:
  - Extend `MappingWriter<SAMMapping>` to optionally open a second output file.
  - Mirror `@SQ` headers to both outputs.
  - When writing records, always emit to the primary output; emit to the no-Y output only if the read ID is not in the filtered set.
- Edge cases to watch:
  - `/dev/stdout`/`stderr` main outputs need a safe derived path for the no-Y file.
  - Low-memory spill path should continue to work; the filter is read-ID based, so it can be applied after merging thread-local Y hit sets.
  - Preserve existing summary/metadata behavior; filter should not change counters.

# Notes / what I learned
- SAM emission is centralized in `MappingWriter<SAMMapping>::AppendMapping`, and headers in `MappingWriter<SAMMapping>::OutputHeader`.
- Read IDs are stable across alignments (`SequenceBatch` assigns incremental ids), so they are reliable for post-hoc filtering.
- Mapping generators are the right place to notice Y hits because they see both primary and secondary alignments before dedup/output.
- Avoid destructive commands when reverting; restoring files via `git show HEAD:<path>` is a safe pattern.

# Tests (staged, minimal where possible)
- Unit-ish writer harness:
  - Build a tiny harness around `MappingWriter<SAMMapping>` with a mocked `SequenceBatch` containing 2–3 contigs (e.g., `chr1`, `chrY`, `chrM`).
  - Feed handcrafted `SAMMapping` records for:
    - read on `chr1` -> present in both outputs.
    - read on `chrY` -> absent from noY output.
    - read with primary on `chr1` and mate on `chrY` (mrid) -> filtered.
    - secondary alignment on Y -> filtered by read id.
  - Assert headers match on both outputs and records are byte-identical where present.
- Y mask detection:
  - Construct `SequenceBatch` names like `chrY`, `Y`, `chr1`, `GL000194.1`.
  - Verify only intended RIDs are flagged; also test a no-Y case to ensure warning path triggers.
- Mapping generator hook:
  - Drive `MappingGenerator<SAMMapping>` with synthetic `MappingInMemory` and a Y mask.
  - Confirm per-thread Y-hit vectors collect read ids for primary and secondary, single and paired.
- End-to-end smoke:
  - Use a 2-contig reference (`chr1`, `chrY`) and a handful of synthetic reads mapping to each.
  - Run in SAM mode with `--emit-noY-bam`; diff the two outputs to ensure only Y reads differ.
  - Run `samtools view` on both outputs to confirm parseable SAM.
- Low-mem/threads:
  - Enable `--low-mem` with a tiny threshold and threads >1 to ensure spill/merge preserves the filter behavior.
- `/dev/stdout` edge:
  - Run with primary output to stdout; verify derived `.noY.sam` path is created and populated.

# Three-stream variant (all / Y-only / noY)
- Add CLI flags for two extra outputs (e.g., `--emit-Y-bam` plus explicit paths) and mirror headers to all three streams.
- Filtering: reuse the Y rid mask; build `reads_with_Y` as the set to include for Y-only; emit:
  - all: everything,
  - Y-only: read ids in `reads_with_Y`,
  - noY: read ids not in `reads_with_Y`.
- Testing: extend the unit harness to assert the Y-only stream contains only Y hits; end-to-end, check `all = Y ∪ noY`.
- Pipes readiness: writing to 3 SAM files is identical to writing to 3 FIFOs; later, swap file paths for FIFOs and start `samtools view` readers first.
