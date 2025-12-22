# Clarifying Questions for noY BAM Feature

Before creating a detailed implementation plan from the [noy_runbook.md](noy_runbook.md), please answer the following questions:

---

## 1. Feature Scope

Which version of the feature should be implemented?

- [ ] **Basic**: Just `--emit-noY-bam` (two outputs: primary/all + noY)
- [ ] **Three-stream**: all + Y-only + noY (as described at the end of the runbook)

---

## 2. Output Format

The flag is named `--emit-noY-bam` but chromap currently outputs SAM (text format). What format should the secondary output(s) be?

- [ ] **SAM only** - Same text format as current primary output (`.noY.sam`)
- [ ] **Binary BAM** - Requires piping through samtools or adding BAM encoding
- [ ] **Match primary** - Whatever format the primary output uses

---

## 3. Testing Scope

Should the implementation plan include the testing strategy from the runbook?

- [ ] **Yes** - Include unit tests, integration tests, and end-to-end tests in the plan
- [ ] **No** - Focus on implementation only; tests will be a separate effort

---

## 4. Path Derivation for Secondary Output

The runbook mentions deriving a default path like `<output>.noY.sam`. Should users be able to specify an explicit path?

- [ ] **Auto-derive only** - Always derive from primary output path (e.g., `out.sam` → `out.noY.sam`)
- [ ] **Optional explicit** - Add an optional `--noY-output` flag for explicit path, fall back to auto-derive
- [ ] **Required explicit** - User must specify the secondary output path explicitly

---

## 5. Filtering Granularity

For paired-end reads where one mate maps to Y and the other doesn't:

- [ ] **Filter the pair** - If either read in a pair touches Y, filter the entire pair from noY output
- [ ] **Filter individually** - Only filter the specific read that maps to Y (may result in orphan reads in noY output)

---

## 6. Secondary Alignments

How should secondary alignments be handled?

- [ ] **Read-ID based** - If a read has ANY alignment (primary or secondary) to Y, filter ALL its alignments from noY output
- [ ] **Per-alignment** - Only filter the specific alignment records that map to Y

*(The runbook suggests read-ID based filtering)*

---

## Additional Context from Codebase Analysis

Key files that will be modified:
- `src/chromap_driver.cc` - CLI flag parsing
- `src/mapping_parameters.h` - Store new parameters
- `src/mapping_writer.h/.cc` - Dual output handling
- `src/mapping_generator.h` - Y-hit detection during mapping
- `src/sequence_batch.h` - Reference contig name access

The `SAMMapping` class has `rid_` (reference ID) and `mrid_` (mate reference ID) fields that can be used to detect Y chromosome alignments.

