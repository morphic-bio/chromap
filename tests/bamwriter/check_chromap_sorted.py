#!/usr/bin/env python3
"""Verify BAM is sorted by chromap's primary sort key: (tid, pos, flag, mtid, mpos, isize)
Note: read_id tie-break is not validated (it's internal to chromap and not stored in BAM tags).
This validates the main coordinate sort order which is sufficient for indexing compatibility."""
import pysam
import sys

def get_sort_key(read):
    """Extract chromap primary sort key from a pysam read (without read_id tie-break)."""
    # Normalize unmapped to sort last
    tid = read.reference_id if read.reference_id >= 0 else 0x7FFFFFFF
    pos = read.reference_start if read.reference_start >= 0 else 0x7FFFFFFF
    
    # Return primary key only (read_id tie-break is internal to chromap)
    return (tid, pos, read.flag, read.next_reference_id, read.next_reference_start, read.template_length)

def check_sorted(bam_path):
    """Check if BAM is sorted by chromap's sort key."""
    try:
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            prev_key = None
            for i, read in enumerate(bam):
                key = get_sort_key(read)
                if prev_key is not None and key < prev_key:
                    print(f"UNSORTED at record {i}: {prev_key} > {key}", file=sys.stderr)
                    print(f"  Previous: tid={prev_key[0]}, pos={prev_key[1]}, flag={prev_key[2]}, mtid={prev_key[3]}, mpos={prev_key[4]}, isize={prev_key[5]}", file=sys.stderr)
                    print(f"  Current:  tid={key[0]}, pos={key[1]}, flag={key[2]}, mtid={key[3]}, mpos={key[4]}, isize={key[5]}", file=sys.stderr)
                    return 1
                prev_key = key
        print("SORTED OK")
        return 0
    except Exception as e:
        print(f"Error checking sortedness: {e}", file=sys.stderr)
        return 1

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: check_chromap_sorted.py <bam_file>", file=sys.stderr)
        sys.exit(1)
    sys.exit(check_sorted(sys.argv[1]))

