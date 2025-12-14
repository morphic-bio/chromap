#!/usr/bin/env python3
"""
Generate deterministic synthetic test data for chromap Y/NoY regression tests.

Creates:
- Reference FASTA with chr1, chr2, and chrY (distinct sequences)
- Single-end FASTQ with reads mapping to each contig + unmapped
- Paired-end FASTQ with various mapping scenarios

Uses fixed seed for reproducibility.
"""

import argparse
import random
import os

# Fixed seed for deterministic output
SEED = 42

# Contig lengths
CHR1_LEN = 1000
CHR2_LEN = 1000
CHRY_LEN = 1000

# Read parameters
READ_LEN = 100
QUAL_CHAR = 'I'  # Phred 40

def generate_random_sequence(length, seed_offset=0):
    """Generate a random DNA sequence with given length."""
    random.seed(SEED + seed_offset)
    return ''.join(random.choices('ACGT', k=length))

def reverse_complement(seq):
    """Return reverse complement of a DNA sequence."""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return ''.join(complement.get(base, 'N') for base in reversed(seq))

def generate_reference(output_dir):
    """Generate reference FASTA with unique sequences per contig."""
    ref_path = os.path.join(output_dir, 'test_ref.fa')
    
    # Generate distinct sequences for each contig using different seed offsets
    chr1_seq = generate_random_sequence(CHR1_LEN, seed_offset=1000)
    chr2_seq = generate_random_sequence(CHR2_LEN, seed_offset=2000)
    chrY_seq = generate_random_sequence(CHRY_LEN, seed_offset=3000)
    
    with open(ref_path, 'w') as f:
        f.write(f">chr1\n{chr1_seq}\n")
        f.write(f">chr2\n{chr2_seq}\n")
        f.write(f">chrY\n{chrY_seq}\n")
    
    print(f"Generated reference: {ref_path}")
    return {'chr1': chr1_seq, 'chr2': chr2_seq, 'chrY': chrY_seq}

def extract_read(contig_seq, position, strand='+'):
    """Extract a read from contig sequence at given position."""
    read_seq = contig_seq[position:position + READ_LEN]
    if strand == '-':
        read_seq = reverse_complement(read_seq)
    return read_seq

def write_fastq_record(f, read_id, seq):
    """Write a single FASTQ record."""
    qual = QUAL_CHAR * len(seq)
    f.write(f"@{read_id}\n{seq}\n+\n{qual}\n")

def generate_single_end_reads(output_dir, ref_seqs):
    """
    Generate single-end FASTQ with:
    - 20 reads mapping to chr1
    - 10 reads mapping to chr2
    - 15 reads mapping to chrY
    - 5 unmapped reads (random sequence)
    """
    se_path = os.path.join(output_dir, 'test_se.fq')
    random.seed(SEED + 100)
    
    read_num = 0
    with open(se_path, 'w') as f:
        # 20 reads to chr1
        for i in range(20):
            pos = random.randint(0, CHR1_LEN - READ_LEN - 1)
            strand = random.choice(['+', '-'])
            seq = extract_read(ref_seqs['chr1'], pos, strand)
            write_fastq_record(f, f"se_chr1_{read_num}", seq)
            read_num += 1
        
        # 10 reads to chr2
        for i in range(10):
            pos = random.randint(0, CHR2_LEN - READ_LEN - 1)
            strand = random.choice(['+', '-'])
            seq = extract_read(ref_seqs['chr2'], pos, strand)
            write_fastq_record(f, f"se_chr2_{read_num}", seq)
            read_num += 1
        
        # 15 reads to chrY
        for i in range(15):
            pos = random.randint(0, CHRY_LEN - READ_LEN - 1)
            strand = random.choice(['+', '-'])
            seq = extract_read(ref_seqs['chrY'], pos, strand)
            write_fastq_record(f, f"se_chrY_{read_num}", seq)
            read_num += 1
        
        # 5 unmapped reads (random sequence unlikely to map)
        for i in range(5):
            # Generate completely random sequence
            seq = ''.join(random.choices('ACGT', k=READ_LEN))
            write_fastq_record(f, f"se_unmapped_{read_num}", seq)
            read_num += 1
    
    print(f"Generated single-end reads: {se_path} ({read_num} reads)")
    return se_path

def generate_paired_end_reads(output_dir, ref_seqs):
    """
    Generate paired-end FASTQ with:
    - 15 pairs both mates to chr1
    - 10 pairs both mates to chr2
    - 10 pairs both mates to chrY
    - 5 pairs with one mate chr1, one mate chrY (split pairs)
    - 3 unmapped pairs
    """
    r1_path = os.path.join(output_dir, 'test_pe_R1.fq')
    r2_path = os.path.join(output_dir, 'test_pe_R2.fq')
    random.seed(SEED + 200)
    
    # Insert size range for paired-end
    INSERT_MIN = 150
    INSERT_MAX = 300
    
    pair_num = 0
    
    with open(r1_path, 'w') as f1, open(r2_path, 'w') as f2:
        # 15 pairs to chr1
        for i in range(15):
            insert_size = random.randint(INSERT_MIN, INSERT_MAX)
            pos1 = random.randint(0, CHR1_LEN - insert_size - 1)
            pos2 = pos1 + insert_size - READ_LEN
            
            seq1 = extract_read(ref_seqs['chr1'], pos1, '+')
            seq2 = extract_read(ref_seqs['chr1'], pos2, '-')
            
            write_fastq_record(f1, f"pe_chr1_{pair_num}/1", seq1)
            write_fastq_record(f2, f"pe_chr1_{pair_num}/2", seq2)
            pair_num += 1
        
        # 10 pairs to chr2
        for i in range(10):
            insert_size = random.randint(INSERT_MIN, INSERT_MAX)
            pos1 = random.randint(0, CHR2_LEN - insert_size - 1)
            pos2 = pos1 + insert_size - READ_LEN
            
            seq1 = extract_read(ref_seqs['chr2'], pos1, '+')
            seq2 = extract_read(ref_seqs['chr2'], pos2, '-')
            
            write_fastq_record(f1, f"pe_chr2_{pair_num}/1", seq1)
            write_fastq_record(f2, f"pe_chr2_{pair_num}/2", seq2)
            pair_num += 1
        
        # 10 pairs to chrY
        for i in range(10):
            insert_size = random.randint(INSERT_MIN, INSERT_MAX)
            pos1 = random.randint(0, CHRY_LEN - insert_size - 1)
            pos2 = pos1 + insert_size - READ_LEN
            
            seq1 = extract_read(ref_seqs['chrY'], pos1, '+')
            seq2 = extract_read(ref_seqs['chrY'], pos2, '-')
            
            write_fastq_record(f1, f"pe_chrY_{pair_num}/1", seq1)
            write_fastq_record(f2, f"pe_chrY_{pair_num}/2", seq2)
            pair_num += 1
        
        # 5 split pairs: R1 to chr1, R2 to chrY
        for i in range(5):
            pos1 = random.randint(0, CHR1_LEN - READ_LEN - 1)
            pos2 = random.randint(0, CHRY_LEN - READ_LEN - 1)
            
            seq1 = extract_read(ref_seqs['chr1'], pos1, '+')
            seq2 = extract_read(ref_seqs['chrY'], pos2, '-')
            
            write_fastq_record(f1, f"pe_split_{pair_num}/1", seq1)
            write_fastq_record(f2, f"pe_split_{pair_num}/2", seq2)
            pair_num += 1
        
        # 3 unmapped pairs
        for i in range(3):
            seq1 = ''.join(random.choices('ACGT', k=READ_LEN))
            seq2 = ''.join(random.choices('ACGT', k=READ_LEN))
            
            write_fastq_record(f1, f"pe_unmapped_{pair_num}/1", seq1)
            write_fastq_record(f2, f"pe_unmapped_{pair_num}/2", seq2)
            pair_num += 1
    
    print(f"Generated paired-end reads: {r1_path}, {r2_path} ({pair_num} pairs)")
    return r1_path, r2_path

def main():
    parser = argparse.ArgumentParser(description='Generate synthetic test data for chromap regression tests')
    parser.add_argument('--output-dir', '-o', default='.', 
                        help='Output directory for generated files')
    args = parser.parse_args()
    
    output_dir = args.output_dir
    os.makedirs(output_dir, exist_ok=True)
    
    print(f"Generating test data in {output_dir}...")
    print(f"Using seed: {SEED}")
    print()
    
    # Generate reference
    ref_seqs = generate_reference(output_dir)
    
    # Generate single-end reads
    generate_single_end_reads(output_dir, ref_seqs)
    
    # Generate paired-end reads
    generate_paired_end_reads(output_dir, ref_seqs)
    
    print()
    print("Test data generation complete!")
    print()
    print("Expected mapping counts:")
    print("  Single-end:")
    print("    chr1: 20, chr2: 10, chrY: 15, unmapped: 5")
    print("  Paired-end:")
    print("    chr1: 15 pairs, chr2: 10 pairs, chrY: 10 pairs")
    print("    split (chr1+chrY): 5 pairs, unmapped: 3 pairs")

if __name__ == '__main__':
    main()

