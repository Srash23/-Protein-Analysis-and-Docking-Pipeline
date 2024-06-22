#!/usr/bin/env python3

import os
import sys
import csv
from Bio import AlignIO
from Bio.Align import AlignInfo

def run_clustal_omega(concatenated_file, output_file="aligned_sequences.fasta"):
    command = f"clustalo -i {concatenated_file} -o {output_file} --auto -v --force"
    print("Running command:", command)
    result = os.system(command)
    if result != 0:
        raise RuntimeError("Clustal Omega failed to run properly")
    return output_file

def analyze_alignment(output_file):
    alignment = AlignIO.read(output_file, "fasta")
    summary_align = AlignInfo.SummaryInfo(alignment)
    consensus = summary_align.dumb_consensus()

    # Write results to a CSV file
    with open('alignment_results.csv', 'w', newline='') as csvfile:
        fieldnames = ['Consensus Sequence']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerow({'Consensus Sequence': str(consensus)})

    print("Consensus Sequence written to CSV.")

def main():
    if len(sys.argv) < 3:
        print("Usage: python run_msa.py <input_concatenated_fasta> <output_aligned_fasta>")
        sys.exit(1)

    input_concatenated_fasta = sys.argv[1]
    output_aligned_fasta = sys.argv[2]
    output_file = run_clustal_omega(input_concatenated_fasta, output_aligned_fasta)
    analyze_alignment(output_file)

if __name__ == "__main__":
    main()
