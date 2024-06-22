#!/usr/bin/env python3
import sys

def concatenate_fasta_files(input_files, output_file):
    """Concatenate multiple .fasta files into one."""
    with open(output_file, 'w') as outfile:
        for file_name in input_files:
            with open(file_name, 'r') as infile:
                outfile.write(infile.read())

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python concat_fasta.py output.fasta input1.fasta input2.fasta ...")
        sys.exit(1)

    output_fasta = sys.argv[1]
    input_fastas = sys.argv[2:]
    concatenate_fasta_files(input_fastas, output_fasta)
    print(f"All files have been concatenated into {output_fasta}")
