import torch
import esm
import os
import argparse

def load_sequence_from_fasta(fasta_file):
    with open(fasta_file, 'r') as file:
        lines = file.readlines()
        sequence = ''.join([line.strip() for line in lines if not line.startswith('>')])
    return sequence

def main(torch_home, fasta_file):
    print("Script started...")

    # Set the environment variable for torch
    os.environ["TORCH_HOME"] = torch_home

    # Load the ESMFold model
    model = esm.pretrained.esmfold_v1()
    model = model.eval().cuda()
    print("Model loaded...")

    # Load sequence from fasta file
    sequence = load_sequence_from_fasta(fasta_file)

    # Infer the 3D structure
    with torch.no_grad():
        output = model.infer_pdb(sequence)

    # Write the result to a file
    with open("isoform_1.pdb", "w") as f:
        f.write(output)

    print("Predictions done...")