import sys
import csv
from Bio import SeqIO
from Bio.SeqUtils import ProtParam


def analyze_protein(sequence):
    protein_analysis = ProtParam.ProteinAnalysis(str(sequence))
    properties = {
        "Molecular Weight": protein_analysis.molecular_weight(),
        "Isoelectric Point": protein_analysis.isoelectric_point(),
        "Aromaticity": protein_analysis.aromaticity(),
        "Instability Index": protein_analysis.instability_index(),
        "Amino Acid Percentage": protein_analysis.get_amino_acids_percent(),
    }
    return properties


def process_files_and_output_csv(fasta_files, output_file):
    with open(output_file, 'w', newline='') as csvfile:
        fieldnames = ['Protein ID', 'Molecular Weight', 'Isoelectric Point', 'Aromaticity', 'Instability Index'] + [f"AA_{aa}" for aa in "ACDEFGHIKLMNPQRSTVWY"]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        
        for fasta_file in fasta_files:
            for record in SeqIO.parse(fasta_file, "fasta"):
                properties = analyze_protein(record.seq)
                row = {
                    'Protein ID': record.id,
                    'Molecular Weight': properties['Molecular Weight'],
                    'Isoelectric Point': properties['Isoelectric Point'],
                    'Aromaticity': properties['Aromaticity'],
                    'Instability Index': properties['Instability Index'],
                }
                row.update({f"AA_{aa}": properties['Amino Acid Percentage'].get(aa, 0) for aa in "ACDEFGHIKLMNPQRSTVWY"})
                writer.writerow(row)

def main():
    if len(sys.argv) < 3:
        print("Usage: python script_name.py output.csv input1.fasta input2.fasta ...")
        sys.exit(1)
    
    output_file = sys.argv[1]
    fasta_files = sys.argv[2:]
    process_files_and_output_csv(fasta_files, output_file)

if __name__ == "__main__":
    main()


