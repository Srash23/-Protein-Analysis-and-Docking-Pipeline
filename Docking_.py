import os
import sys

def run_vina(receptor_pdbqt, ligand_pdbqt, config_file, output_file):
    # Ensure filenames are properly quoted to handle special characters
    receptor_pdbqt = f"'{receptor_pdbqt}'"
    ligand_pdbqt = f"'{ligand_pdbqt}'"
    config_file = f"'{config_file}'"
    output_file = f"'{output_file}'"

    command = f"vina --receptor {receptor_pdbqt} --ligand {ligand_pdbqt} --config {config_file} --out {output_file}"
    os.system(command)

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python run_docking.py <receptor_pdbqt> <ligand_pdbqt> <config_file> <output_file>")
        sys.exit(1)

    receptor_pdbqt = sys.argv[1]
    ligand_pdbqt = sys.argv[2]
    config_file = sys.argv[3]
    output_file = sys.argv[4]
    run_vina(receptor_pdbqt, ligand_pdbqt, config_file, output_file)