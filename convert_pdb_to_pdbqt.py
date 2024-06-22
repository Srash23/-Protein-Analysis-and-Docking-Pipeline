#!/usr/bin/env python3
import os
import sys
import subprocess
from rdkit import Chem
from rdkit.Chem import AllChem

def is_obabel_installed():
    """Check if Open Babel is correctly installed and accessible."""
    try:
        result = subprocess.run(['obabel', '-V'], capture_output=True, text=True)
        if result.returncode == 0:
            print("Open Babel Version:", result.stdout.strip())
            return True
        else:
            print("Open Babel not installed or not found in PATH.")
    except Exception as e:
        print("Failed to execute Open Babel:", str(e))
    return False

def prepare_receptor(filepath):
    if not is_obabel_installed():
        print("Please ensure that Open Babel is installed and accessible.")
        return
    
    try:
        mol = Chem.MolFromPDBFile(filepath, removeHs=False)
        if mol is None:
            raise ValueError("Could not parse the molecule from the file.")
        mol = Chem.AddHs(mol, explicitOnly=False, addCoords=True)
        AllChem.ComputeGasteigerCharges(mol)
        temp_filepath = "temp_receptor.pdb"
        Chem.MolToPDBFile(mol, temp_filepath)
        output_path = filepath.replace('.pdb', '.pdbqt')
        # Call obabel using subprocess to better handle errors
        result = subprocess.run(['obabel', temp_filepath, '-xr', '-opdbqt', '-O', output_path], capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"obabel failed: {result.stderr}")
        os.remove(temp_filepath)
        print(f"Receptor prepared and written to {output_path}")
    except Exception as e:
        print(f"Failed to prepare receptor due to: {str(e)}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python convert_pdb_to_pdbqt.py <input_pdb> <output_pdbqt>")
        sys.exit(1)
    input_pdb, output_pdbqt = sys.argv[1:3]
    prepare_receptor(input_pdb)
