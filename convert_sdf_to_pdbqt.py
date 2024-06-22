#!/usr/bin/env python
import os
import sys
import subprocess

def convert_sdf_to_pdbqt(input_sdf, output_pdbqt):
    command = ["obabel", "-isdf", input_sdf, "-opdbqt", "-O", output_pdbqt, "-h"]
    try:
        result = subprocess.run(command, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"Error converting SDF to PDBQT: {result.stderr}")
        else:
            print(f"Conversion successful: {output_pdbqt}")
    except Exception as e:
        print(f"Exception during conversion: {str(e)}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python convert_sdf_to_pdbqt.py <input_sdf> <output_pdbqt>")
        sys.exit(1)
    input_sdf, output_pdbqt = sys.argv[1:3]
    convert_sdf_to_pdbqt(input_sdf, output_pdbqt)
