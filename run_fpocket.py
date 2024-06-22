#!/usr/bin/env python3
import os
import sys
import subprocess

def run_fpocket(pdb_files, output_directory):
    fpocket_path = '/Users/srashtibajpai/Desktop/fpocket/bin/fpocket'
    if not os.path.exists(fpocket_path):
        raise FileNotFoundError("fpocket executable not found at " + fpocket_path)

    results = {}

    for pdb_file in pdb_files:
        output_dir = os.path.join(output_directory, os.path.basename(pdb_file) + "_out")
        if os.path.exists(output_dir):
            print(f"Output directory {output_dir} already exists, skipping {pdb_file}.")
            continue

        command = [fpocket_path, '-f', pdb_file, '-o', output_dir]
        print("Running command:", ' '.join(command))
        try:
            result = subprocess.run(command, capture_output=True, text=True, check=True)
            results[pdb_file] = output_dir
            print(f"fpocket analysis for {pdb_file} completed successfully.")
        except subprocess.CalledProcessError as e:
            print(f"Error running fpocket on {pdb_file}: {e.stderr}")

    return results

def main():
    if len(sys.argv) < 3:
        print("Usage: python run_fpocket.py <output_directory> <pdb_file1> <pdb_file2> ...")
        sys.exit(1)

    output_directory = sys.argv[1]
    pdb_files = sys.argv[2:]
    results = run_fpocket(pdb_files, output_directory)

    if results:
        for pdb_file, dir in results.items():
            print(f"Results for {pdb_file} are saved in {dir}")

if __name__ == "__main__":
    main()