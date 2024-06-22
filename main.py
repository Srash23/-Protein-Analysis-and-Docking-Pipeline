import os
import subprocess
from termcolor import colored

def run_command(command, env_path=None, print_output=True):
    print(colored(f"Starting: {' '.join(command)}", 'grey'))
    try:
        result = subprocess.run(command, capture_output=True, text=True, check=True)
        if print_output:
            print(colored(result.stdout, 'black'))
        if result.stderr:
            print(colored(result.stderr, 'red'))
    except subprocess.CalledProcessError as e:
        print(colored(f"Error in {' '.join(command)}: {e.stderr}", 'red'))
        return False
    return True

def main():
    print(colored("Welcome to the Protein Analysis and Docking Pipeline......", 'magenta'))
    
    docking_env = input(colored("Enter the path to the docking environment: ", 'blue'))
    analysis_env = input(colored("Enter the path for general analysis (Python environment): ", 'blue'))
    scripts_directory = input(colored("Enter the directory path where your scripts are located: ", 'blue'))
    result_directory = input(colored("Enter the directory path where you want to save the results: ", 'blue'))
    fpocket_path = input(colored("Enter the path to the fpocket executable: ", 'blue'))

    # Protein analysis
    print(colored("Starting protein analysis...", 'yellow'))
    fasta_files = input(colored("Enter the paths to the FASTA files for analysis separated by a space: ", 'blue')).split()
    analyze_script = os.path.join(scripts_directory, 'analyze_protein.py')
    output_csv = os.path.join(result_directory, "protein_analysis.csv")
    analyze_command = [os.path.join(analysis_env, 'bin', 'python'), analyze_script, output_csv] + fasta_files
    run_command(analyze_command)
    print(colored(f"Protein analysis results saved to {output_csv}.", 'green'))
    
    concatenated_fasta = os.path.join(result_directory, "concatenated.fasta")

    # Concatenation of FASTA files
    print(colored("Concatenating FASTA files...", 'yellow'))
    concatenate_command = [os.path.join(scripts_directory, 'concat_fasta.py'), concatenated_fasta] + fasta_files
    run_command(concatenate_command)
    print(colored("Concatenation complete", 'green'))
    
    aligned_fasta = os.path.join(result_directory, "aligned_sequences.fasta")

    # Multiple Sequence Alignment
    print(colored("Running multiple sequence alignment...", 'yellow'))
    python_executable = os.path.join(analysis_env, 'bin', 'python')  # Path to the Python executable in the analysis environment
    msa_script = os.path.join(scripts_directory, 'run_msa.py')
    msa_command = [python_executable, msa_script, concatenated_fasta, aligned_fasta]

    if not run_command(msa_command):
        print(colored("Multiple sequence alignment failed.", 'red'))
    else:
        print(colored("Multiple Sequence Alignment complete", 'green'))

    # fpocket processing without creating a specific "fpocket_results" directory
    print(colored("Running fpocket...", 'yellow'))
    pdb_files_fpocket = input(colored("Enter the paths to the PDB files for fpocket analysis separated by a space: ", 'blue')).split()
    for pdb_file in pdb_files_fpocket:
        output_dir = os.path.join(result_directory, os.path.basename(pdb_file).replace('.pdb', '_fpocket_out'))
        fpocket_command = ["python", os.path.join(scripts_directory, 'run_fpocket.py'), fpocket_path, pdb_file, output_dir]
        run_command(fpocket_command)
    print(colored("fpocket complete...", 'green'))

    # Convert PDB to PDBQT
    print(colored("Converting .pdb files to .pdbqt files...", 'yellow'))
    
    pdb_files = input(colored("Enter the paths to the PDB files to convert to PDBQT separated by a space: ", 'blue')).split()
    
    for pdb_file in pdb_files:
        pdbqt_file = pdb_file.replace('.pdb', '.pdbqt')
        # Construct the command as a flat list of strings
        convert_command = [os.path.join(docking_env, 'bin', 'python'), os.path.join(scripts_directory, 'convert_pdb_to_pdbqt.py'), pdb_file, pdbqt_file]
        run_command(convert_command)  # Note that docking_env is no longer passed here.

    print(colored("Converted .pdb files to .pdbqt files...", 'green'))
    

    
    # Convert SDF to PDBQT
    print(colored("Converting .sdf files to .pdbqt files...", 'yellow'))
    sdf_files = input(colored("Enter the paths to the SDF files to convert to PDBQT separated by a space: ", 'blue')).split()
    for sdf_file in sdf_files:
        pdbqt_file = sdf_file.replace('.sdf', '.pdbqt')
        # Ensure the command is a flat list of arguments
        convert_command = [os.path.join(scripts_directory, 'convert_sdf_to_pdbqt.py'), sdf_file, pdbqt_file]
        run_command(convert_command)
    print(colored("Converted .sdf files to .pdbqt files...", 'green'))


if __name__ == "__main__":
    main()