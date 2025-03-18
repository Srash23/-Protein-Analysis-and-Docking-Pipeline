# Protein Analysis and Docking Pipeline

## Introduction

The Protein Analysis and Docking Pipeline is an integrated suite of bioinformatics tools designed to simplify and enhance protein structural and functional analysis. It provides an end-to-end computational framework for analyzing protein sequences, predicting structures, visualizing ligand-binding pockets, and performing docking studies. This pipeline supports a wide range of molecular biology and biochemical research applications.

By leveraging tools such as ESMfold, ProtParam, Clustal-Omega, fPockets, and AutoDock Vina, the pipeline enables researchers to:
1. Predict the 3D structure of proteins from amino acid sequences.
2. Compute physicochemical properties of proteins.
3. Perform multiple sequence alignment to identify conserved regions.
4. Identify potential ligand-binding sites in protein structures.
5. Prepare input files for docking simulations.
6. Run small-molecule blind docking between a protein and a ligand of interest.

## Installation and Usage (For End-Users)

To use this pipeline, ensure that the following dependencies are installed:

### Required Packages

**ESMfold –** Predicts 3D structures of proteins. [Installation Guide](https://github.com/facebookresearch/esm)

**ProtParam –** Computes physicochemical properties of proteins. [Documentation](https://biopython.org/docs/1.76/api/Bio.SeqUtils.ProtParam.html)

**Clustal-Omega –** Performs multiple sequence alignments. [Installation Guide](https://github.com/facebookresearch/esm)

**fPockets –** Identifies ligand-binding pockets in protein structures.

**AutoDock Vina –** Performs small-molecule docking. [Installation Guide](https://github.com/facebookresearch/esm)

### Python Dependencies:
1. pip install termcolor biopython
2. conda install -c conda-forge rdkit

### Other Required Installations
1. Open Babel: Used for molecular file conversions. Installation Guide

### Required Input Files
**1. Protein Sequence Files (.fasta) –** Used for structure prediction, alignment, and parametrization.

**2. Protein Structure Files (.pdb) –** Input for binding pocket identification and docking preparation.

**3. Ligand Structure Files (.sdf) –** Converted to .pdbqt format for docking.
   
_'''These files can be downloaded from:
UniProt Database (https://www.uniprot.org/),

NCBI Protein Database (https://www.ncbi.nlm.nih.gov/protein/),

PubChem (https://pubchem.ncbi.nlm.nih.gov/)'''_

## Pipeline Workflow
<img width="1497" alt="Screenshot 2025-03-18 at 5 13 35 PM" src="https://github.com/user-attachments/assets/1d6a8ebc-8147-470a-a284-ddfd89279406" />

### Step 1: Running ESMfold for 3D Structure Prediction
conda activate esmfold

python esm_fold.py "/path/to/torch_home" "/path/to/protein.fasta"

_**Input:** FASTA sequence

**Output:** .pdb file of the predicted 3D protein structure

### Step 2: Structural Analysis Using PyMOL
Open pymol_.py and modify the script to include the path to the generated .pdb file.

Run the script using:

pymol

run /path/to/pymol_.py

_**Output:** RMSD scores, B-factor analysis, electrostatic surface visualization._

### Step 3: Running the Main Pipeline
To run the core pipeline, execute:

python /path/to/main.py

_You will be prompted to enter:_

_1.Path to the docking environment_

_2.Path to Python analysis environment_

_3.Directory for scripts and results_

_4.Paths to .fasta and .pdb files_

_5.Paths for ligand file preparation_

### Step 4: Running fPockets for Binding Site Identification

python /path/to/run_fpocket.py /path/to/fpocket /path/to/protein.pdb /output/directory

_**Output:** Binding pocket detection results._

### Step 5: Preparing Files for Docking
1. Convert .pdb to .pdbqt

python /path/to/convert_pdb_to_pdbqt.py /path/to/protein.pdb /output/path/protein.pdbqt

2. Convert .sdf to .pdbqt

python /path/to/convert_sdf_to_pdbqt.py /path/to/ligand.sdf /output/path/ligand.pdbqt

_**Output:** Ligand and protein files formatted for docking.
_

### Step 6: Running AutoDock Vina for Docking

python /path/to/Docking_.py /path/to/receptor.pdbqt /path/to/ligand.pdbqt /path/to/config.txt /output/path/docking_results.pdbqt

_**Output:** Predicted binding conformations and docking scores._

## Key Insights
1. ESMfold successfully predicts 3D protein structures, allowing researchers to analyze molecular interactions in-depth.
2. ProtParam identifies important physicochemical properties such as molecular weight, isoelectric point, and instability index.
3. Clustal-Omega helps align multiple protein sequences, uncovering conserved functional domains.
4. fPockets efficiently detects potential ligand-binding sites, providing crucial data for drug discovery.
5. AutoDock Vina docking simulations predict ligand-receptor interactions, aiding rational drug design.

## Conclusion

The Protein Analysis and Docking Pipeline offers a powerful suite of bioinformatics tools to predict, analyze, and model protein structures and interactions. It integrates cutting-edge computational techniques to streamline protein-ligand docking studies, ultimately advancing research in structural biology and drug discovery.

