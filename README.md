# Protein Structure Prediction and Docking Pipeline

This repository presents a modular bioinformatics pipeline designed to streamline protein analysis from sequence to docking. It combines deep learning-based structure prediction, sequence analysis, binding pocket identification, and small-molecule docking in a single, reproducible framework.

## Project Overview
This pipeline empowers researchers to:
- Predict 3D structures from FASTA sequences using ESMfold
- Compute physicochemical properties with ProtParam
- Align sequences using Clustal Omega
- Detect binding pockets via fPockets
- Convert and prepare molecular files for docking
- Perform protein-ligand docking using AutoDock Vina

## Required Tools & Dependencies

### Tools & Resources
| Tool | Purpose | Link |
|------|---------|------|
| **ESMfold** | 3D protein structure prediction | [GitHub](https://github.com/facebookresearch/esm) |
| **ProtParam (via Biopython)** | Compute protein properties | [Docs](https://biopython.org/docs/1.76/api/Bio.SeqUtils.ProtParam.html) |
| **Clustal Omega** | Multiple sequence alignment | [Clustal Omega](https://www.ebi.ac.uk/Tools/msa/clustalo/) |
| **fPockets** | Ligand-binding site prediction | - |
| **AutoDock Vina** | Protein-ligand docking | [Vina GitHub](https://github.com/ccsb-scripps/AutoDock-Vina) |
| **Open Babel** | Molecule format conversion | [Open Babel Docs](https://openbabel.org/docs/index.html) |

### Python Dependencies
```bash
pip install termcolor biopython
conda install -c conda-forge rdkit
```

## Required Input Files
| File Type | Format | Use |
|-----------|--------|-----|
| Protein sequence | `.fasta` | Structure prediction, alignment |
| Protein structure | `.pdb` | Pocket detection, docking input |
| Ligand structure  | `.sdf` | Docking input after conversion |

_Get sequences from [UniProt](https://www.uniprot.org/), [NCBI](https://www.ncbi.nlm.nih.gov/protein/), or [PubChem](https://pubchem.ncbi.nlm.nih.gov/) as needed._

## Workflow Overview
<img width="1019" alt="Screenshot 2025-04-21 at 4 57 16 PM" src="https://github.com/user-attachments/assets/aa91c53e-cf41-4bcc-a4b2-a3aac63caf76" />

## Step-by-Step Execution

### Structure Prediction (ESMfold)
```bash
conda activate esmfold
python esm_fold.py "/path/to/torch_home" "/path/to/protein.fasta"
```
**Input:** FASTA file → **Output:** PDB file

### Structural Visualization (PyMOL)
Edit `pymol_.py` to include your `.pdb` file path, then run:
```bash
pymol
run /path/to/pymol_.py
```
_Output:_ RMSD, B-factors, electrostatic visuals

### Main Analysis Pipeline
```bash
python /path/to/main.py
```
Prompts for paths to: docking env, analysis env, scripts, sequence, structure, ligand

### Binding Pocket Identification (fPockets)
```bash
python run_fpocket.py /path/to/fpocket /path/to/protein.pdb /output/folder
```

### File Conversion for Docking
```bash
python convert_pdb_to_pdbqt.py /path/to/protein.pdb /output/protein.pdbqt
python convert_sdf_to_pdbqt.py /path/to/ligand.sdf /output/ligand.pdbqt
```

### Run Docking (AutoDock Vina)
```bash
python Docking_.py /path/to/protein.pdbqt /path/to/ligand.pdbqt /path/to/config.txt /output/results.pdbqt
```

## Key Features & Insights
- **ESMfold** generates accurate 3D structures from sequence alone
- **ProtParam** provides essential molecular descriptors
- **Clustal Omega** highlights conserved domains
- **fPockets** detects ligand binding sites with spatial accuracy
- **AutoDock Vina** predicts binding affinity and conformations

## Applications
- Drug target validation
- Structure-function annotation
- Ligand docking for rational drug design
- Structure-based functional prediction

## Tech Stack
- **Languages**: Python, Shell
- **Libraries**: Biopython, RDKit, PyMOL
- **Tools**: ESMfold, AutoDock Vina, Open Babel

## License
MIT License – open to modification, reuse, and extension.
