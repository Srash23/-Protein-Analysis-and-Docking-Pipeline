from pymol import cmd, stored

def load_and_align_proteins():
    # Hardcoded file paths for testing
    pdb_files = [
        "/Users/srashtibajpai/Desktop/isoform_1.pdb",
        "/Users/srashtibajpai/Desktop/isoform_3.pdb"
    ]

    # Assuming the first file is the reference
    reference = pdb_files[0].split('/')[-1][:-4]

    # Load proteins and align to reference
    for pdb_file in pdb_files:
        protein_name = pdb_file.split('/')[-1][:-4]
        cmd.load(pdb_file, protein_name)
        if protein_name != reference:
            alignment = cmd.align(protein_name, reference)
            with open(f"/Users/srashtibajpai/Desktop/{protein_name}_alignment_scores.txt", "w") as f:
                f.write(f"Alignment Score for {protein_name} vs {reference}: {alignment[0]}\n")

    # Calculate and save RMSD results
    for pdb_file in pdb_files[1:]:
        target_name = pdb_file.split('/')[-1][:-4]
        rmsd = cmd.align(target_name, reference)
        with open(f"/Users/srashtibajpai/Desktop/{target_name}_rmsd.txt", "w") as f:
            f.write(f"RMSD between {reference} and {target_name}: {rmsd[0]} Å\n")

    # Analyze B-factors
    analyze_b_factors(pdb_files)

    # Hydrophobicity and electrostatic surface analysis
    analyze_surface_properties(pdb_files)

def analyze_b_factors(pdb_files):
    for pdb_file in pdb_files:
        protein_name = pdb_file.split('/')[-1][:-4]
        stored.b_factors = []
        cmd.iterate(f"{protein_name} and name CA", "stored.b_factors.append(b)")
        average_b_factor = sum(stored.b_factors) / len(stored.b_factors) if stored.b_factors else 0
        with open(f"/Users/srashtibajpai/Desktop/{protein_name}_b_factors.txt", "w") as f:
            f.write(f"Average B-factor for {protein_name}: {average_b_factor:.2f}\n")

def analyze_surface_properties(pdb_files):
    for pdb_file in pdb_files:
        protein_name = pdb_file.split('/')[-1][:-4]
        # Hydrophobic surface
        cmd.show('surface', protein_name)
        cmd.set('surface_type', 1, protein_name)
        cmd.png(f"/Users/srashtibajpai/Desktop/{protein_name}_hydrophobic_surface.png")
        
        # Electrostatic surface
        cmd.set('surface_type', 2, protein_name)
        cmd.png(f"/Users/srashtibajpai/Desktop/{protein_name}_electrostatic_surface.png")

load_and_align_proteins()
