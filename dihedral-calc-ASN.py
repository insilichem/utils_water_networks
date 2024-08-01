import MDAnalysis as mda
import numpy as np
import os

"""
CHI1 AND CHI2 OF GLUTAMATE: https://www.researchgate.net/publication/44651362_Beyond_rotamers_A_generative_probabilistic_model_of_side_chains_in_proteins/figures?lo=1 

"""

# Define the folder containing PDB files
pdb_folder = "./"

# Create an empty list to store results
results = []

# Iterate over all PDB files in the folder
for pdb_file in os.listdir(pdb_folder):
    if pdb_file.endswith(".pdb"):
        pdb_path = os.path.join(pdb_folder, pdb_file)
        
        # Load the PDB file
        u = mda.Universe(pdb_path)
        
        # Select atoms for the φ and ψ dihedral angle calculations
        atom_selection_chi1 = "(resid 219 and name CA) or (resid 219 and name CB) or (resid 219 and name CG) or (resid 219 and name ND2)"
        atom_selection_chi2 = "(resid 219 and name CB) or (resid 219 and name CG) or (resid 219 and name ND2) or (resid 219 and name OD1)"
      
        
        # Get the selected atoms
        chi1_atoms = u.select_atoms(atom_selection_chi1)
        chi2_atoms = u.select_atoms(atom_selection_chi2)

        
        # Get positions of the selected atoms
        chi1_positions = chi1_atoms.positions
        chi2_positions = chi2_atoms.positions

        
        # Calculate the φ and ψ dihedral angles
        chi1_dihedral_angle = np.rad2deg(mda.lib.distances.calc_dihedrals(chi1_positions[:-3], chi1_positions[1:-2], chi1_positions[2:-1], chi1_positions[3:]))
        chi2_dihedral_angle = np.rad2deg(mda.lib.distances.calc_dihedrals(chi2_positions[:-3], chi2_positions[1:-2], chi2_positions[2:-1], chi2_positions[3:]))

        

        # Append the results to the list
        results.append((pdb_file, chi1_dihedral_angle[0], chi2_dihedral_angle[0]))

# Save the results to a text file
output_file = "dihedral_angles_results-ASN219.txt"
with open(output_file, "w") as f:
    #f.write("PDB_File\tchi1_Dihedral_Angle\tchi2_Dihedral_Angle\n")
    for result in results:
        f.write(f"{result[0]},{result[1]},{result[2]}\n")

print("Results saved to:", output_file)

