import MDAnalysis as mda

def load_trajectory(dcd_file, prmtop_file):
    universe = mda.Universe(prmtop_file, dcd_file)
    return universe

def save_selection(trajectory, selection, output_prefix):
    for i, ts in enumerate(trajectory.trajectory[::10]):
        frame_number = i * 2
        selected_atoms = trajectory.select_atoms(selection)
        output_file = f"{output_prefix}_{frame_number}.pdb"
        selected_atoms.write(output_file)
        print(f"Saved selection for frame {frame_number} to {output_file}")

# Example usage:
dcd_file = "../wt-cov-complex-1000ns-rep.dcd"
prmtop_file = "../new-wt-cov-complex-solv.prmtop"
output_prefix = "selected_atoms"
selection = "protein or (resname ROH 0GB 3GB COV or (same resid as ((((type OW and (around 3.5 resid 495 and name OE1)) or (type OW and (around 3.5 resid 495 and name OE2))) and (type OW and (around 4.5 resid 289 and name C1))))))"  # Example selection, modify as needed

trajectory = load_trajectory(dcd_file, prmtop_file)
save_selection(trajectory, selection, output_prefix)

