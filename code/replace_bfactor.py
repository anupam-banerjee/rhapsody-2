import numpy as np
import pandas as pd
import sys
import os

# Load predictions
predictions = pd.read_csv('predictions.txt', delimiter='\t', header=None, names=['file_name', 'residue_index', 'WT_residue', 'mutant_residue', 'probability'])

# Compute average probability for each residue
average_probabilities = predictions.groupby('residue_index')['probability'].mean().reset_index()

# Load PDB file
def load_pdb(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    return lines

# Update B-factor in PDB file
def update_bfactor(pdb_lines, average_probabilities):
    residue_prob_map = {int(row['residue_index']): row['probability'] for _, row in average_probabilities.iterrows()}
    updated_lines = []
    for line in pdb_lines:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            res_num = int(line[22:26].strip())
            if res_num in residue_prob_map:
                b_factor = f"{residue_prob_map[res_num]:.2f}".rjust(6)
                updated_line = f"{line[:60]}{b_factor}{line[66:]}"
                updated_lines.append(updated_line)
            else:
                updated_lines.append(line)
        else:
            updated_lines.append(line)
    return updated_lines

# Save updated PDB file
def save_pdb(file_path, pdb_lines):
    with open(file_path, 'w') as file:
        file.writelines(pdb_lines)

# Main function to handle file input and output
def main(input_pdb):
    # Ensure the file exists
    if not os.path.isfile(input_pdb):
        print(f"Error: File {input_pdb} not found.")
        return

    # Load PDB file
    pdb_lines = load_pdb(input_pdb)

    # Update B-factors in PDB file
    updated_pdb_lines = update_bfactor(pdb_lines, average_probabilities)

    # Create output file name
    base_name, _ = os.path.splitext(input_pdb)
    output_pdb = f"{base_name}_bfac.pdb"

    # Save updated PDB file
    save_pdb(output_pdb, updated_pdb_lines)

    print(f"Updated PDB file saved as {output_pdb}")

# Check if the script is run directly and handle command-line arguments
if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python replace_bfactor.py <input_pdb_file>")
    else:
        input_pdb_file = sys.argv[1]
        main(input_pdb_file)

