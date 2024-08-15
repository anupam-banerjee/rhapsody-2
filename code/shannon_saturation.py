import sys
import os
import prody as pr
from collections import Counter
import math
from utilities.paths import Shannon_gap_folder, Shannon_nogap_folder, MSA_folder

# Amino acid codes and names
aa_codes = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
aa_names = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', 'GLY': 'G', 'HIS': 'H',
    'ILE': 'I', 'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q',
    'ARG': 'R', 'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
}

# Amino acid groups excluding gap '-'
amino_acid_groups_no_gaps1 = {
    'hydrophobic': {'A', 'V', 'I', 'L', 'M', 'F', 'W', 'P'},
    'polar': {'S', 'T', 'Y', 'N', 'Q', 'C', 'G'},
    'positive': {'K', 'R', 'H'},
    'negative': {'D', 'E'}
}

amino_acid_groups_no_gaps2 = {
    'smallest': {'G', 'A', 'S'},
    'small': {'T', 'C', 'P', 'D', 'N', 'V'},
    'medium': {'I', 'E', 'Q', 'L', 'H', 'M'},
    'large': {'F', 'Y', 'W', 'R', 'K'}
}

# Define amino acid groups including gap '-' (for gaps scenario)
amino_acid_groups_with_gaps1 = {
    'hydrophobic': {'A', 'V', 'I', 'L', 'M', 'F', 'W', 'P'},
    'polar': {'S', 'T', 'Y', 'N', 'Q', 'C', 'G'},
    'positive': {'K', 'R', 'H'},
    'negative': {'D', 'E'},
    'gap': {'-'}
}

amino_acid_groups_with_gaps2 = {
    'smallest': {'G', 'A', 'S'},
    'small': {'T', 'C', 'P', 'D', 'N', 'V'},
    'medium': {'I', 'E', 'Q', 'L', 'H', 'M'},
    'large': {'F', 'Y', 'W', 'R', 'K'},
    'gap': {'-'}
}

def extract_sequence_from_pdb(pdb_file):
    """Extract sequence from PDB file using ProDy."""
    structure = pr.parsePDB(pdb_file, subset='ca')
    resnames = structure.getResnames()
    sequence = ''.join([aa_names.get(resname, 'X') for resname in resnames])
    return sequence

def calculate_shannon_entropy_no_gaps(sequences, position, group1, group2):
    # Count the frequency of each amino acid at the specified position, ignoring gaps
    column = [seq[position] for seq in sequences if seq[position] != '-']
    freq_counts = Counter(column)

    # Calculate the total number of sequences excluding gaps
    total_seqs = len(column)

    # Avoid division by zero if column contains only gaps
    if total_seqs == 0:
        return 0.0, 0.0, 0.0

    # Calculate the Shannon entropy
    shannon_entropy = 0
    for count in freq_counts.values():
        probability = count / total_seqs
        shannon_entropy -= probability * math.log2(probability)

    # Calculate the group-based Shannon entropy for group1
    shannon_group1_entropy = 0
    for amino_acids in group1.values():
        group_count = sum(freq_counts[aa] for aa in amino_acids if aa in freq_counts)
        if group_count != 0:
            group_proba = group_count / total_seqs
            shannon_group1_entropy -= group_proba * math.log2(group_proba)

    # Calculate the group-based Shannon entropy for group2
    shannon_group2_entropy = 0
    for amino_acids in group2.values():
        group_count = sum(freq_counts[aa] for aa in amino_acids if aa in freq_counts)
        if group_count != 0:
            group_proba = group_count / total_seqs
            shannon_group2_entropy -= group_proba * math.log2(group_proba)

    return shannon_entropy, shannon_group1_entropy, shannon_group2_entropy

def calculate_shannon_entropy_with_gaps(sequences, position, group1, group2):
    """Calculate Shannon entropy at a specific position, including gaps."""
    column = [seq[position] for seq in sequences]
    freq_counts = Counter(column)
    total_seqs = len(sequences)

    shannon_entropy = -sum((count / total_seqs) * math.log2(count / total_seqs) for count in freq_counts.values())

    shannon_group1_entropy = -sum(
        (sum(freq_counts[aa] for aa in amino_acids) / total_seqs) * math.log2(
            sum(freq_counts[aa] for aa in amino_acids) / total_seqs)
        for amino_acids in group1.values() if sum(freq_counts[aa] for aa in amino_acids) > 0
    )

    shannon_group2_entropy = -sum(
        (sum(freq_counts[aa] for aa in amino_acids) / total_seqs) * math.log2(
            sum(freq_counts[aa] for aa in amino_acids) / total_seqs)
        for amino_acids in group2.values() if sum(freq_counts[aa] for aa in amino_acids) > 0
    )

    return shannon_entropy, shannon_group1_entropy, shannon_group2_entropy

def calculate_mean_stddev(values):
    mean = sum(values) / len(values)
    variance = sum((x - mean) ** 2 for x in values) / len(values)
    stddev = math.sqrt(variance)
    return mean, stddev

def calculate_z_score(value, mean, stddev):
    if stddev == 0:
        return 0
    return (value - mean) / stddev

def generate_shannon_files(msa_file, pdb_name):
    """Generate Shannon entropy files (with and without gaps) and save them."""
    with open(msa_file, 'r') as f:
        sequences = [line.strip() for line in f if line.strip()]

    output_file_gap = os.path.join(Shannon_gap_folder, f"{pdb_name}.shagap")
    output_file_nogap = os.path.join(Shannon_nogap_folder, f"{pdb_name}.shanogap")

    # Calculate entropy values for all positions
    all_entropies_with_gaps = [
        calculate_shannon_entropy_with_gaps(sequences, pos - 1, amino_acid_groups_with_gaps1, amino_acid_groups_with_gaps2)
        for pos in range(1, len(sequences[0]) + 1)
    ]
    all_entropies_no_gaps = [
        calculate_shannon_entropy_no_gaps(sequences, pos - 1, amino_acid_groups_no_gaps1, amino_acid_groups_no_gaps2)
        for pos in range(1, len(sequences[0]) + 1)
    ]

    overall_entropies_with_gaps = [e[0] for e in all_entropies_with_gaps]
    group1_entropies_with_gaps = [e[1] for e in all_entropies_with_gaps]
    group2_entropies_with_gaps = [e[2] for e in all_entropies_with_gaps]

    overall_mean_with_gaps, overall_stddev_with_gaps = calculate_mean_stddev(overall_entropies_with_gaps)
    group1_mean_with_gaps, group1_stddev_with_gaps = calculate_mean_stddev(group1_entropies_with_gaps)
    group2_mean_with_gaps, group2_stddev_with_gaps = calculate_mean_stddev(group2_entropies_with_gaps)

    overall_entropies_no_gaps = [e[0] for e in all_entropies_no_gaps]
    group1_entropies_no_gaps = [e[1] for e in all_entropies_no_gaps]
    group2_entropies_no_gaps = [e[2] for e in all_entropies_no_gaps]

    overall_mean_no_gaps, overall_stddev_no_gaps = calculate_mean_stddev(overall_entropies_no_gaps)
    group1_mean_no_gaps, group1_stddev_no_gaps = calculate_mean_stddev(group1_entropies_no_gaps)
    group2_mean_no_gaps, group2_stddev_no_gaps = calculate_mean_stddev(group2_entropies_no_gaps)

    # Add header to the files
    header = (
        "Position|Overall Entropy|Z-score Overall Entropy|"
        "Group1 Entropy|Z-score Group1 Entropy|"
        "Group2 Entropy|Z-score Group2 Entropy\n"
    )

    with open(output_file_gap, 'w') as f_gap, open(output_file_nogap, 'w') as f_nogap:
        # Write the headers
        f_gap.write(header)
        f_nogap.write(header)

        for position in range(len(sequences[0])):
            entropy_with_gaps, group1_with_gaps, group2_with_gaps = all_entropies_with_gaps[position]
            zscore_overall_with_gaps = calculate_z_score(entropy_with_gaps, overall_mean_with_gaps, overall_stddev_with_gaps)
            zscore_group1_with_gaps = calculate_z_score(group1_with_gaps, group1_mean_with_gaps, group1_stddev_with_gaps)
            zscore_group2_with_gaps = calculate_z_score(group2_with_gaps, group2_mean_with_gaps, group2_stddev_with_gaps)

            entropy_no_gaps, group1_no_gaps, group2_no_gaps = all_entropies_no_gaps[position]
            zscore_overall_no_gaps = calculate_z_score(entropy_no_gaps, overall_mean_no_gaps, overall_stddev_no_gaps)
            zscore_group1_no_gaps = calculate_z_score(group1_no_gaps, group1_mean_no_gaps, group1_stddev_no_gaps)
            zscore_group2_no_gaps = calculate_z_score(group2_no_gaps, group2_mean_no_gaps, group2_stddev_no_gaps)

            f_gap.write(f"{position + 1}|{entropy_with_gaps:.3f}|{zscore_overall_with_gaps:.3f}|{group1_with_gaps:.3f}|{zscore_group1_with_gaps:.3f}|{group2_with_gaps:.3f}|{zscore_group2_with_gaps:.3f}\n")
            f_nogap.write(f"{position + 1}|{entropy_no_gaps:.3f}|{zscore_overall_no_gaps:.3f}|{group1_no_gaps:.3f}|{zscore_group1_no_gaps:.3f}|{group2_no_gaps:.3f}|{zscore_group2_no_gaps:.3f}\n")

def read_shannon_entropy_file(filepath, residue_index):
    """Read a specific line from the Shannon entropy file corresponding to the residue index."""
    try:
        with open(filepath, 'r') as file:
            # Skip the header line
            header = file.readline().strip()
            #print(f"Skipping header: {header}")  # Debugging line to ensure header is skipped
            
            for i, line in enumerate(file):
                if i == residue_index - 1:  # Residue index is 1-based
                    #print(f"Reading line {i+2} from {filepath}: {line.strip()}")  # +2 because of header
                    return line.strip().split('|')
    except FileNotFoundError:
        print(f"File not found: {filepath}")
        return None
    return None


def generate_saturation_mutagenesis_shannon(pdb_file, sequence, msa_file, output_file):
    pdb_name = os.path.splitext(os.path.basename(pdb_file))[0]
    structure = pr.parsePDB(pdb_file, subset='ca')
    resnums = structure.getResnums()
    resnames = structure.getResnames()

    shagap_file = os.path.join(Shannon_gap_folder, f"{pdb_name}.shagap")
    shanogap_file = os.path.join(Shannon_nogap_folder, f"{pdb_name}.shanogap")

    # Check if the Shannon files exist; if not, generate them
    if not os.path.isfile(shagap_file) or not os.path.isfile(shanogap_file):
        print(f"Shannon files not found for {pdb_name}. Generating now...")
        generate_shannon_files(msa_file, pdb_name)

    with open(output_file, 'w') as outfile:
        # Write the headers to the TSV file
        headers = [
            "PDB", "ResidueNumber", "WTResidue", "MutResidue",
            "ShannonWithGaps", "ZScoreOverallWithGaps", "Group1WithGaps",
            "ZScoreGroup1WithGaps", "Group2WithGaps", "ZScoreGroup2WithGaps",
            "ShannonNoGaps", "ZScoreOverallNoGaps", "Group1NoGaps",
            "ZScoreGroup1NoGaps", "Group2NoGaps", "ZScoreGroup2NoGaps"
        ]
        outfile.write("\t".join(headers) + "\n")

        for i, (resnum, resname) in enumerate(zip(resnums, resnames)):
            print(f"Processing residue {i+1}/{len(resnums)}")
            wt_residue = aa_names.get(resname, 'X')
            if wt_residue == 'X':
                print(f"Unknown residue {resname} at position {resnum}")
                continue

            for mut_residue in aa_codes:
                # Read the corresponding lines from the Shannon entropy files
                shagap_entry = read_shannon_entropy_file(shagap_file, resnum)
                shanogap_entry = read_shannon_entropy_file(shanogap_file, resnum)

                if shagap_entry and shanogap_entry:
                    outfile.write(f"{pdb_name}\t{resnum}\t{wt_residue}\t{mut_residue}\t")

                    # Populate Shannon entropy with gaps
                    outfile.write(f"{shagap_entry[1]}\t{shagap_entry[2]}\t{shagap_entry[3]}\t{shagap_entry[4]}\t{shagap_entry[5]}\t{shagap_entry[6]}\t")

                    # Populate Shannon entropy without gaps
                    outfile.write(f"{shanogap_entry[1]}\t{shanogap_entry[2]}\t{shanogap_entry[3]}\t{shanogap_entry[4]}\t{shanogap_entry[5]}\t{shanogap_entry[6]}\n")

                else:
                    print(f"Error in mapping {pdb_name} at residue {resnum}")
                    print(f"shagap_entry: {shagap_entry}, shanogap_entry: {shanogap_entry}")
                    outfile.write(f"{pdb_name}\t{resnum}\t{wt_residue}\t{mut_residue}\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN\n")


def main():
    if len(sys.argv) != 3:
        print("Usage: python shannon_saturation_script.py -f <pdb_file>")
        return

    if sys.argv[1] != '-f':
        print("Invalid argument. Use -f to specify the PDB file.")
        return

    pdb_file = sys.argv[2]
    pdb_name = os.path.splitext(os.path.basename(pdb_file))[0]

    if not os.path.isfile(pdb_file):
        print(f"PDB file {pdb_file} does not exist.")
        return

    shagap_file = os.path.join(Shannon_gap_folder, f"{pdb_name}.shagap")
    shanogap_file = os.path.join(Shannon_nogap_folder, f"{pdb_name}.shanogap")

    if os.path.isfile(shagap_file) and os.path.isfile(shanogap_file):
        print(f"Precomputed Shannon files found for {pdb_name}.")
        msa_file = None
    else:
        msa_file = os.path.join(MSA_folder, f"{pdb_name}.seqmsa")
        if not os.path.isfile(msa_file):
            print(f"MSA file {msa_file} does not exist.")
            return

    output_file = f"saturation_mutagenesis_shannon.tsv"

    try:
        print("Extracting sequence from PDB file...")
        sequence = extract_sequence_from_pdb(pdb_file)

        generate_saturation_mutagenesis_shannon(pdb_file, sequence, msa_file, output_file)
        print(f"Shannon entropy features saved in {output_file}")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    main()

