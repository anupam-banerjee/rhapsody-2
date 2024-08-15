import sys
import os
import math
import json
from collections import Counter
from Bio.Blast import NCBIWWW, NCBIXML
import prody as pr
from utilities.paths import PSIC_folder, MSA_folder, BLAST_folder

# Define amino acid groups including gap '-' as a subgroup
amino_acid_groups1 = {
    'hydrophobic': {'A', 'V', 'I', 'L', 'M', 'F', 'W', 'P'},
    'polar': {'S', 'T', 'Y', 'N', 'Q', 'C', 'G'},
    'positive': {'K', 'R', 'H'},
    'negative': {'D', 'E'},
    'other': {'B', 'J', 'Z', 'X'},  # Non-standard or ambiguous amino acids
    'gap': {'-'}  # Gap character
}

amino_acid_groups2 = {
    'smallest': {'G', 'A', 'S'},
    'small': {'T', 'C', 'P', 'D', 'N', 'V'},
    'medium': {'I', 'E', 'Q', 'L', 'H', 'M'},
    'large': {'F', 'Y', 'W', 'R', 'K'},
    'other': {'B', 'J', 'Z', 'X'},  # Non-standard or ambiguous amino acids
    'gap': {'-'}  # Gap character
}

aa_codes = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
aa_names = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', 'GLY': 'G', 'HIS': 'H',
    'ILE': 'I', 'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q',
    'ARG': 'R', 'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
}

def extract_sequence_from_pdb(pdb_file):
    """Extract sequence from PDB file using ProDy."""
    structure = pr.parsePDB(pdb_file, subset='ca')
    resnames = structure.getResnames()
    sequence = ''.join([aa_names.get(resname, 'X') for resname in resnames])
    return sequence

def validate_seqmsa_and_json_with_pdb(seqmsa_file, json_file, pdb_sequence):
    """Validate that the sequence lengths in position_freqs in JSON match the PDB sequence length, 
    and only if necessary look for the seqmsa file for further validation."""

    # First, attempt to validate using the JSON file
    print(f"Loading and processing {json_file}...")
    try:
        with open(json_file, 'r') as f:
            json_data = json.load(f)
    except FileNotFoundError:
        print(f"Error: JSON file {json_file} not found.")
        return False
    except Exception as e:
        print(f"Error processing JSON file {json_file}: {e}")
        return False
    
    position_freqs_length = len(json_data.get('position_freqs', []))
    print(f"Length of position_freqs in {json_file}: {position_freqs_length}")

    # Compare lengths with the PDB sequence length
    pdb_sequence_length = len(pdb_sequence)
    print(f"Length of PDB sequence: {pdb_sequence_length}")

    if position_freqs_length == pdb_sequence_length:
        print("Validation successful with JSON file: The lengths match with the PDB sequence.")
        return True

    print(f"Validation failed: Length of position_freqs ({position_freqs_length}) does not match PDB sequence length ({pdb_sequence_length}).")
    print(f"Proceeding to validate MSA file as fallback...")

    # If JSON validation fails, proceed to validate the MSA file
    if not os.path.isfile(seqmsa_file):
        print(f"Error: MSA file {seqmsa_file} not found.")
        return False

    print(f"Loading and processing {seqmsa_file}...")
    try:
        with open(seqmsa_file, 'r') as f:
            sequences = [line.strip() for line in f if line.strip()]
    except Exception as e:
        print(f"Error processing MSA file {seqmsa_file}: {e}")
        return False

    if not sequences:
        print(f"Error: No sequences found in {seqmsa_file}")
        return False

    seqmsa_length = len(sequences[0])  # All sequences in seqmsa should have the same length
    print(f"Length of sequences in {seqmsa_file} (with gaps): {seqmsa_length}")

    if seqmsa_length != pdb_sequence_length:
        print(f"Validation failed: Length of seqmsa ({seqmsa_length}) does not match PDB sequence length ({pdb_sequence_length}).")
        print(f"The MSA file {seqmsa_file} saved with the same name does not correspond to the present input.")
        return False

    print("Validation successful with MSA file: The lengths match with the PDB sequence.")
    return True


def validate_blast_xml_with_pdb(blast_file, pdb_sequence_length):
    """Validate that the BLAST XML file corresponds to the current PDB sequence."""
    print(f"Checking BLAST file {blast_file}...")

    # Check file size or length by reading the XML data
    with open(blast_file, 'r') as f:
        xml_data = f.read()

    if len(xml_data) == 0:
        print(f"Error: BLAST file {blast_file} is empty.")
        return False

    # Simple check by validating the length of the input sequence in the BLAST file (if included)
    if len(xml_data) < pdb_sequence_length:
        print(f"Validation failed: The length of the BLAST XML data is less than the PDB sequence length ({pdb_sequence_length}).")
        print(f"The BLAST XML file {blast_file} saved with the same name does not correspond to the present input.")
        return False

    print("Validation successful: The BLAST XML file appears to match the PDB sequence.")
    return True

def perform_blast_search(sequence, expect=0.001, hitlist_size=50000, blast_file=None):
    """Perform BLAST search for the given sequence and save the results."""
    print("Performing BLAST search...")

    # Warning message
    print("Warning: As Rhapsody-2 requires all hits below the e-value, a hitlist_size of 50,000 is considered.")
    print("This can incur significant runtimes, and the BLAST search might fail.")

    result_handle = NCBIWWW.qblast("blastp", "nr", sequence, expect=expect, hitlist_size=hitlist_size)
    if blast_file:
        with open(blast_file, "w") as out_handle:
            out_handle.write(result_handle.read())
    else:
        with open("blast_results.xml", "w") as out_handle:
            out_handle.write(result_handle.read())
    result_handle.close()

def parse_blast_msa(blast_file, seq, msa_out_path):
    """Parse BLAST results and write MSA file, even from a potentially corrupted XML file."""
    from Bio.Blast import NCBIXML  # for parsing XML
    import tempfile
    
    print("Parsing BLAST results and creating MSA file...")

    seq_length = len(seq)
    with open(msa_out_path, "w") as fp:
        try:
            with open(blast_file, "r") as blast_handle:
                blast_records = NCBIXML.parse(blast_handle)
                
                for blast_record in blast_records:
                    for ali in blast_record.alignments:
                        for hsp in ali.hsps:
                            query_cover = (hsp.query_end - hsp.query_start + 1) / seq_length * 100
                            if query_cover >= 30 and hsp.expect <= 0.001:
                                pre_seq = '-' * (hsp.query_start - 1)
                                pos_seq = '-' * max(0, (seq_length - hsp.query_end))
                                seq_match = ''.join(hsp.sbjct[count] for count, ele in enumerate(hsp.query) if ele != '-')
                                aligned_seq = pre_seq + seq_match + pos_seq
                                fp.write(aligned_seq + '\n')
        except Exception as e:
            print(f"An error occurred during parsing: {e}")
            print("Continuing with the valid entries parsed so far.")

    print(f"MSA file created at {msa_out_path}")



def precompute_frequencies(msa_file):
    with open(msa_file, 'r') as f:
        sequences = [line.strip() for line in f if line.strip()]

    if len(sequences) < 3:
        print(f"Insufficient records found in MSA file: {msa_file}")
        return None

    seq_len = len(sequences[0])
    position_freqs = [Counter() for _ in range(seq_len)]
    overall_freqs = Counter()

    for seq in sequences:
        for i, aa in enumerate(seq):
            position_freqs[i][aa] += 1
            overall_freqs[aa] += 1

    total_seqs = len(sequences)
    total_aa_counts = sum(overall_freqs.values())
    
    overall_freqs = {aa: count / total_aa_counts for aa, count in overall_freqs.items()}
    position_freqs = [{aa: count / total_seqs for aa, count in pos_freq.items()} for pos_freq in position_freqs]

    return overall_freqs, position_freqs

def save_precomputed_frequencies(output_folder, msa_filename, precomputed_data):
    overall_freqs, position_freqs = precomputed_data
    data = {
        'overall_freqs': overall_freqs,
        'position_freqs': position_freqs
    }
    output_path = os.path.join(output_folder, os.path.basename(msa_filename).replace('.pdb', '') + '.seqmsa.json')
    with open(output_path, 'w') as f:
        json.dump(data, f)

def calculate_psic(precomputed_folder, msa_filename, position, amino_acid, groups):
    precomputed_path = os.path.join(precomputed_folder, os.path.basename(msa_filename).replace('.pdb', '') + '.seqmsa.json')
    
    try:
        with open(precomputed_path, 'r') as f:
            data = json.load(f)
    except Exception as e:
        print(f"Error reading precomputed JSON file: {e}")
        return float('NaN'), float('NaN'), 'unknown'
    
    overall_freqs = data.get('overall_freqs', {})
    position_freqs = data.get('position_freqs', [])
    
    if not position_freqs or position < 1 or position > len(position_freqs):
        print("Position out of range or no position frequencies available")
        return float('NaN'), float('NaN'), 'unknown'
    
    position -= 1  # Adjust to 0-based index
    column_freqs = position_freqs[position]
    
    total_seqs = sum(column_freqs.values())
    observed_freq = column_freqs.get(amino_acid, 0) / total_seqs if total_seqs > 0 else 0
    
    total_aa_counts = sum(overall_freqs.values())
    expected_freq = overall_freqs.get(amino_acid, 0) / total_aa_counts if total_aa_counts > 0 else 0
    
    if observed_freq == 0 or expected_freq == 0:
        psic_score = float('-inf')
    else:
        psic_score = math.log(observed_freq / expected_freq)
    
    amino_acid_group = get_amino_acid_group(amino_acid, groups)
    
    if amino_acid_group == 'unknown':
        print(f"Unknown amino acid: {amino_acid}")
        return float('NaN'), float('NaN'), 'unknown'
    
    group_freq = sum(column_freqs.get(aa, 0) for aa in groups[amino_acid_group]) / total_seqs if total_seqs > 0 else 0
    group_expected_freq = sum(overall_freqs.get(aa, 0) for aa in groups[amino_acid_group]) / total_aa_counts if total_aa_counts > 0 else 0
    
    if group_freq == 0 or group_expected_freq == 0:
        group_psic_score = float('-inf')
    else:
        group_psic_score = math.log(group_freq / group_expected_freq)

    return psic_score, group_psic_score, amino_acid_group

def get_amino_acid_group(amino_acid, groups):
    for group, amino_acids in groups.items():
        if amino_acid in amino_acids:
            return group
    return 'unknown'

def generate_saturation_mutagenesis_features(pdb_file, sequence, precomputed_folder, output_file):
    print("Generating saturation mutagenesis features...")
    pdb_file_name = os.path.basename(pdb_file)
    structure = pr.parsePDB(pdb_file, subset='ca')
    resnums = structure.getResnums()
    resnames = structure.getResnames()
    
    output_data = []
    for resnum, resname in zip(resnums, resnames):
        wt_residue = aa_names.get(resname, 'X')
        if wt_residue != 'X':
            for mut_residue in aa_codes:
                wt_psic, wt_group_psic, wt_group = calculate_psic(precomputed_folder, pdb_file, resnum, wt_residue, amino_acid_groups1)
                wt_size_psic, _, _ = calculate_psic(precomputed_folder, pdb_file, resnum, wt_residue, amino_acid_groups2)
                mut_psic, mut_group_psic, mut_group = calculate_psic(precomputed_folder, pdb_file, resnum, mut_residue, amino_acid_groups1)
                mut_size_psic, _, _ = calculate_psic(precomputed_folder, pdb_file, resnum, mut_residue, amino_acid_groups2)
                wt_mut_psic = wt_psic - mut_psic
                wt_mut_group_psic = wt_group_psic - mut_group_psic
                wt_mut_size_psic = wt_size_psic - mut_size_psic
                output_data.append([
                    pdb_file_name,
                    resnum,
                    wt_residue,
                    mut_residue,
                    f"{wt_psic:.3f}",
                    f"{wt_group_psic:.3f}",
                    f"{wt_size_psic:.3f}",
                    f"{mut_psic:.3f}",
                    f"{mut_group_psic:.3f}",
                    f"{mut_size_psic:.3f}",
                    f"{wt_mut_psic:.3f}",
                    f"{wt_mut_group_psic:.3f}",
                    f"{wt_mut_size_psic:.3f}"
                ])

    with open(output_file, "w") as f:
        f.write("PDB_File\tResidue_Number\tWT_Residue\tMutation\tWT_Seq_PSIC\tWT_PhyChem_Group_Seq_PSIC\tWT_Size_Group_Seq_PSIC\tMut_Seq_PSIC\tMut_PhyChem_Group_Seq_PSIC\tMut_Size_Group_Seq_PSIC\tWT-Mut_Seq_PSIC\tWT-Mut_PhyChem_Group_Seq_PSIC\tWT-Mut_Size_Group_Seq_PSIC\n")
        for row in output_data:
            f.write("\t".join(map(str, row)) + "\n")

def main():
    if len(sys.argv) != 3:
        print("Usage: python script.py -f <pdb_file>")
        return

    if sys.argv[1] != '-f':
        print("Invalid argument. Use -f to specify the PDB file.")
        return

    pdb_file = sys.argv[2]
    pdb_name = os.path.splitext(os.path.basename(pdb_file))[0]

    if not os.path.isfile(pdb_file):
        print(f"PDB file {pdb_file} does not exist.")
        return

    precomputed_folder = PSIC_folder  # Use the PSIC folder path
    precomputed_file = os.path.join(precomputed_folder, f"{pdb_name}.seqmsa.json")
    output_file = f"saturation_mutagenesis_psic.tsv"
    blast_file = os.path.join(BLAST_folder, f"{pdb_name}_blast.xml")
    msa_folder_path = MSA_folder  # Use the MSA folder path
    msa_out_path = os.path.join(msa_folder_path, f"{pdb_name}.seqmsa")

    if not os.path.exists(precomputed_folder):
        os.makedirs(precomputed_folder)

    try:
        print("Extracting sequence from PDB file...")
        sequence = extract_sequence_from_pdb(pdb_file)

        # Check if the JSON (PSIC) file exists and validate it
        if os.path.isfile(precomputed_file):
            print(f"Validating existing positional frequencies file with PDB sequence...")
            valid = validate_seqmsa_and_json_with_pdb(msa_out_path, precomputed_file, sequence)
            if valid:
                print("JSON validation successful, skipping MSA and BLAST steps.")
            else:
                raise ValueError("Validation failed: The positional frequencies file does not correspond to the present input.")
        else:
            # If JSON file does not exist, check for MSA file and proceed
            print("PSIC file not found, checking MSA file...")
            if os.path.isfile(msa_out_path):
                print(f"Validating existing MSA file with PDB sequence...")
                valid = validate_seqmsa_and_json_with_pdb(msa_out_path, precomputed_file, sequence)
                if not valid:
                    raise ValueError("Validation failed: The MSA file saved with the same name does not correspond to the present input.")
            else:
                # If MSA file is also not found, perform BLAST search to generate MSA
                print("MSA file not found. Performing BLAST search to generate MSA file...")
                perform_blast_search(sequence, blast_file=blast_file)
                parse_blast_msa(blast_file, sequence, msa_out_path)
                print("Precomputing frequencies...")
                precomputed_data = precompute_frequencies(msa_out_path)
                save_precomputed_frequencies(precomputed_folder, pdb_name, precomputed_data)

            # Check or generate the BLAST file if necessary
            if os.path.isfile(blast_file):
                print(f"Validating existing BLAST XML file...")
                valid = validate_blast_xml_with_pdb(blast_file, len(sequence))
                if not valid:
                    raise ValueError(f"Validation failed: The BLAST XML file {blast_file} saved with the same name does not correspond to the present input.")
            else:
                print("BLAST XML file not found. Performing BLAST search...")
                perform_blast_search(sequence, blast_file=blast_file)

        # Generate saturation mutagenesis features
        generate_saturation_mutagenesis_features(pdb_file, sequence, precomputed_folder, output_file)
        print(f"Features for saturation mutagenesis saved in {output_file}")
    except FileNotFoundError as fnf_error:
        print(f"File not found: {fnf_error}")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    main()

