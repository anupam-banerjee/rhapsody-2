import os
import subprocess
import logging
import pandas as pd
import argparse
from utilities.paths import Stride_folder, stride

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Define the amino acid codes and names
aa_codes = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
aa_names = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', 'GLY': 'G', 'HIS': 'H',
    'ILE': 'I', 'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q',
    'ARG': 'R', 'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
}

# Sander values for RSA calculation
sander_values = {
    'ALA': 106, 'CYS': 135, 'ASP': 163, 'GLU': 194, 'PHE': 197, 'GLY': 84, 'HIS': 184,
    'ILE': 169, 'LYS': 205, 'LEU': 164, 'MET': 188, 'ASN': 157, 'PRO': 136, 'GLN': 198,
    'ARG': 248, 'SER': 130, 'THR': 142, 'VAL': 142, 'TRP': 227, 'TYR': 222
}

def generate_saturation_mutagenesis_input(pdb_file_path, output_file_path):
    """
    Generates input for saturation mutagenesis by creating all possible single amino acid mutations.
    """
    try:
        structure = pr.parsePDB(pdb_file_path, subset='ca')
        if structure is None:
            logging.error(f"Error: Structure is None for file {pdb_file_path}")
            return
        resnums = structure.getResnums()
        resnames = structure.getResnames()
    except Exception as e:
        logging.error(f"Error parsing PDB file {pdb_file_path}: {e}")
        return
    
    with open(output_file_path, 'w') as f:
        for resnum, resname in zip(resnums, resnames):
            wt_residue = aa_names.get(resname, 'X')
            if wt_residue != 'X':
                for mut in aa_codes:
                    f.write(f"{pdb_file_path}\t{resnum}\t{wt_residue}\t{mut}\n")

def generate_stride_file(pdb_file, stride_folder, stride_path):
    """
    Generate STRIDE file for the given PDB file if it doesn't exist.
    """
    stride_output_file = os.path.join(stride_folder, f"{os.path.basename(pdb_file).replace('.pdb', '')}.strd")
    
    if not os.path.isfile(stride_output_file):
        logging.info(f"STRIDE file not found for {pdb_file}. Generating...")
        try:
            with open(stride_output_file, 'w') as outfile:
                subprocess.run([stride_path, pdb_file], stdout=outfile)
        except Exception as e:
            logging.error(f"Error generating STRIDE file for {pdb_file}: {e}")
            return None
    else:
        logging.info(f"STRIDE file found for {pdb_file}.")
    
    return stride_output_file

def calculate_rsa(residue, stride_file):
    """
    Calculate the RSA for a given residue using the STRIDE file.
    """
    with open(stride_file, 'r') as fp:
        for line in fp:
            if line.startswith('ASG'):
                line_data = line.split()
                res_name = line_data[1]
                res_num = line_data[3]
                solvent_access = float(line_data[9])

                if res_num == str(residue['resnum']):
                    rsa = solvent_access / sander_values.get(res_name, 1)
                    return round(rsa, 2)
    
    return 'NaN'

def calculate_blosum62(wt_residue, mut_residue):
    """
    Calculate the BLOSUM62 score using an external script.
    """
    try:
        blosum_score = subprocess.check_output(['bash', 'blosum.sh', wt_residue, mut_residue])
        return blosum_score.strip().decode('utf-8')
    except subprocess.CalledProcessError as e:
        logging.error(f"Error calculating BLOSUM62 score: {e}")
        return 'NaN'

def analyze_residue(row, output_file, stride_folder, stride_path):
    """
    Analyze a single residue, including RSA and BLOSUM62 calculations.
    """
    pdb_path = row['pdb_file_path']
    
    # Generate STRIDE file if not present
    stride_file = generate_stride_file(pdb_path, stride_folder, stride_path)
    if not stride_file:
        logging.error(f"Failed to generate STRIDE file for {pdb_path}")
        result = [row['pdb_file_path'], row['resnum'], row['wt_residue'], row['mut'], "NaN", "NaN"]
        with open(output_file, 'a') as f:
            f.write("\t".join(map(str, result)) + "\n")
        return
    
    # Calculate RSA
    rsa = calculate_rsa(row, stride_file)
    
    # Calculate BLOSUM62 score
    blosum_score = calculate_blosum62(row['wt_residue'], row['mut'])
    
    # Prepare the result row
    results = [row['pdb_file_path'], row['resnum'], row['wt_residue'], row['mut'], rsa, blosum_score]
    
    # Write the result to the output file
    with open(output_file, 'a') as f:
        f.write("\t".join(map(str, results)) + "\n")

def main():
    """
    Main function to process input and perform analysis.
    """
    parser = argparse.ArgumentParser(description='Process PDB file or input file for residue analysis.')
    parser.add_argument('-o', '--output', required=True, help='Output file path')
    parser.add_argument('-i', '--input', help='Input file path (optional)')
    parser.add_argument('-f', '--pdbfile', help='PDB file path (optional)')
    parser.add_argument('-d', '--pdbfolder', help='Folder path containing PDB files (optional)')

    args = parser.parse_args()
    
    if not (args.input or args.pdbfile):
        parser.error('Either --input or --pdbfile must be specified.')

    if args.pdbfile:
        saturation_mutagenesis_input_path = 'saturation_mutagenesis_input.csv'
        generate_saturation_mutagenesis_input(args.pdbfile, saturation_mutagenesis_input_path)
        input_file_path = saturation_mutagenesis_input_path
    else:
        input_file_path = args.input

    df = pd.read_csv(input_file_path, names=['pdb_file_path', 'resnum', 'wt_residue', 'mut'], sep='\t')

    # Write the header to the output file
    columns = ['pdb_file_path', 'resnum', 'wt_residue', 'mut', 'RSA', 'BLOSUM62']
    with open(args.output, 'w') as f:
        f.write("\t".join(columns) + "\n")
    
    for _, row in df.iterrows():
        logging.info(f"Processing {row['pdb_file_path']} residue {row['resnum']}...")
        analyze_residue(row, args.output, stride_folder=Stride_folder, stride_path=stride)

if __name__ == '__main__':
    main()

