import argparse
import pandas as pd
import prody as pr
import os
import subprocess
import logging
import re
import numpy as np
import math
import shutil
from dynamics_utils import cal_MSF, cal_effecsens, cal_meancorr, cal_hingedist

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Define the amino acid codes
aa_codes = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
aa_names = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', 'GLY': 'G', 'HIS': 'H',
    'ILE': 'I', 'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q',
    'ARG': 'R', 'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
}

def preprocess_pdb(pdb_path):
    """
    Extracts ATOM lines from a PDB file, saving them back to the original file after using a temporary file.
    """
    temp_pdb_path = os.path.join(os.getcwd(), "temp.pdb")
    
    try:
        # Copy the original PDB file to temp.pdb
        shutil.copy(pdb_path, temp_pdb_path)
        
        # Use grep to extract ATOM lines from temp.pdb and save to the original PDB file
        with open(pdb_path, 'w') as outfile:
            subprocess.run(['grep', '^ATOM', temp_pdb_path], stdout=outfile)
        
        # Remove the temporary file
        os.remove(temp_pdb_path)
        
        return pdb_path
    except Exception as e:
        logging.error(f"Error processing PDB file {pdb_path}: {e}")
        return None

def generate_saturation_mutagenesis_input(pdb_file_path, output_file_path):
    """
    Generates input for saturation mutagenesis.
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

def analyze_residue(row, output_file, pdb_folder=None):
    """
    Analyzes a single residue.
    """
    pdb_path = row['pdb_file_path']
    preprocessed_pdb_path = preprocess_pdb(pdb_path)
    if not preprocessed_pdb_path or not os.path.exists(preprocessed_pdb_path):
        logging.error(f"Error: PDB file not found or preprocessing failed: {pdb_path}")
        result = [row['pdb_file_path'], row['resnum'], row['wt_residue'], row['mut']] + ["NA"] * 80
        with open(output_file, 'a') as f:
            f.write("\t".join(map(str, result)) + "\n")
        return
    
    results = [row['pdb_file_path'], row['resnum'], row['wt_residue'], row['mut']]
    chain = 'A'  # Assuming chain A; modify as needed
    mut_res = row['resnum']
    
    try:
        results.extend(calculate_effectiveness_sensitivity(preprocessed_pdb_path, chain, mut_res))
        results.extend(calculate_mean_correlation(preprocessed_pdb_path, chain, mut_res))
        results.extend(calculate_MSF(preprocessed_pdb_path, chain, mut_res))
        results.extend(calculate_hinge_distance(preprocessed_pdb_path, chain, mut_res))
    except Exception as e:
        logging.error(f"Error processing residue {mut_res} in file {pdb_path}: {e}")
        results.extend(["NA"] * (84 - len(results)))  # Ensure the length of results matches the expected columns
    finally:
        if os.path.exists(preprocessed_pdb_path):
            #os.remove(preprocessed_pdb_path)
            pass
    
    logging.info(f"Length of results for residue {mut_res}: {len(results)}")
    
    # Write the result to the output file
    with open(output_file, 'a') as f:
        f.write("\t".join(map(str, results)) + "\n")

def calculate_effectiveness_sensitivity(preprocessed_pdb_path, chain, mut_res):
    """
    Calculates effectiveness and sensitivity.
    """
    results = []
    modes = ['soft', 'all', 5]
    enms = ['ANM', 'GNM']
    for enm in enms:
        for mode in modes:
            try:
                effec_sens = cal_effecsens(preprocessed_pdb_path, chain, mut_res, n_modes=mode, enm=enm)
                results.extend([f'{val:.3f}' for val in effec_sens])
            except Exception as e:
                logging.error(f"Error calculating effectiveness and sensitivity for {enm} mode {mode}: {e}")
                results.extend(["NA"] * 6)  # Extend with NA if there's an error
    return results

def calculate_mean_correlation(preprocessed_pdb_path, chain, mut_res):
    """
    Calculates mean correlation.
    """
    results = []
    modes = ['all', 'soft', 5]
    enms = ['GNM', 'ANM']
    for enm in enms:
        for mode in modes:
            try:
                mean_corr = cal_meancorr(preprocessed_pdb_path, chain, mut_res, n_modes=mode, enm=enm)
                results.extend([f'{val:.3f}' for val in mean_corr])
            except Exception as e:
                logging.error(f"Error calculating mean correlation for {enm} mode {mode}: {e}")
                results.extend(["NA"] * 3)  # Extend with NA if there's an error
    return results

def calculate_MSF(preprocessed_pdb_path, chain, mut_res):
    """
    Calculates MSF.
    """
    results = []
    msf_modes = ['hot', 'soft', 'all', 5]
    enms = ['ANM', 'GNM']
    for enm in enms:
        for mode in msf_modes:
            try:
                msf = cal_MSF(preprocessed_pdb_path, chain, mut_res, n_modes=mode, enm=enm)
                results.extend([f'{val:.3f}' for val in msf])
            except Exception as e:
                logging.error(f"Error calculating MSF for {enm} mode {mode}: {e}")
                results.extend(["NA"] * 3)  # Extend with NA if there's an error
    return results

def calculate_hinge_distance(preprocessed_pdb_path, chain, mut_res):
    """
    Calculates hinge distance binary.
    """
    results = []
    hinge_modes = [5, 'soft']
    for mode in hinge_modes:
        try:
            hinge_dist = cal_hingedist(preprocessed_pdb_path, chain, mut_res, n_modes=mode, enm='GNM')
            results.append(f'{hinge_dist:.3f}')
        except Exception as e:
            logging.error(f"Error calculating hinge distance for mode {mode}: {e}")
            results.append("NA")  # Extend with NA if there's an error
    return results

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
    columns = [
        'pdb_file_path', 'resnum', 'wt_residue', 'mut',
        'res_effectiveness (ANM soft)', 'res_zscore_effectiveness (ANM soft)', 'res_norm_zscore_effectiveness (ANM soft)', 
        'res_sensitivity (ANM soft)', 'res_zscore_sensitivity (ANM soft)', 'res_norm_zscore_sensitivity (ANM soft)', 
        'res_effectiveness (ANM all)', 'res_zscore_effectiveness (ANM all)', 'res_norm_zscore_effectiveness (ANM all)', 
        'res_sensitivity (ANM all)', 'res_zscore_sensitivity (ANM all)', 'res_norm_zscore_sensitivity (ANM all)', 
        'res_effectiveness (ANM top5)', 'res_zscore_effectiveness (ANM top5)', 'res_norm_zscore_effectiveness (ANM top5)', 
        'res_sensitivity (ANM top5)', 'res_zscore_sensitivity (ANM top5)', 'res_norm_zscore_sensitivity (ANM top5)', 
        'res_effectiveness (GNM soft)', 'res_zscore_effectiveness (GNM soft)', 'res_norm_zscore_effectiveness (GNM soft)', 
        'res_sensitivity (GNM soft)', 'res_zscore_sensitivity (GNM soft)', 'res_norm_zscore_sensitivity (GNM soft)', 
        'res_effectiveness (GNM all)', 'res_zscore_effectiveness (GNM all)', 'res_norm_zscore_effectiveness (GNM all)', 
        'res_sensitivity (GNM all)', 'res_zscore_sensitivity (GNM all)', 'res_norm_zscore_sensitivity (GNM all)', 
        'res_effectiveness (GNM top5)', 'res_zscore_effectiveness (GNM top5)', 'res_norm_zscore_effectiveness (GNM top5)', 
        'res_sensitivity (GNM top5)', 'res_zscore_sensitivity (GNM top5)', 'res_norm_zscore_sensitivity (GNM top5)', 
        'res_avg_meancorr (GNM all)', 'res_avg_meancorr_zscore (GNM all)', 'res_avg_meancorr_norm_zscore (GNM all)', 
        'res_avg_meancorr (GNM soft)', 'res_avg_meancorr_zscore (GNM soft)', 'res_avg_meancorr_norm_zscore (GNM soft)', 
        'res_avg_meancorr (GNM top5)', 'res_avg_meancorr_zscore (GNM top5)', 'res_avg_meancorr_norm_zscore (GNM top5)', 
        'res_avg_meancorr (ANM all)', 'res_avg_meancorr_zscore (ANM all)', 'res_avg_meancorr_norm_zscore (ANM all)', 
        'res_avg_meancorr (ANM soft)', 'res_avg_meancorr_zscore (ANM soft)', 'res_avg_meancorr_norm_zscore (ANM soft)', 
        'res_avg_meancorr (ANM top5)', 'res_avg_meancorr_zscore (ANM top5)', 'res_avg_meancorr_norm_zscore (ANM top5)', 
        'res_MSF (ANM hot)', 'res_MSF_zscore (ANM hot)', 'res_MSF_norm_zscore (ANM hot)', 
        'res_MSF (ANM soft)', 'res_MSF_zscore (ANM soft)', 'res_MSF_norm_zscore (ANM soft)', 
        'res_MSF (ANM all)', 'res_MSF_zscore (ANM all)', 'res_MSF_norm_zscore (ANM all)', 
        'res_MSF (ANM top5)', 'res_MSF_zscore (ANM top5)', 'res_MSF_norm_zscore (ANM top5)', 
        'res_MSF (GNM hot)', 'res_MSF_zscore (GNM hot)', 'res_MSF_norm_zscore (GNM hot)', 
        'res_MSF (GNM soft)', 'res_MSF_zscore (GNM soft)', 'res_MSF_norm_zscore (GNM soft)', 
        'res_MSF (GNM all)', 'res_MSF_zscore (GNM all)', 'res_MSF_norm_zscore (GNM all)', 
        'res_MSF (GNM top5)', 'res_MSF_zscore (GNM top5)', 'res_MSF_norm_zscore (GNM top5)', 
        'res_hinge_binary (GNM top5)', 'res_hinge_binary (GNM soft)'
    ]
    with open(args.output, 'w') as f:
        f.write("\t".join(columns) + "\n")
    
    for _, row in df.iterrows():
        logging.info(f"Processing {row['pdb_file_path']} residue {row['resnum']}...")
        analyze_residue(row, args.output, pdb_folder=args.pdbfolder)

if __name__ == '__main__':
    main()

