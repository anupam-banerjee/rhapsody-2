import os
import shutil
import subprocess
import time
import sys
import pandas as pd
import re
from collections import defaultdict

def run_command(command):
    try:
        subprocess.run(command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Command failed: {e}")

def copy_pdb_file(pdb_file, dest_dir):
    if os.path.abspath(pdb_file) != os.path.join(dest_dir, os.path.basename(pdb_file)):
        shutil.copy(pdb_file, dest_dir)
        print(f"File {pdb_file} copied to {dest_dir}")

def check_file_exists(file_path, timeout=300):
    start_time = time.time()
    while not os.path.exists(file_path):
        if time.time() - start_time > timeout:
            return False
        time.sleep(5)
    return True

def background_file_check_and_abort(process, file_path, timeout=300):
    if not check_file_exists(file_path, timeout):
        process.terminate()
        print(f"{os.path.basename(file_path)} not created within {timeout} seconds. Aborting.")
        return False
    return True

def concatenate_files(output_file, *files):
    dataframes = []
    for file in files:
        df = pd.read_csv(file, sep="\t")
        if 'dyn' in file or 'shannon' in file or 'rsa_blosum' in file:
            df = df.iloc[:, 4:]
        dataframes.append(df)
    concatenated_df = pd.concat(dataframes, axis=1)
    concatenated_df.to_csv(output_file, sep="\t", index=False)
    print(f"Consolidated file saved as {output_file}")

def delete_files(*files):
    for file in files:
        if os.path.exists(file):
            os.remove(file)
            print(f"Deleted {file}")
        else:
            print(f"{file} not found")

def renumber_pdb_residues(input_pdb, original_pdb_backup):
    forward_map = {}
    residue_counter = 1
    current_residue = None
    residue_seen = set()

    shutil.copy(input_pdb, original_pdb_backup)

    with open(input_pdb, 'r') as infile:
        lines = infile.readlines()

    with open(input_pdb, 'w') as outfile:
        for line in lines:
            if line.startswith(('ATOM', 'HETATM')):
                res_num = line[22:26].strip()

                if res_num != current_residue:
                    if res_num not in residue_seen:
                        forward_map[residue_counter] = res_num
                        residue_seen.add(res_num)
                        current_residue = res_num
                        residue_counter += 1

                new_res_num_str = f"{residue_counter - 1:4}"
                new_line = line[:22] + new_res_num_str + line[26:]
                outfile.write(new_line)
            else:
                outfile.write(line)

    return forward_map

def reverse_map_predictions(results_file, forward_map, output_file):
    df = pd.read_csv(results_file, sep="\t", header=None)

    def map_residue_number(res_num):
        res_num_int = int(res_num)
        original_res = forward_map.get(res_num_int, '')
        return original_res if original_res else 'NA'

    df[1] = df[1].apply(map_residue_number)

    df.to_csv(output_file, sep="\t", header=False, index=False)

def main(pdb_file):
    current_dir = os.getcwd()
    original_pdb_backup = pdb_file.replace('.pdb', '_orinum.pdb')
    forward_map = renumber_pdb_residues(pdb_file, original_pdb_backup)

    copy_pdb_file(pdb_file, current_dir)

    pdb_filename = os.path.basename(pdb_file)
    run_command(f"python dyn_saturation.py -f \"{pdb_filename}\" -o saturation_mutagenesis_dyn.tsv")

    psic_process = subprocess.Popen(f"python psic_saturation.py -f \"{pdb_filename}\"", shell=True)
    if not background_file_check_and_abort(psic_process, "saturation_mutagenesis_psic.tsv"):
        return
    psic_process.wait()

    shannon_process = subprocess.Popen(f"python shannon_saturation.py -f \"{pdb_filename}\"", shell=True)
    if not background_file_check_and_abort(shannon_process, "saturation_mutagenesis_shannon.tsv"):
        return
    shannon_process.wait()

    run_command(f"python rsa_blosum_saturation.py -f \"{pdb_filename}\" -o saturation_mutagenesis_rsa_blosum.tsv")

    concatenate_files("saturation_mutagenesis_red.tsv", "saturation_mutagenesis_psic.tsv", 
                      "saturation_mutagenesis_shannon.tsv", "saturation_mutagenesis_dyn_selected21.tsv", 
                      "saturation_mutagenesis_rsa_blosum.tsv")

    run_command("python predict_red.py")
    shutil.move("predictions.txt", "predictions_nonmapped.txt")

    reverse_map_predictions("predictions_nonmapped.txt", forward_map, "predictions.txt")

    run_command("python heatmap_withavg.py")
    run_command(f"python replace_bfactor.py \"{original_pdb_backup}\"")
    
    # Restore original PDB file
    shutil.move(original_pdb_backup, pdb_filename)

    # Clean up intermediate files
    delete_files("saturation_mutagenesis_dyn.tsv",
                 "saturation_mutagenesis_psic.tsv",
                 "saturation_mutagenesis_shannon.tsv",
                 "saturation_mutagenesis_rsa_blosum.tsv",
                 "saturation_mutagenesis_red.tsv",
                 "predictions_nonmapped.txt")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <pdb_file>")
        sys.exit(1)
    main(sys.argv[1])

