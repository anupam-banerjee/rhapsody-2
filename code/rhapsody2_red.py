import os
import shutil
import subprocess
import time
import threading
import sys
import pandas as pd

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
        print(f"{os.path.basename(file_path)} not created within {timeout} seconds. Aborting and proceeding with dynamics only predictions.")
        return False
    return True

def concatenate_files(output_file, *files):
    dataframes = []
    for file in files:
        df = pd.read_csv(file, sep="\t")
        if 'dyn' in file or 'shannon' in file or 'rsa_blosum' in file:
            df = df.iloc[:, 4:]  # Select columns from the 5th to the last
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

def main(pdb_file):
    current_dir = os.getcwd()
    
    # Copy pdb file to the current directory if it's not there
    copy_pdb_file(pdb_file, current_dir)
    
    pdb_filename = os.path.basename(pdb_file)
    
    # Run the dyn_saturation.py command
    run_command(f"python dyn_saturation.py -f \"{pdb_filename}\" -o saturation_mutagenesis_dyn.tsv")
    
    dyn_old_file = "saturation_mutagenesis_dyn.tsv"
    dyn_file = "saturation_mutagenesis_dyn_selected21.tsv"
    psic_file = "saturation_mutagenesis_psic.tsv"
    shannon_file = "saturation_mutagenesis_shannon.tsv"
    rsa_blosum_file = "saturation_mutagenesis_rsa_blosum.tsv"
    red_file = "saturation_mutagenesis_red.tsv"
    max_min_file = "max_min_values.txt"
    input_csv = "saturation_mutagenesis_input.csv"
    
    # Start running psic_saturation.py
    psic_process = subprocess.Popen(f"python psic_saturation.py -f \"{pdb_filename}\"", shell=True)
    
    # Check for the creation of the PSIC file in the background
    if not background_file_check_and_abort(psic_process, psic_file, timeout=300):
        run_command("python predict_dyn.py")
        run_command("python heatmap_withavg.py")
        run_command(f"python replace_bfactor.py \"{pdb_filename}\"")
        # Delete the intermediate files
        delete_files(dyn_file, max_min_file, input_csv)
        return
    
    psic_process.wait()  # Ensure the process has finished

    # Start running shannon_saturation.py
    shannon_process = subprocess.Popen(f"python shannon_saturation.py -f \"{pdb_filename}\"", shell=True)
    
    # Check for the creation of the Shannon file in the background
    if not background_file_check_and_abort(shannon_process, shannon_file, timeout=300):
        run_command("python predict_dyn.py")
        run_command("python heatmap_withavg.py")
        run_command(f"python replace_bfactor.py \"{pdb_filename}\"")
        delete_files(dyn_file, max_min_file, input_csv, psic_file)
        return
    
    shannon_process.wait()  # Ensure the process has finished

    # Run rsa_blosum_saturation.py after successful creation of all previous files
    run_command(f"python rsa_blosum_saturation.py -f \"{pdb_filename}\" -o {rsa_blosum_file}")
    
    # Create the concatenated file
    concatenate_files(red_file, psic_file, shannon_file, dyn_file, rsa_blosum_file)
    
    print("All dynamics based and evolutionary features computed. Predicting pathogenicity based on reduced model.")
    run_command("python predict_red.py")
    run_command("python heatmap_withavg.py")
    run_command(f"python replace_bfactor.py \"{pdb_filename}\"")
    
    # Delete the intermediate files
    delete_files(dyn_file, psic_file, shannon_file, rsa_blosum_file, max_min_file, input_csv, red_file, dyn_old_file)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <pdb_file>")
        sys.exit(1)
    
    pdb_file = sys.argv[1]
    main(pdb_file)

