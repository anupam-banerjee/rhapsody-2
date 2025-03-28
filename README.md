Rhapsody-2

Overview

This repository contains Python scripts designed to predict the pathogenicity probabilities for saturation mutagenesis, covering all possible single amino acid variations at every position in target proteins. Additionally, the Rhapsody-2 SAV-DB is provided in the database folder.

Features

DIAMOND BLAST Integration: The scripts utilize DIAMOND to perform BLAST searches on the nr database with an e-value of 0.001, specifically targeting human protein structures resolved by AlphaFold2. All hits meeting this e-value threshold are used to calculate PSIC and Shannon entropy features. To facilitate integration, the pipeline approximates this process by conducting a BLAST search over the internet using the same e-value and a maximum hit size of 50,000. This step can be time-consuming depending on the protein size and server load. A timeout mechanism ensures that if BLAST results are not generated within 5 minutes, the script defaults to dynamics-only predictions. Precomputed PSIC and Shannon entropy files for human proteins will be made available in future releases to expedite this process.

### Optional: Manual BLAST Search

Alternatively, users can perform the BLAST search manually at [NCBI BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE=Proteins) using the **nr** database. Set the **E-value threshold** to `0.001` and adjust the **Max target sequences** to `5000` under the *Algorithm parameters* section. 

After running the search:

1. Download the results in **XML format**.
2. Rename the file to `{pdb_name}_blast.xml`.
3. Place it in the `lib/blast` directory.

This allows users to skip running the BLAST search programmatically over the internet. However, note that the web-based BLAST has a hit limit of 5000, which may not be sufficient when compared to a DIAMOND search with the same E-value threshold. In such cases, the limited output may not fully approximate the results obtained via DIAMOND.


File Management and Cleanup: The rhapsody2_red.py script ensures that all necessary files are present, copying them to the current directory if needed, and removing intermediate files after processing is complete.
Feature Calculation: The script computes various dynamics and evolutionary features, including PSIC and Shannon entropy, which are essential for pathogenicity predictions.
Consolidation and Prediction: After feature computation, the results are consolidated into a single file, which is then used for final pathogenicity predictions.
Error Handling: The script includes safeguards to abort processes if necessary files are not generated within a designated timeout, ensuring that the pipeline can still produce results even if some steps fail.
Requirements

Python Packages
The pipeline relies on several Python packages, including ProDy, NumPy, Pandas, Biopython, Matplotlib, and Seaborn. These packages are essential for tasks such as protein data handling, running BLAST searches, data processing, and visualization.

Additional Requirements
The script also requires the Stride software for secondary structure assignments. Stride can be downloaded from the official website: Stride. Ensure that the Stride executable is properly configured and accessible on your system.

Setting Up Paths
Before running the pipeline, paths need to be configured correctly. Specifically, you must provide the path to the lib folder and specify the path to the Stride executable in the paths.py file located in the ./code/utilities/ directory. Proper configuration of these paths is crucial for the pipeline's functionality.

Precomputed Files
To save time and computational resources, precomputed PSIC and Shannon entropy files, along with sample input files, are provided in the lib folder. These files are currently used for computations, allowing you to bypass the time-consuming BLAST search steps.

Usage

To run the script, provide the path to your PDB file as an argument when executing the script.

php
Copy code
python rhapsody2_red.py <pdb_file>
Example:

Copy code
python rhapsody2_red.py AF-A2A3L6-F1-model_v4.pdb
This will start the pipeline, which copies the PDB file to the current working directory if it is not already there, and then proceeds with the following steps: dynamics saturation, PSIC computation, Shannon entropy calculation, RSA and BLOSUM analysis, feature consolidation, and pathogenicity prediction.

Script Workflow

Copy PDB File: The script begins by copying the specified PDB file to the current working directory if it is not already present.
Run Dynamics Saturation: The script executes the dyn_saturation.py script to compute saturation mutagenesis dynamics.
PSIC Computation: The script initiates the psic_saturation.py script in the background, monitoring the creation of the output file. If the file is not generated within 300 seconds, the script aborts the process and defaults to dynamics-only predictions.
Shannon Entropy Calculation: Following a similar procedure, the script executes the shannon_saturation.py script and monitors the output file for Shannon entropy calculations.
RSA and BLOSUM Analysis: Upon successful completion of the previous steps, the script runs the rsa_blosum_saturation.py script for RSA and BLOSUM analysis.
Feature Consolidation: Once all features have been computed, the script consolidates the results into a single file.
Pathogenicity Prediction: The consolidated file is then used to predict the pathogenicity of single amino acid variants, with the results being integrated back into the PDB file.
Cleanup: Finally, the script deletes any intermediate files to conserve space and maintain a clean working directory.
This workflow ensures comprehensive analysis while maintaining efficient resource use and error management.

Results

The results of the saturation mutagenesis predictions will be stored in a file named predictions.txt. Additionally, the predictions will be visualized as heatmaps, generated in chunks of 100 residues. These heatmaps will help in understanding the pathogenicity distribution across the protein structure. The average pathogenicity values for each residue will also be written to a modified PDB file named input_file_bfac.pdb, where the B-factor column will be replaced with the computed pathogenicity values.

This comprehensive workflow ensures thorough analysis and effective visualization of the pathogenicity predictions, providing valuable insights into the impact of single amino acid variations on protein function.
