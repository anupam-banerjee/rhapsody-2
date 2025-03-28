# Rhapsody-2

## Overview

This repository contains Python scripts designed to predict pathogenicity probabilities for saturation mutagenesis, covering all possible single amino acid variations at every position in target proteins. Additionally, the Rhapsody-2 SAV-DB is provided in the `database` folder.

## Features

### DIAMOND BLAST Integration

The scripts utilize DIAMOND to perform BLAST searches on the **nr** database with an E-value threshold of `0.001`, specifically targeting human protein structures resolved by AlphaFold2. All hits meeting this threshold are used to calculate PSIC and Shannon entropy features.

To approximate this process, the pipeline can also perform a BLAST search over the internet using the same E-value and a maximum hit size of 50,000. This step may be time-consuming depending on protein size and server load. A timeout mechanism ensures that if BLAST results are not generated within 5 minutes, the script defaults to dynamics-only predictions.

Precomputed PSIC and Shannon entropy files for human proteins will be made available in future releases to expedite this process.

### Optional: Manual BLAST Search

Alternatively, users can perform the BLAST search manually via [NCBI BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE=Proteins) using the **nr** database. Set the following parameters in the *Algorithm parameters* section:

- **E-value threshold**: `0.001`
- **Max target sequences**: `5000`

After running the search:

1. Download the results in **XML format**.
2. Rename the file to `{pdb_name}_blast.xml`.
3. Place it in the `lib/blast` directory.

This bypasses the need to run BLAST over the internet within the pipeline. However, note that the web-based BLAST has a maximum hit limit of 5000, which may be insufficient in cases where DIAMOND yields a (much) larger set of hits based on the same E-value. Therefore, this manual method might not fully approximate the results of a DIAMOND search.

### File Management and Cleanup

The `rhapsody2_red.py` script ensures all necessary files are available, copying them to the current directory if needed and removing intermediate files after processing.

### Feature Calculation

The script computes various dynamics and evolutionary features, including PSIC and Shannon entropy, essential for pathogenicity predictions.

### Consolidation and Prediction

All computed features are consolidated into a single file, which is then used for final pathogenicity predictions.

### Error Handling

A timeout safeguard aborts steps if required files are not generated within a set time, ensuring the pipeline can fall back to dynamics-only predictions when necessary.

## Requirements

### Python Packages

The pipeline requires the following Python libraries:

- `ProDy`
- `NumPy`
- `Pandas`
- `Biopython`
- `Matplotlib`
- `Seaborn`

### Additional Requirements

- **Stride**  
  The pipeline depends on [Stride](http://webclu.bio.wzw.tum.de/stride/) for secondary structure assignments. Ensure the Stride executable is downloaded and its path is properly set.

## Setting Up Paths

Before running the pipeline, update the paths in the `paths.py` file located in `./code/utilities/`:

- Path to the `lib` folder
- Path to the `Stride` executable

## Precomputed Files

To reduce compute time, precomputed PSIC and Shannon entropy files are provided in the `lib` folder. These allow the pipeline to skip the BLAST search step for known structures.

## Usage

To run the pipeline, provide the path to your PDB file:

```bash
python rhapsody2_red.py <pdb_file>
```

**Example:**

```bash
python rhapsody2_red.py AF-A2A3L6-F1-model_v4.pdb
```

## Script Workflow

1. **Copy PDB File**: Ensures input file is in the working directory.
2. **Run Dynamics Saturation**: Executes `dyn_saturation.py`.
3. **PSIC Computation**: Runs `psic_saturation.py`. If no output is detected within 5 minutes, it skips to dynamics-only.
4. **Shannon Entropy Calculation**: Executes `shannon_saturation.py` with similar timeout handling.
5. **RSA and BLOSUM Analysis**: Runs `rsa_blosum_saturation.py`.
6. **Feature Consolidation**: Combines all computed features into a unified file.
7. **Pathogenicity Prediction**: Predicts mutation impacts and writes results back into the PDB.
8. **Cleanup**: Deletes intermediate files to keep things tidy.

## Results

- Results are saved in a file named `predictions.txt`.
- Heatmaps are generated (in 100-residue chunks) to visualize mutation impact.
- A modified PDB file (`input_file_bfac.pdb`) is created with the B-factor column overwritten by predicted pathogenicity values, providing a structure-based visualization.

## Precomputed features and pathogenicity predictions in future releases

Precomputed PSIC and Shannon entropy files for a majority of human proteins will be made available for download in future releases to streamline the pipeline and reduce runtime. Additionally, precomputed pathogenicity predictions generated using Rhapsody-2 (red) will also be provided, offering users immediate access to mutation impact insights across a wide range of human proteins.
