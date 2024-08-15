import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.gridspec as gridspec


plt.rcParams.update({'font.size': 16})

# Use Arial font
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'

# Load predictions
predictions = pd.read_csv('predictions.txt', delimiter='\t', header=None, names=['file_name', 'residue_index', 'WT_residue', 'mutant_residue', 'probability'])

# Define a custom colormap with blue to red gradient
custom_cmap = LinearSegmentedColormap.from_list('custom_bwr', ['#89CFF0', '#FF6961'])

# Function to create heat maps
def create_heatmaps(predictions, chunk_size=100):
    # Create a pivot table with residue_index and mutant_residue
    pivot_table = predictions.pivot(index='mutant_residue', columns='residue_index', values='probability')
    
    # Set wild type residues substituted with themselves to NaN
    for idx, row in predictions.iterrows():
        if row['WT_residue'] == row['mutant_residue']:
            pivot_table.at[row['mutant_residue'], row['residue_index']] = np.nan
    
    # Get unique residue indices and sort them
    residue_indices = sorted(predictions['residue_index'].unique())
    
    # Determine the number of chunks
    num_chunks = (len(residue_indices) + chunk_size - 1) // chunk_size
    
    # Amino acids order
    amino_acids = list("ACDEFGHIKLMNPQRSTVWY")
    
    # Iterate over each chunk to create heat maps
    for i in range(num_chunks):
        start_idx = i * chunk_size
        end_idx = min((i + 1) * chunk_size, len(residue_indices))
        chunk_residues = residue_indices[start_idx:end_idx]
        
        # Extract the relevant part of the pivot table
        chunk_data = pivot_table[chunk_residues]
        
        # Reindex to ensure all amino acids are present
        chunk_data = chunk_data.reindex(amino_acids)
        
        # Calculate the average for each column excluding NaN values
        column_averages = chunk_data.mean(axis=0)
        
        # Create a figure with subplots
        fig = plt.figure(figsize=(22, 14))
        gs = gridspec.GridSpec(2, 1, height_ratios=[10, 2])
        
        # Create heat map in the top subplot
        ax1 = plt.subplot(gs[0])
        sns.heatmap(chunk_data, cmap=custom_cmap, cbar_kws={'label': 'Probability'}, vmin=0, vmax=1, mask=chunk_data.isna(), linewidths=.6, ax=ax1)
        
        # Adjust font size for colorbar label and ticks
        colorbar = ax1.collections[0].colorbar
        colorbar.ax.tick_params(labelsize=20)  # Font size for colorbar ticks
        colorbar.set_label('Probability', size=24)  # Font size for colorbar label
        
        # Set y-axis labels for the heatmap
        ax1.set_yticks(np.arange(0.5, len(amino_acids), 1))
        ax1.set_yticklabels(amino_acids, rotation=0, fontsize=20)
        
        wt_residue_map = predictions.drop_duplicates(subset='residue_index').set_index('residue_index')['WT_residue'].to_dict()
        chunk_labels = [f"{wt_residue_map[index]}{index}" for index in chunk_residues]
        
        # Set titles and labels for the heatmap
        ax1.set_title(f'Heat Map for Residues {chunk_labels[0]} to {chunk_labels[-1]}', fontsize=28)
        ax1.set_ylabel('Amino acid type', fontsize=24)
        ax1.set_xlabel('Residue Index', fontsize=24)
        
        # Create a dictionary to map residue_index to WT_residue
        wt_residue_map = predictions.drop_duplicates(subset='residue_index').set_index('residue_index')['WT_residue'].to_dict()
        
        # Generate chunk_labels by appending WT_residue before the residue index
        chunk_labels = [f"{wt_residue_map[index]}{index}" for index in chunk_residues]
        
        # Plot the average values as a line plot in the bottom subplot
        ax2 = plt.subplot(gs[1])
        ax2.plot(np.arange(len(chunk_residues)) + 0.5, column_averages, marker='o', color='black', linestyle='-', linewidth=2, label='Average')
        
        # Set x-axis labels for the average plot
        ax2.set_xticks(np.arange(len(chunk_residues)) + 0.5)
        ax2.set_xticklabels([chunk_labels[i] if i % 2 == 0 else '' for i in range(len(chunk_labels))], rotation=90, fontsize=16)
        
               
        # Set y-axis limits and labels for the average plot
        ax2.set_ylim(0, 1)
        ax2.set_ylabel('Avg. Pathogenicity', fontsize=22)
        ax2.set_xlabel('Residue',fontsize=22)
        
        # Adjust the layout to make space for the colorbar and ensure alignment
        gs.tight_layout(fig)
        
        # Manually adjust the width of the bottom plot
        pos1 = ax1.get_position()  # Get the original position
        pos2 = ax2.get_position()  # Get the original position
        new_pos2 = [pos1.x0, pos2.y0, pos1.width, pos2.height]  # Set new position with the same width as ax1
        ax2.set_position(new_pos2)  # Apply the new position
        
        # Remove the extra space on the sides
        ax2.set_xlim(ax1.get_xlim())
        
        # Save the figure
        plt.savefig(f'heatmap_{chunk_residues[0]}_{chunk_residues[-1]}.png', bbox_inches='tight', dpi=300)
        plt.close()

# Create heat maps
create_heatmaps(predictions)

