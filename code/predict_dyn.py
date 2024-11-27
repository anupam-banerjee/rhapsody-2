import numpy as np
import pandas as pd
from xgboost import XGBClassifier

# Function to load and process training data, handling Inf and -Inf values
def load_training_data(file_path):
    try:
        data = np.genfromtxt(file_path, delimiter='\t')
        # Check for Inf and -Inf values and replace them with the max and min of the respective column
        with np.errstate(invalid='ignore'):
            inf_mask = np.isinf(data)
            max_values = np.nanmax(np.where(np.isfinite(data), data, np.nan), axis=0)
            min_values = np.nanmin(np.where(np.isfinite(data), data, np.nan), axis=0)
            data[inf_mask & (data > 0)] = np.take(max_values, np.where(inf_mask & (data > 0))[1])
            data[inf_mask & (data < 0)] = np.take(min_values, np.where(inf_mask & (data < 0))[1])
        X = data[:, :-1]  # Features
        y = data[:, -1]   # Labels
        return X, y, max_values, min_values
    except Exception as e:
        print(f"Error loading data from {file_path}: {e}")
        return None, None, None, None

# Function to load and process test data, handling Inf and -Inf values using training data statistics
def load_test_data(file_path, max_values, min_values):
    try:
        data = pd.read_csv(file_path, delimiter='\t', header=None)
        metadata = data.iloc[1:, :4].values  # First four columns are metadata
        features = data.iloc[1:, 4:].values.astype(np.float32)  # Features start from the 5th column

        # Check for Inf and -Inf values and replace them with the max and min of the respective column from training data
        with np.errstate(invalid='ignore'):
            inf_mask = np.isinf(features)
            features[inf_mask & (features > 0)] = np.take(max_values, np.where(inf_mask & (features > 0))[1])
            features[inf_mask & (features < 0)] = np.take(min_values, np.where(inf_mask & (features < 0))[1])
        return metadata, features
    except Exception as e:
        print(f"Error loading data from {file_path}: {e}")
        return None, None

# Load training data
X_train, y_train, max_values, min_values = load_training_data('train_dyn.txt')

# Initialize XGBoost Classifier
xgb_classifier = XGBClassifier(
    max_depth=11, 
    n_estimators=750, 
    learning_rate=0.02, 
    random_state=42,  
    device='cuda',  # Ensure GPU usage if available
    missing=np.nan,  # Set the missing parameter to handle NaN values
    colsample_bytree=0.7,
    colsample_bylevel=0.8,
    colsample_bynode=0.9,
    gamma=1,
    min_child_weight=1,
    reg_alpha=0.5,
    reg_lambda=5,
    scale_pos_weight=2  # Adjust this based on the class imbalance in your dataset
)

# Fit the classifier
xgb_classifier.fit(X_train, y_train)

# Save max and min values to a file for future use
with open('max_min_values.txt', 'w') as f:
    f.write("Max Values:\n")
    np.savetxt(f, max_values, delimiter='\t')
    f.write("Min Values:\n")
    np.savetxt(f, min_values, delimiter='\t')

# Function to predict probabilities for test data
def predict_probabilities(input_file, output_file):
    # Load and process test data
    metadata, X_test = load_test_data(input_file, max_values, min_values)
    if X_test is None:
        print("Failed to load test data.")
        return

    # Predict probabilities
    y_probs = xgb_classifier.predict_proba(X_test)[:, 1]

    # Prepare output with metadata and predicted probabilities
    output_data = []
    for i in range(metadata.shape[0]):
        file_name, residue_index, WT_residue, mutant_residue = metadata[i]
        if WT_residue == mutant_residue:
            prediction = 0.0
        else:
            prediction = y_probs[i]
        output_data.append([file_name, residue_index, WT_residue, mutant_residue, prediction])

    # Write output to file
    with open(output_file, 'w') as f:
        for row in output_data:
            f.write('\t'.join(map(str, row)) + '\n')

# Predict probabilities for test data
predict_probabilities('saturation_mutagenesis_dyn_selected21.tsv', 'predictions.txt')

