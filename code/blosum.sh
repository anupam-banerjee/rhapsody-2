#!/bin/bash

# Ensure two arguments are provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <amino_acid1> <amino_acid2>"
    exit 1
fi

# Assign arguments to variables
aa1=$1
aa2=$2

# BLOSUM62 matrix file
matrix_file="blosum62.txt"

# Function to get the score
blosum62_score() {
    awk -v aa1="$1" -v aa2="$2" '
    BEGIN {
        # Read the first line for amino acid labels
        getline < "'"$matrix_file"'"
        split($0, labels, " ")
        for (i in labels) {
            if (labels[i] == aa1) aa1_idx = i
            if (labels[i] == aa2) aa2_idx = i
        }
        if (!aa1_idx || !aa2_idx) {
            print -4
            exit
        }
    }
    NR > 1 {
        split($0, row, " ")
        if (row[1] == aa1) {
            print row[aa2_idx]
            exit
        }
    }
    ' $matrix_file
}

# Get the score
score=$(blosum62_score $aa1 $aa2)
#echo "The BLOSUM62 score for $aa1 and $aa2 is $score"
echo "$score"

