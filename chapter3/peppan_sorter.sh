#!/bin/bash

# Define the list of permutations
permutations=(
    'agal_pneumo' 'agal_iniae' 'agal_uberis' 'agal_equi' 'agal_suis'
    'pneumo_iniae' 'pneumo_uberis' 'pneumo_equi' 'pneumo_suis'
    'iniae_uberis' 'iniae_equi' 'iniae_suis' 'uberis_equi' 'uberis_suis'
    'equi_suis' 'agal_pneumo_iniae' 'agal_pneumo_uberis' 'agal_pneumo_equi'
    'agal_pneumo_suis' 'agal_iniae_uberis' 'agal_iniae_equi' 'agal_iniae_suis'
    'agal_uberis_equi' 'agal_uberis_suis' 'agal_equi_suis' 'pneumo_iniae_uberis'
    'pneumo_iniae_equi' 'pneumo_iniae_suis' 'pneumo_uberis_equi'
    'pneumo_uberis_suis' 'pneumo_equi_suis' 'iniae_uberis_equi'
    'iniae_uberis_suis' 'iniae_equi_suis' 'uberis_equi_suis'
    'agal_pneumo_iniae_uberis' 'agal_pneumo_iniae_equi' 'agal_pneumo_iniae_suis'
    'agal_pneumo_uberis_equi' 'agal_pneumo_uberis_suis' 'agal_pneumo_equi_suis'
    'agal_iniae_uberis_equi' 'agal_iniae_uberis_suis' 'agal_iniae_equi_suis'
    'agal_uberis_equi_suis' 'pneumo_iniae_uberis_equi' 'pneumo_iniae_uberis_suis'
    'pneumo_iniae_equi_suis' 'pneumo_uberis_equi_suis' 'iniae_uberis_equi_suis'
    'agal_pneumo_iniae_uberis_equi' 'agal_pneumo_iniae_uberis_suis'
    'agal_pneumo_iniae_equi_suis' 'agal_pneumo_uberis_equi_suis'
    'agal_iniae_uberis_equi_suis' 'pneumo_iniae_uberis_equi_suis'
)

# Loop through each permutation
for permutation in "${permutations[@]}"; do
    # Create a directory for the permutation
    mkdir -p "$permutation"
    
    # Move files with the exact prefix into the corresponding folder
    for file in ${permutation}*; do
        # Check if file starts with the exact permutation and is not a directory
        if [[ -f "$file" && "$file" == $permutation* && "$file" != *"$permutation"_* ]]; then
            mv "$file" "$permutation/"
        fi
    done
done

echo "Files have been successfully organized into their respective folders."

