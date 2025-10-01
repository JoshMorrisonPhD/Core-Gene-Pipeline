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

# Loop through each permutation and run the PEPPAN command
for permutation in "${permutations[@]}"; do
    # Replace underscores with commas and add braces
    formatted_permutation="{${permutation//_/\,}}*.gff"
    peppan_prefix="-p $permutation"
    command="PEPPAN $peppan_prefix $formatted_permutation"
    
    # Print the command (or execute it if you prefer)
    echo "$command"
    
    # Execute the command and wait for it to finish
    eval $command
    
    # Check the exit status of the command
    if [ $? -ne 0 ]; then
        echo "Command failed: $command"
        exit 1
    fi
done

echo "All permutations have been processed."

