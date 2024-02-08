#!/bin/bash

# Loop through each subdirectory
find . -type d | while read -r subfolder; do
  # Find a file that begins with "GCF_" in the current subfolder
  gcf_file=$(find "$subfolder" -maxdepth 1 -type f -name "GCF_*" -print -quit)

  # If a "GCF_" file exists in this subfolder
  if [[ -n $gcf_file ]]; then
    # Extract the full path without extension
    filename_no_ext="${gcf_file%.*}"

    # Check if the genomic.gff file exists in this subfolder
    if [[ -f "${subfolder}/genomic.gff" ]]; then
      # Rename the genomic.gff file to the new name based on the GCF_ file
      mv "${subfolder}/genomic.gff" "${filename_no_ext}.gff"
    fi
  fi
done

echo "All eligible genomic.gff files have been renamed."
