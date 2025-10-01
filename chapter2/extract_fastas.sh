#!/bin/bash

# this code is for extracting all the .fna files from a downloaded NCBI dataset to a new singular folder.

# Define the destination directory
DEST_DIR="Agal_full_genome_fastas"

# Create the destination directory if it doesn't exist
mkdir -p "$DEST_DIR"

# Find all .fna files that do not start with "cds" and copy them to the destination directory
find . -type f -name "*.fna" ! -name "cds*" -exec cp {} "$DEST_DIR" \;

echo "All eligible .fna files have been copied to $DEST_DIR."
