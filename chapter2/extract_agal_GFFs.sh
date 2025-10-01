#!/bin/bash

# this code is for extracting all the .gff files from the NCBI dataset to a new singular folder

# Define the destination directory
DEST_DIR="Agal_GFFs"

# Create the destination directory if it doesn't exist
mkdir -p "$DEST_DIR"

# Find all .gff files in subdirectories and copy them to the destination directory
find . -type f -name "*.gff" -exec cp {} "$DEST_DIR" \;

echo "All .gff files have been copied to $DEST_DIR."
