#!/bin/bash

#This code is for appending refseq .gff files to the start of .fna files, to immitate prokka .gff formatting.

# Define the source and target directories
FNA_DIR="/home/josh/Documents/Agalactiae/Agal_full_genome_fastas"  # Replace with the path to your .fna files
GFF_DIR="/home/josh/Documents/Agalactiae/Agal_GFFs"  # Replace with the path to your .gff files
OUTPUT_DIR="/home/josh/Documents/Agalactiae/annotated_genomes"  # Replace with the path for the new combined files

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Check if the FNA_DIR exists
if [ ! -d "$FNA_DIR" ]; then
  echo "Error: Source directory for .fna files does not exist."
  exit 1
fi

# Check if the GFF_DIR exists
if [ ! -d "$GFF_DIR" ]; then
  echo "Error: Source directory for .gff files does not exist."
  exit 1
fi

# Loop through all .fna files in the source directory
for fna_file in "$FNA_DIR"/*.fna; do
  # Extract the base file name without the extension
  base_name=$(basename "$fna_file" .fna)
  
  # Define the path for the corresponding .gff file
  gff_file="$GFF_DIR/${base_name}.gff"

  # Define the path for the new output file
  output_file="$OUTPUT_DIR/${base_name}.gff"

  # Check if the corresponding .gff file exists
  if [ -f "$gff_file" ]; then
    # Combine the .gff and .fna files, and save the result to the output file
    cat "$gff_file" "$fna_file" > "$output_file"
  else
    echo "Warning: No corresponding .gff file found for $fna_file"
  fi
done

echo "All matching .fna files have been appended to their corresponding .gff files in $OUTPUT_DIR."

