# The inputs for this are the .txt file outputs from gene_categoriser.R They need to be moved into the species folders that are within the present working directory and renamed "$speciesname_core_peppan_gene_locuses.txt".
# This script also supposes that within each species folder is a folder called "annotated_genomes" containing the "peppan_out" outputs folder. This script needs to be run from within the directory above all the species name directories:
# --PWD
# -----$Speciesname (contains the $speciesname_core_peppan_gene_locuses.txt file)
# -------annotated_genomes
# ---------peppan_out (contains the PEPPAN.PEPPAN.GFF file)
# The output of this script is a list of genes that can be matched to a reference genome multifasta to extract a multifasta file of core genes.
#This script was only used for the single species core genes, as it requires a reference genome to run on.

import re
import os
from os import path
import sys
import fnmatch
import collections

# Function to check if the line is for reference_species and extract old locus tag
def extract_old_locus_tag_for_species(line):
    old_tag_match = re.search(r"old_locus_tag=([^;\n:]+)", line)
    return old_tag_match.group(1) if old_tag_match else None

species_folder_names = ["Agalactiae", "Iniae", "All", "Equi", "Pneumo", "Suis", "Uberis"]
species_to_reference_strain_id = { "Agalactiae": "Agal_01173",
                    "Equi": "Equi_4047", 
                    "Iniae": "Iniae_GCF_000300915", 
                    "Pneumo": "Pneumo_TIGR4", 
                    "Suis": "Suis_P17", 
                    "Uberis": "Uberis_0140J",
                    "All": "Uberis_0140J"
                    }

for species in species_folder_names:
    print(f'Processing species {species} and mapping peppan and reference locus tags...')
    # input paths
    species_folder = path.join(path.curdir, species)
    if not path.exists(species_folder):
        sys.exit(f'Can not find {species} species folder, script cannot run successfully')

    peppan_folder = path.join(species_folder, "annotated_genomes", "peppan_out")
    
    if not path.exists(peppan_folder):
        sys.exit(f'Can not find the peppan output folder at {peppan_folder} for the species {species}')
    
    core_peppan_gene_locuses_file_name = None
    for file in os.listdir(species_folder):
        if fnmatch.fnmatch(file, '*_core_peppan_gene_locuses.txt'):
            core_peppan_gene_locuses_file_name = file
            break
    if not core_peppan_gene_locuses_file_name:
        sys.exit("Did you forget to run the R code to generate the output?")

    core_peppan_gene_locuses_file_path = path.join(species_folder, core_peppan_gene_locuses_file_name)

    # Path to the .gff file
    peppan_gff_file_path = path.join(peppan_folder, "PEPPAN.PEPPAN.gff")

    # Paths for the output .txt files
    output_reference_tags_file_path = path.join(species_folder, 'core_' + species + '_locus_tags.txt')
    output_peppan_tags_not_found_file_path = path.join(species_folder, "peppan_locus_tags_not_found.txt")
    output_duplicate_tags_file_path = path.join(species_folder, "duplicate_tags_found.txt")


    # Mapping of new locus tags to old locus tags for the specified species
    species_specific_mapping = {}
    with open(peppan_gff_file_path, 'r') as gff_file:
        for line in gff_file:
            # Only pick lines starting with reference_strain_id (e.g Equi_4047)
            if line.startswith(species_to_reference_strain_id[species]):
                # Does the line have a peppan locustag?

                ortholog_group_content = re.search(r"ortholog_group:([^;]+)", line)
                if ortholog_group_content:
                    # Step 2: Split the content by ',GCF_' to handle multiple entries
                    # Prepend 'GCF_' to each split part except the first one to restore the cut-off part
                    entries = ['GCF_' + entry if i != 0 else entry for i, entry in enumerate(ortholog_group_content.group(1).split(',GCF_'))]

                    # Step 3: Extract the "RS tag" from each segment
                    new_tag_matches = [re.search(r":([^:]+):", entry).group(1) for entry in entries if re.search(r":([^:]+):", entry)]

                    # Create an entry in the species specific lookup table if we have a peppan locus tag (old_locus_tag can either exist or not)
                    for new_tag_match in new_tag_matches:
                        old_tag = extract_old_locus_tag_for_species(line)
                        new_tag = new_tag_match
                        species_specific_mapping[new_tag] = old_tag

    print(f'Entries in species_specific_mapping (Full data dictionary): {len(species_specific_mapping)}')

    # Map the core genes (new locus tags) to old locus tags for the specified species
    with open(core_peppan_gene_locuses_file_path, 'r') as core_peppan_gene_locuses_file:
         # Strip quotes and newline characters
        file_content = [line.strip('"\n') for line in core_peppan_gene_locuses_file.readlines()]

        # List of dictionary values (old_locus_tag) for peppan locus tag keys found in dictionary
        core_genes_species_old_tags = [species_specific_mapping.get(locus_tag) for locus_tag in file_content if locus_tag in species_specific_mapping]

        # List of dictionary values (old_locus_tag) without None values
        core_genes_old_tags_values = [tag for tag in core_genes_species_old_tags if tag]

        duplicate_core_genes_old_tags_values = [duplicate for duplicate, count in collections.Counter(core_genes_old_tags_values).items() if count > 1]

        duplicate_core_genes_old_tags = {}
        for duplicate_locus_tag in duplicate_core_genes_old_tags_values:
            keys = [key for key, v in species_specific_mapping.items() if v == duplicate_locus_tag]
            duplicate_core_genes_old_tags[duplicate_locus_tag] = keys
            
        unique_core_genes_old_tags = set(core_genes_old_tags_values)

        # List of peppan locus tag keys for which dictionary has no value for that key
        new_tags_not_found = [gene for gene in file_content if species_specific_mapping.get(gene) is None]
            
    # Count and a sample of unique old locus tags
    print(f'Number of lines in core_peppan_gene_locuses_txt: {len(file_content)}')
    print(f'Entries in core_genes_species_old_tags: {len(core_genes_species_old_tags)}')
    print(f'Number of duplicate old_locus_tags: {len(duplicate_core_genes_old_tags_values)}')
    print(f'Entries in new_tags_not_found: {len(new_tags_not_found)}')
    print(f'Entries in unique_core_genes_old_tags: {len(unique_core_genes_old_tags)}')
    
    # Writing the old locus tags that are present to a file
    with open(output_reference_tags_file_path, 'w') as file:
        for tag in unique_core_genes_old_tags:
            file.write(f'{tag}\n')

    with open(output_peppan_tags_not_found_file_path, 'w') as file:
        for tag in new_tags_not_found:
            file.write(f'{tag}\n')

    with open(output_duplicate_tags_file_path, 'w') as file:
        for seq_value, duplicate_keys in duplicate_core_genes_old_tags.items():
            file.write(f'{seq_value}\t')
            for duplicate_key in duplicate_keys:
                file.write(f'{duplicate_key}\t')
            file.write("\n")
    print("\n")
