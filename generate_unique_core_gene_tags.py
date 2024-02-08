import re
import os
from os import path
import sys
import fnmatch

# Function to check if the line is for reference_species and extract old locus tag
def extract_old_locus_tag_for_species(line, reference_strain_id):
    if not reference_strain_id:
        sys.exit("Must supply reference_strain_id")

    if line.startswith(reference_strain_id):
        #print(f'extract_old_locus_tag_for_species-> found line starting with reference_strain_id: {reference_strain_id}')
        old_tag_match = re.search(r"old_locus_tag=([^;\n:]+)", line)
        #if old_tag_match and old_tag_match.group(1):
        #    print(f'old_tag_match: {old_tag_match.group(1)}')
        #else:
        #    print("old_tag_match: None")
        return old_tag_match.group(1) if old_tag_match else None
    return None

species_folder_names = ["Agalactiae", "Equi", "Pneumo", "Suis", "Uberis"] # Add Iniae later, TODO: process "All" later
species_to_reference_strain_id = { "Agalactiae": "Agal_01173",
                    "Equi": "Equi_4047", 
                    "Pneumo": "Pneumo_TIGR4", 
                    "Suis": "Suis_P17", 
                    "Uberis": "Uberis_0140J"
                    }
                    #"Iniae: "STUFF",
                    #"All": "Agal_01173" }

for species in species_folder_names:
    print(f'Processing species {species} and mapping peppan and reference locus tags...')
    # input paths
    species_folder = path.join(path.curdir, species)
    #print(f'species_folder: {species_folder}')
    if not path.exists(species_folder):
        sys.exit(f'Can not find {species} species folder, script cannot run successfully')

    peppan_folder = path.join(species_folder, "annotated_genomes", "peppan_out")
    #print(f'peppan_folder: {peppan_folder}')
    
    if not path.exists(peppan_folder):
        sys.exit(f'Can not find the peppan output folder at {peppan_folder} for the species {species}')
    
    core_peppan_gene_locuses_file_name = None
    for file in os.listdir(peppan_folder):
        if fnmatch.fnmatch(file, '*_core_peppan_gene_locuses.txt'):
            core_peppan_gene_locuses_file_name = file
            break
    if not core_peppan_gene_locuses_file_name:
        sys.exit("Did you forget to run the R code to generate the output?")
    p#rint(f'core_peppan_gene_locuses_file_name: {core_peppan_gene_locuses_file_name}')

    core_peppan_gene_locuses_file_path = path.join(peppan_folder, core_peppan_gene_locuses_file_name)
    #print(f'core_peppan_gene_locuses_file_path: {core_peppan_gene_locuses_file_path}')

    # Path to the .gff file
    peppan_gff_file_path = path.join(peppan_folder, "PEPPAN.PEPPAN.gff")
    #print(f'peppan_gff_file_path: {peppan_gff_file_path}')

    # Paths for the output .txt files
    output_reference_tags_file_path = path.join(species_folder, 'core_' + species + '_locus_tags.txt')
    output_peppan_tags_not_found_file_path = path.join(species_folder, "peppan_locus_tags_not_found.txt")

    #print(f'output_reference_tags_file_path: {output_reference_tags_file_path}')
    #print(f'output_peppan_tags_not_found_file_path: {output_peppan_tags_not_found_file_path}')

    # Mapping of new locus tags to old locus tags for the specified species
    species_specific_mapping = {}
    with open(peppan_gff_file_path, 'r') as gff_file:
        for line in gff_file:
            old_tag = extract_old_locus_tag_for_species(line, species_to_reference_strain_id[species])
            new_tag_match = re.search(r"ortholog_group:[^:]+:([^:;\n]+)", line)  # Match the new locus tag
            if old_tag and new_tag_match:
                new_tag = new_tag_match.group(1)
                print(f'found new tag {new_tag} -> points to old tag {old_tag}')
                species_specific_mapping[new_tag] = old_tag
    print(f'Entries in species_specific_mapping (Full data dictionary): {len(species_specific_mapping)}')

    # Map the core genes (new locus tags) to old locus tags for the specified species
    with open(core_peppan_gene_locuses_file_path, 'r') as core_peppan_gene_locuses_file:
        file_content = [line.strip('"\n') for line in core_peppan_gene_locuses_file.readlines()]  # Strip quotes and newline characters
        core_genes_species_old_tags = [species_specific_mapping.get(gene) for gene in file_content if gene in species_specific_mapping]
        print(f'Entries in core_genes_species_old_tags: {len(core_genes_species_old_tags)}')

    # Filter out None values and create a unique set of old locus tags
    unique_core_genes_old_tags = set(tag for tag in core_genes_species_old_tags if tag)
    print(f'Entries in unique_core_genes_old_tags: {len(unique_core_genes_old_tags)}')

    # Count and a sample of unique old locus tags
    #print(f'Sample output length: {len(unique_core_genes_old_tags)}, sample: {list(unique_core_genes_old_tags)[:10]}')

    # Writing the old locus tags that are present to a file
    with open(output_reference_tags_file_path, 'w') as file:
        for tag in unique_core_genes_old_tags:
            file.write(f'{tag}\n')

    # Finding the new locus tags that couldn't be found and writing them to a file
    with open(core_peppan_gene_locuses_file_path, 'r') as core_peppan_gene_locuses_file:
        file_content = [line.strip('"\n') for line in core_peppan_gene_locuses_file.readlines()]  # Strip quotes and newline characters
        new_tags_not_found = [gene for gene in core_peppan_gene_locuses_file.readlines() if gene not in species_specific_mapping]
        with open(output_peppan_tags_not_found_file_path, 'w') as file:
            for tag in new_tags_not_found:
                file.write(f'{tag}\n')
