# move annotated .gff files from subdirectories to a single directory

import os
import shutil

# Define the top-level directory containing the subdirectories
top_level_dir = '/home/josh/Documents/EverySpecies'

# Define the subdirectories (you can also dynamically list these if you prefer)
subdirs = ['agal_annotated_genomes', 'equi_annotated_genomes', 'iniae_annotated_genomes',
           'pneumo_annotated_genomes', 'suis_annotated_genomes', 'uberis_annotated_genomes']

for subdir in subdirs:
    subdir_path = os.path.join(top_level_dir, subdir)
    species_name = subdir.split('_')[0]
    
    # List all files in the subdirectory
    for filename in os.listdir(subdir_path):
        file_path = os.path.join(subdir_path, filename)
        if os.path.isfile(file_path):
            # Create the new filename with the species prefix
            new_filename = f"{species_name}_{filename}"
            new_file_path = os.path.join(top_level_dir, new_filename)
            
            # Move the file to the top-level directory with the new name
            shutil.move(file_path, new_file_path)

print("Files have been successfully moved and renamed.")

