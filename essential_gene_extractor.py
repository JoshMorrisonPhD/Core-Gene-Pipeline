import pandas as pd
import argparse
import os

def filter_locus_tags(input_file, output_dir):
    # Load the spreadsheet
    df = pd.read_excel(input_file)

    # Filter for locus tags with 0 test_num_insertions_mapped_per_feat and type CDS
    cds_essential_locus_tags = df[(df['test_num_insertions_mapped_per_feat'] == 0) & (df['type'] == 'CDS')]['locus_tag']

    # Filter for locus tags with 0 test_num_insertions_mapped_per_feat but type is not CDS
    non_cds_essential_locus_tags = df[(df['test_num_insertions_mapped_per_feat'] == 0) & (df['type'] != 'CDS')]['locus_tag']

    # Ensure the output directory exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Save the lists to separate .txt files
    cds_file_path = os.path.join(output_dir, 'cds_essential_locus_tags.txt')
    non_cds_file_path = os.path.join(output_dir, 'non_cds_essential_locus_tags.txt')

    # Write CDS essential locus tags to file and print count
    with open(cds_file_path, 'w') as file:
        file.write("\n".join(cds_essential_locus_tags))
    cds_count = len(cds_essential_locus_tags)

    # Write non-CDS essential locus tags to file and print count
    with open(non_cds_file_path, 'w') as file:
        file.write("\n".join(non_cds_essential_locus_tags))
    non_cds_count = len(non_cds_essential_locus_tags)

    # Print files saved, counts, and total results
    print(f"Files saved:\n{cds_file_path} with {cds_count} results\n{non_cds_file_path} with {non_cds_count} results")

    # Calculate and print the total unique results found
    total_unique_results = len(set(cds_essential_locus_tags).union(set(non_cds_essential_locus_tags)))
    print(f"Total unique results found: {total_unique_results}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter locus tags based on criteria.")
    parser.add_argument("-i", "--input", required=True, help="Input Excel file path.")
    parser.add_argument("-o", "--output", required=True, help="Output directory for the text files.")
    
    args = parser.parse_args()
    
    filter_locus_tags(args.input, args.output)
