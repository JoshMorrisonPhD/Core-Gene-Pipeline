# input is a blastp result file, either .txt or .html
# output is 3 .txt files
# generates recap .txt files of blastp runs:  genes that had hits,  genes without hits and a recap.txt showing numbers of each


import os
import argparse
from collections import defaultdict

def analyze_and_save_blast_results(file_path):
    gene_hits_info = defaultdict(list)
    gene_no_hits = set()
    base_name = os.path.splitext(os.path.basename(file_path))[0]

    with open(file_path, 'r') as file:
        lines = file.readlines()

    for line in lines:
        if line.startswith('# Query:'):
            current_gene = line.split('[locus_tag=')[-1].split(']')[0]
            gene_name = line.split('[gene=')[-1].split(']')[0] if '[gene=' in line else ""
            if gene_name == current_gene:  # No gene name found
                gene_name = ""  # Keep gene_name blank if not found
            gene_no_hits.add(current_gene)  # Assume no hits initially
        elif line.startswith('# 0 hits found'):
            continue  # Skip processing for genes with no hits
        elif not line.startswith('#') and line.strip():
            parts = line.strip().split('\t')
            if len(parts) > 1:
                subject_id = parts[1]
                gene_hits_info[current_gene].append((gene_name, current_gene, subject_id))
                if current_gene in gene_no_hits:
                    gene_no_hits.remove(current_gene)  # Remove from no hits if a hit is found

    output_dir = os.path.dirname(file_path)
    recap_path = os.path.join(output_dir, f"{base_name}_recap.txt")
    hits_path = os.path.join(output_dir, f"{base_name}_genes_with_hits.txt")
    zero_hits_path = os.path.join(output_dir, f"{base_name}_genes_with_zero_hits.txt")

    # Write detailed hits information
    with open(hits_path, 'w') as hits_file:
        for gene, hits in gene_hits_info.items():
            for hit in hits:
                gene_name, locus_tag, subject_id = hit
                hits_file.write(f"{gene_name}\t{locus_tag}\t{subject_id}\n")

    # Write genes with zero hits
    with open(zero_hits_path, 'w') as zero_hits_file:
        for gene in gene_no_hits:
            zero_hits_file.write(f"{gene}\n")

    # Write recap content and print summary to the command line
    recap_content = f"Total number of genes with one or more hits: {len(gene_hits_info)}\n" \
                    f"Total number of genes with zero hits: {len(gene_no_hits)}\n"
    with open(recap_path, 'w') as recap_file:
        recap_file.write(recap_content)
    
    print(f"Results for {base_name}:")
    print(recap_content)

def main():
    parser = argparse.ArgumentParser(description="Analyze BLAST results and generate summary files.")
    parser.add_argument("-i", "--input_folder", type=str, required=True, help="Input folder containing BLAST result files.")

    args = parser.parse_args()

    for file_name in os.listdir(args.input_folder):
        if file_name.endswith('.txt') and not any(substring in file_name for substring in ["_genes_with_hits", "_genes_with_zero_hits", "_recap"]):
            file_path = os.path.join(args.input_folder, file_name)
            analyze_and_save_blast_results(file_path)

if __name__ == "__main__":
    main()
