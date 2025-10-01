# inputs for this are a .txt line seperated keyfile of locus tags and a the ALLELE.fna output file from a PEPPAN run
# outputs are a multifasta with the sequences for each line in the keyfile that a match was found for
# this script was used to generate multiFASTA files for all the Core/Shell/Cloud genes for each species or intersection of species

from Bio import SeqIO
import re
import argparse

def load_keys(keyfile):
    """Load locus tags from the keyfile."""
    with open(keyfile, "r") as kf:
        return [line.strip() for line in kf if line.strip()]

def parse_fasta(fasta_file):
    """Parse the allele fasta file and return a dictionary of sequences by locus tag."""
    sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        header = record.description  # Full header
        parts = header.split(":")
        if len(parts) < 2:
            continue  # Skip malformed headers
        locus_tag = parts[1].split()[0]  # Extract locus tag
        sequences[header] = (locus_tag, record)
    return sequences

def filter_best_matches(keys, sequences):
    """Filter sequences to retain the best match based on criteria."""
    best_matches = {}
    keys_with_matches = set()
    
    for key in keys:
        matches = {hdr: seq for hdr, (locus, seq) in sequences.items() if locus.startswith(key)}
        
        if not matches:
            continue  # No match found
        
        keys_with_matches.add(key)
        # Sort matches by criteria: prioritize `_1` suffix and shortest header length
        sorted_matches = sorted(
            matches.items(),
            key=lambda x: (not re.search(r"_1(\s|$)", x[0]), len(x[0]))
        )
        
        best_header, best_seq = sorted_matches[0]
        best_matches[best_header] = best_seq
    
    return best_matches, keys_with_matches

def write_fasta(output_file, best_matches):
    """Write filtered sequences to a multifasta file."""
    with open(output_file, "w") as out_f:
        SeqIO.write(best_matches.values(), out_f, "fasta")

def main():
    parser = argparse.ArgumentParser(description="Filter FASTA sequences based on keyfile.")
    parser.add_argument("-k", "--keyfile", required=True, help="Path to keyfile.")
    parser.add_argument("-f", "--fasta", required=True, help="Path to input FASTA file.")
    parser.add_argument("-o", "--output", required=True, help="Path to output FASTA file.")
    args = parser.parse_args()
    
    keys = load_keys(args.keyfile)
    sequences = parse_fasta(args.fasta)
    best_matches, keys_with_matches = filter_best_matches(keys, sequences)
    write_fasta(args.output, best_matches)
    
    total_keys = len(keys)
    matched_keys = len(keys_with_matches)
    unmatched_keys = total_keys - matched_keys
    
    print(f"Filtered FASTA saved to {args.output} with {len(best_matches)} sequences.")
    print(f"Keys with matches: {matched_keys}")
    print(f"Keys without matches: {unmatched_keys}")

if __name__ == "__main__":
    main()

