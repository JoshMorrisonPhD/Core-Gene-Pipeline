# input is a nucleotide fasta or multifasta file .fna
# output is a n amino acid fasta file .faa
# this script translates nucleotide sequences to amino acid sequences

import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def translate_dna_to_protein(input_file, output_file):
    # List to hold the translated protein sequences
    protein_sequences = []

    # Iterate over each DNA sequence in the input file
    for dna_record in SeqIO.parse(input_file, "fasta"):
        # Translate the DNA sequence to a protein sequence
        protein_seq = dna_record.seq.translate(to_stop=True)
        # Create a new SeqRecord for the protein sequence
        protein_record = SeqRecord(protein_seq, id=dna_record.id, description="translated protein")
        # Append the protein record to the list
        protein_sequences.append(protein_record)

    # Write the protein sequences to the output file
    SeqIO.write(protein_sequences, output_file, "fasta")
    print(f"Protein sequences have been written to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Translate DNA sequences in a multi-FASTA file to protein sequences.")
    parser.add_argument("-i", "--input", required=True, help="Path to the input multi-FASTA file containing DNA sequences.")
    parser.add_argument("-o", "--output", required=True, help="Path to the output FASTA file to save translated protein sequences.")

    args = parser.parse_args()

    translate_dna_to_protein(args.input, args.output)

