import pandas as pd
import argparse

def process_spreadsheet(input_file, output_file):
    # Load the spreadsheet
    data = pd.read_excel(input_file)

    # Extract Locus_tag values where test_num_insertions_mapped_per_feat is 0
    locus_tags_with_zero = data[data['test_num_insertions_mapped_per_feat'] == 0]['locus_tag']

    # Convert the series to a list and then to a line-separated string
    locus_tags_with_zero_str = "\n".join(locus_tags_with_zero.tolist())

    # Save the complete list of Locus_tag values to a text file
    with open(output_file, 'w') as file:
        file.write(locus_tags_with_zero_str)
    print(f"Output saved to {output_file}")

if __name__ == "__main__":
    # Initialize the argument parser
    parser = argparse.ArgumentParser(description="Extract Locus_tag values with 0 insertions from a spreadsheet.")
    
    # Add arguments for input file path and output file name
    parser.add_argument("input_file", help="Path to the input Excel file")
    parser.add_argument("output_file", help="Path for the output text file")

    # Parse the command line arguments
    args = parser.parse_args()

    # Process the spreadsheet with the provided arguments
    process_spreadsheet(args.input_file, args.output_file)

