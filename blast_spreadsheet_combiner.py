import pandas as pd
import sys
import argparse

def get_args():
    try:
        parser = argparse.ArgumentParser(
            description="Generates an xlsx spreadsheet from a given blast result file"
        )
        parser.add_argument("-s", "--species", required=True, action="store", help="The species name to generate a master spreadsheet for")
        parser.add_argument("-o","--outfile", action="store", required=True, help="Output file to the master spreadsheet in.")
        if len(sys.argv) == 1:
            parser.print_help(sys.stderr)
            sys.exit(1)

    except NameError:
        sys.stderr.write(
            "An exception occured with argument parsing. Check your provided options."
        )
        sys.exit(1)

    return parser.parse_args()

def main():
    args = get_args()

    # Create an empty DataFrame with the required columns
    combined_df = pd.DataFrame(columns=[
        'Query', 'Database', 'Hits found', 'query id', 'subject id',
        'alignment length', 'query length', 'subject length',
        'q. start', 'q. end', 's. start', 's. end', 'evalue'
    ])

    species_name = args.species
    excel_paths = [
        f'{species_name}_essential_vs_iniaecoredb.xlsx',
        f'{species_name}_essential_vs_pneumocoredb.xlsx',
        f'{species_name}_essential_vs_suiscoredb.xlsx',
        f'{species_name}_essential_vs_uberiscoredb.xlsx',
        f'{species_name}_essential_vs_agalcoredb.xlsx',
        f'{species_name}_essential_vs_equicoredb.xlsx'
    ]

    # Process each file and concatenate the data into the combined DataFrame
    for excel_path in excel_paths:
        # Read the Excel file into a DataFrame
        temp_df = pd.read_excel(excel_path)
        
        # Determine the database name from the file name
        database_name = excel_path.split('_vs_')[1].replace('.xlsx', '')
        
        # Add the database name as a new column to the DataFrame
        temp_df['Database'] = database_name
        
        # Combine with the main DataFrame
        combined_df = pd.concat([combined_df if not combined_df.empty else None, temp_df])

    # Replace NaN values with empty strings
    combined_df.fillna('', inplace=True)

    # Reset index in the combined DataFrame
    combined_df.reset_index(drop=True, inplace=True)

    # Save the combined DataFrame to a new Excel file
    combined_excel_path = args.outfile
    combined_df.to_excel(combined_excel_path, index=False)

if __name__ == "__main__":
    main()