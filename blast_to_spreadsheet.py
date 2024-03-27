import pandas as pd
import sys
import argparse

fields = ['query id', 'subject id', 'alignment length', 'query length', 'subject length', 'q. start', 'q. end', 's. start', 's. end', 'evalue']

def get_args():
    try:
        parser = argparse.ArgumentParser(
            description="Generates an xlsx spreadsheet from a given blast result file"
        )
        parser.add_argument(
            "-f",
            "--file",
            action="store",
            help="The filename or filenames comma separated",
        )
        if len(sys.argv) == 1:
            parser.print_help(sys.stderr)
            sys.exit(1)

    except NameError:
        sys.stderr.write(
            "An exception occured with argument parsing. Check your provided options."
        )
        sys.exit(1)

    return parser.parse_args()


# Define a function to process BLAST results and create a DataFrame
def process_blast_results(file_path):
    with open(file_path, 'r') as file:
        blast_results = file.read()

    data = []
    for line in blast_results.split('\n'):
        if line.startswith('# Query:'):
            current_query = line.split()[2]
        elif line.startswith('# Database:'):
            current_database = line.split(': ')[1]
        elif line.startswith('# 0 hits found') or line.startswith('# 1 hits found'):
            hits_found = int(line.split()[1])
            if hits_found == 0:
                data.append({
                    'Query': current_query,
                    'Database': current_database,
                    'Hits found': hits_found,
                    **{field: '' for field in fields}
                })
        elif not line.startswith('#'):
            hit_values = line.split('\t')
            if hit_values and len(hit_values) > 1:
                hit_data = dict(zip(fields, hit_values))
                hit_data['Query'] = current_query
                hit_data['Database'] = current_database
                hit_data['Hits found'] = 1
                data.append(hit_data)

    return pd.DataFrame(data)

def main():
    args = get_args()
    file_paths = []
    if not args.file:
        sys.stderr.write("No file or comma separated list of files provided. Exiting.")
        sys.exit(1)
    else:
        file_paths = args.file.split(",")

    # Initialize dictionary to hold DataFrame and corresponding excel path for each file
    dfs_and_paths = {}

    # Process each file
    for file_path in file_paths:
        df = process_blast_results(file_path)
        excel_path = file_path.replace('.txt', '.xlsx')
        df.to_excel(excel_path, index=False)
        dfs_and_paths[file_path] = excel_path
        print(f'Wrote {excel_path} to disk...')

if __name__ == "__main__":
    main()