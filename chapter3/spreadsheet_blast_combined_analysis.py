# input is the species name + its $speciesname_combined.xlsx spreadsheet
# output is a $speciesname_hits_matrix.xlsx spreadsheet which we can manually parse for the queries with 0 hits to add to no_results_all_species.xlsx, as well as some visualisations (heatmap + network graph)
# this script was only used to be able to get the queries with 0 hits for no_results_all_species.xlsx

import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import argparse
from collections import defaultdict

def get_args():
    try:
        parser = argparse.ArgumentParser(
            description="Generates an xlsx spreadsheet and visualizations from a given BLAST result file"
        )
        parser.add_argument(
            "-s",
            "--species",
            action="store",
            required=True,
            help="The name of the species, e.g. equi or agal",
        )
        parser.add_argument(
            "-f",
            "--file",
            action="store",
            required=True,
            help="The filename of the combined Excel file",
        )
        if len(sys.argv) == 1:
            parser.print_help(sys.stderr)
            sys.exit(1)

    except NameError:
        sys.stderr.write("An exception occurred with argument parsing. Check your provided options.\n")
        sys.exit(1)

    return parser.parse_args()

def GenerateQueriesPerNumberOfDatabasesWithHits(species_name, combined_excel_file):
    # Load the spreadsheet
    data = pd.read_excel(combined_excel_file)

    # Count number of unique databases
    num_databases = data['Database'].nunique()
    print(f"Comparing against {num_databases} databases.")

    # Create a matrix: 1 if a query had any hit in a database, 0 otherwise
    hits_matrix = pd.pivot_table(
        data,
        values='Hits found',
        index='Query',
        columns='Database',
        aggfunc=lambda x: int(any(x > 0)),
        fill_value=0
    )

    # Count how many databases each query hit
    hits_count_per_query = hits_matrix.sum(axis=1)

    # Count how many queries have hits in 0 to N databases
    hits_distribution = hits_count_per_query.value_counts().sort_index()
    hits_distribution = hits_distribution.reindex(range(0, num_databases + 1), fill_value=0)
    print("Hits found")
    print(hits_distribution)

    # Save query lists grouped by number of databases with hits
    queries_by_hits_category = defaultdict(list)
    for query, hits in hits_count_per_query.items():
        queries_by_hits_category[hits].append(query)

    file_paths = {}
    for hits, queries in queries_by_hits_category.items():
        file_name = f'{species_name}_queries_with_{hits}_databases_hits.txt'
        file_paths[hits] = file_name
        with open(file_name, 'w') as file:
            for query in queries:
                file.write(query + '\n')

def GenerateNetworkGraphAndHeatmap(species_name, combined_excel_file):
    df = pd.read_excel(combined_excel_file)

    # Create binary matrix of hits
    hits_matrix = pd.pivot_table(
        df,
        values='Hits found',
        index='Query',
        columns='Database',
        aggfunc=lambda x: int(any(x > 0)),
        fill_value=0
    )

    # Save to Excel
    hits_matrix.to_excel(f'{species_name}_hits_matrix.xlsx', index=True)
    print(f'Matrix saved to {species_name}_hits_matrix.xlsx')

    # Construct graph edges based on shared queries
    edges = []
    for i, db1 in enumerate(hits_matrix.columns):
        for j, db2 in enumerate(hits_matrix.columns):
            if j <= i:
                continue
            shared_queries = (hits_matrix[db1] + hits_matrix[db2] == 2).sum()
            if shared_queries > 0:
                edges.append((db1, db2, shared_queries))

    # Create network graph
    G = nx.Graph()
    G.add_weighted_edges_from(edges)

    plt.figure(figsize=(10, 8))
    ax = plt.gca()
    pos = nx.spring_layout(G, seed=42)
    edges_drawn, weights = zip(*nx.get_edge_attributes(G, 'weight').items())
    nx.draw(
        G, pos, ax=ax,
        node_color='lightblue',
        edgelist=edges_drawn,
        edge_color=weights,
        width=4,
        edge_cmap=plt.cm.Blues,
        with_labels=True,
        font_weight='bold'
    )
    plt.title('Network Graph of Shared Queries Between Databases')
    plt.colorbar(
        plt.cm.ScalarMappable(cmap=plt.cm.Blues),
        label='Number of Shared Queries',
        ax=ax
    )
    plt.savefig(f'{species_name}_shared_essential_core_networkgraph.png')
    print(f'Wrote {species_name}_shared_essential_core_networkgraph.png to disk...')

    # Create heatmap
    plt.figure(figsize=(40, 60))
    sns.heatmap(hits_matrix, cmap="YlGnBu", cbar_kws={'label': 'Hit Present (1) or Absent (0)'})
    plt.title('Heatmap of Query Hits Across Databases')
    plt.xlabel('Database')
    plt.ylabel('Query')
    plt.xticks(rotation=45)
    plt.savefig(f'{species_name}_shared_essential_core_heatmap.png', dpi=300)
    print(f'Wrote {species_name}_shared_essential_core_heatmap.png to disk...')

if __name__ == "__main__":
    args = get_args()
    GenerateQueriesPerNumberOfDatabasesWithHits(args.species, args.file)
    GenerateNetworkGraphAndHeatmap(args.species, args.file)
