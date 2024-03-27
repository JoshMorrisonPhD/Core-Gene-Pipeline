import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import argparse

def get_args():
    try:
        parser = argparse.ArgumentParser(
            description="Generates an xlsx spreadsheet from a given blast result file"
        )
        parser.add_argument(
            "-s",
            "--species",
            action="store",
            required=True,
            help="The name of the species, e.g equi or agal",
        )
        parser.add_argument(
            "-f",
            "--file",
            action="store",
            required=True,
            help="The filename of the combined excel file",
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

def GenerateQueriesPerNumberOfDatabasesWithHits(species_name, combined_excel_file):
        # Load the spreadsheet
        data = pd.read_excel(combined_excel_file)

        # Count hits found for each query
        hits_count_per_query = data.groupby('Query')['Hits found'].sum()

        # Count how many queries have hits in 0 to 6 databases
        hits_distribution = hits_count_per_query.value_counts().sort_index()

        # Since we know there are 6 databases, ensure all possible counts (0-6) are represented in the final result
        hits_distribution = hits_distribution.reindex(range(0, 7), fill_value=0)

        # Display the hits_distribution
        print(hits_distribution)

        # Create a dictionary to hold lists of queries for each category (0-6 databases with hits)
        queries_by_hits_category = {i: [] for i in range(7)}

        # Populate the dictionary with queries
        for query, hits in hits_count_per_query.items():
                queries_by_hits_category[hits].append(query)

        # Save each list to a separate .txt file
        file_paths = {}
        for hits, queries in queries_by_hits_category.items():
                file_name = f'{species_name}_queries_with_{hits}_databases_hits.txt'
                file_paths[hits] = file_name
                with open(file_name, 'w') as file:
                        for query in queries:
                                file.write(query + '\n')


def GenerateNetworkGraphAndHeatmap(species_name, combined_excel_file):
    # Create a matrix indicating whether each query has a hit in each database
    # This involves transforming the original dataframe to have queries as rows and databases as columns

    df = pd.read_excel(combined_excel_file)
    # Create a matrix indicating whether each query has a hit in each database
    hits_matrix = pd.pivot_table(df, values='Hits found', index='Query', columns='Database', aggfunc=lambda x: int(any(x>0)), fill_value=0)

    # Use Pandas to_excel() function to save the DataFrame
    hits_matrix.to_excel(f'{species_name}_hits_matrix.xlsx', index=True)
    # print(f'Wrote the matrix hits_matrix.xlsx to disk...')
    print(f'Matrix saved to {combined_excel_file}')
       
    # Constructing edges based on shared queries between databases
    edges = []
    for i, db1 in enumerate(hits_matrix.columns):
        for j, db2 in enumerate(hits_matrix.columns):
            if j <= i:  # Avoid double counting and self-loops
                continue
            # Count the number of shared queries between db1 and db2
            shared_queries = (hits_matrix[db1] + hits_matrix[db2] == 2).sum()
            if shared_queries > 0:
                edges.append((db1, db2, shared_queries))

    # Create a graph from the edges
    G = nx.Graph()
    G.add_weighted_edges_from(edges)

    # Draw the network graph
    plt.figure(figsize=(10, 8))
    ax = plt.gca()  # Get the current axes
    pos = nx.spring_layout(G, seed=42)  # For consistent layout
    edges, weights = zip(*nx.get_edge_attributes(G, 'weight').items())
    nx.draw(G, pos, ax=ax, node_color='lightblue', edgelist=edges, edge_color=weights, width=4, edge_cmap=plt.cm.Blues,
        with_labels=True, font_weight='bold')
    plt.title('Network Graph of Shared Queries Between Databases')
    # Specify ax argument when creating the colorbar
    plt.colorbar(plt.cm.ScalarMappable(cmap=plt.cm.Blues), label='Number of Shared Queries', ax=ax)

    plt.savefig(f'{species_name}_shared_essential_core_networkgraph.png')
    print(f'Wrote shared_essential_core_networkgraph.png to disk...')

    # Creating the heatmap
    plt.figure(figsize=(40, 60))
    sns.heatmap(hits_matrix, cmap="YlGnBu", cbar_kws={'label': 'Hit Present (1) or Absent (0)'})
    plt.title('Heatmap of Query Hits Across Databases')
    plt.xlabel('Database')
    plt.ylabel('Query')
    plt.xticks(rotation=45)
    plt.savefig(f'{species_name}_shared_essential_core_heatmap.png', dpi=300)
    print(f'Wrote shared_essential_core_heatmap.png to disk...')

if __name__ == "__main__":
    args = get_args()
    GenerateQueriesPerNumberOfDatabasesWithHits(args.species, args.file)
    GenerateNetworkGraphAndHeatmap(args.species, args.file)