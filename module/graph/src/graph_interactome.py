import pandas as pd
from collections import Counter
import itertools
import pickle
import argparse
import networkx as nx
import os
import re

def create_chimera_graph(data_df, gname, nodes_properties, edges_properties):
    
    # Create an empty graph
    G = nx.Graph(name=gname)

    # Iterate over the rows of the dataframe
    for _, row in data_df.iterrows():
        query_id = row['queryId']
        chimera_idx = row['chimeraIdx']

        # Add the node to the graph
        G.add_node(query_id)
        
        # Add the properties to the node
        for att in nodes_properties:
            G.nodes[query_id][att] = row[att]

        # Find other nodes with the same chimera_idx
        nodes_with_same_chimera = data_df[data_df['chimeraIdx'] == chimera_idx]['queryId']

        edges_att = {}
        # Add edges between the current node and other nodes with the same chimera_idx
        for node in nodes_with_same_chimera:
            if node != query_id:
                G.add_edge(query_id, node)
                
                for att in edges_properties:    
                    G.edges[(node, query_id)][att] = row[att]
    
    return G

def convert_case(columns):
    new_columns = [re.sub(r'[_\s\'\-]([a-zA-Z])', lambda x: x.group(1).upper(), col) for col in columns]
    return new_columns

def load_sRNA_interactome_graph(data_csv_path, kmer_sz=4, output_file=None):

    # Extract the name from the basename of data_csv_path without extension
    name = os.path.splitext(os.path.basename(data_csv_path))[0]
    
    data_df = pd.read_csv(data_csv_path)
    data_df['origin'] = 'chimera'
    edges_properties = ['from', 'to', 'ligation from', 'ligation to', "Odds Ratio", "Fisher's exact test p-value", "Number of interactions", 'origin']
    nodes_properties = ['type', 'Strand', 'seq', 'origin']
    
    # data_df = to_numeric_data(data_df, nodes_properties, edges_properties)
    
    new_columns = convert_case(data_df.columns)
    data_df.rename(columns=dict(zip(data_df.columns, new_columns)), inplace=True)
    edges_properties = convert_case(edges_properties)
    nodes_properties = convert_case(nodes_properties)
    
    graph = create_chimera_graph(data_df, name, nodes_properties, edges_properties) 
    if output_file:
        nx.write_gml(graph, output_file)
        
    return graph

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Load sRNA interactome graph")
    parser.add_argument("--data_csv_path", type=str, help="Path to the data CSV file")
    parser.add_argument("--output_path", type=str, help="Path to the data CSV file")
    parser.add_argument("--kmer_sz", type=int, default=4, help="Size of k-mers (default: 4)")
    
    args = parser.parse_args()
    
    load_sRNA_interactome_graph(args.data_csv_path, args.kmer_sz, args.output_path)
