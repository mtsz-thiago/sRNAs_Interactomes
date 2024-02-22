# import
from jaal import Jaal
from jaal.datasets import load_got
import pandas as pd
import argparse


from jaal.datasets import load_got
# load the data
# edge_df, node_df = load_got()

edge_df = pd.read_csv("/workspaces/sRNAs_Interactomes/output/interactome_graphs/EP-OD0.5_w_alignments_edges.csv")
node_df = pd.read_csv("/workspaces/sRNAs_Interactomes/output/interactome_graphs/EP-OD0.5_w_alignments_nodes.csv")    

# create the app and server
app = Jaal(edge_df, node_df).create()
server = app.server