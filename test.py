#!/usr/bin/env python3

import pandas as pd
import numpy as np
from graphdatascience import GraphDataScience

NEO4J_AUTH = ('neo4j','Password')
gds = GraphDataScience('bolt://neo4j:7687', auth=NEO4J_AUTH, database='neo4j')

graphName = 'teste'

chimeras_df = pd.read_csv('/workspaces/sRNAs_Interactomes/output/expected_alignments_results/EP-OD0.5.csv')
chimeras_df.rename(columns={"Fisher\'s exact test p-value": 'fishersPValue'}, inplace=True)
node_columns = ['name','Strand','from','to','type','seq','query_id']
edge_columns = ['ligation from','ligation to','Number of interactions','Odds Ratio','fishersPValue']
nodes_data = chimeras_df[node_columns]
edges_data = chimeras_df[edge_columns]


nodes = nodes_data.assign(
    nodeId=nodes_data.index,
    labels=lambda x: [["SEQUENCE", "CHIMERA"]] * len(nodes_data)
)

def find_pair(i, r, chimeras_df):
    return chimeras_df.query(f"({r['chimera_idx']} == chimera_idx) and ({i} != index)").index[0]

targetNode = [find_pair(i, r, chimeras_df) for i,r in chimeras_df.iterrows() ]
edges = edges_data.assign(
    sourceNode=edges_data.index,
    targetNode=targetNode
)

relationships = edges_data.assign(sourceNodeId=edges.sourceNode, targetNodeId=edges.targetNode, relationshipType="LIGATES")
G = gds.graph.construct(graphName, nodes, relationships)