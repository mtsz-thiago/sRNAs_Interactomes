#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.kmer_sz = 4
params.neo4jURI = "bolt://neo4j:7687"
params.neo4jUser = "neo4j"
params.neo4jPassword = "Password"
params.neo4jDB = "neo4j"

chimerasNodesColumns = ['name','Strand','from','to','type','seq','query_id']
chimerasEdgesColumns = ['ligation from','ligation to','Number of interactions','Odds Ratio','fishersPValue']

def getKeyFromFilePath(filePath) {
    return filePath.getName()[0..3]
}

process dataToChimeraGraph {
    input:
    path dataFile

    output:
    path outputFilename

    script:
    outputFilename = "${dataFile.baseName}.gml"
    """
    python3 $projectDir/module/graph/src/graph_interactome.py --data_csv_path ${dataFile} --kmer_sz $params.kmer_sz --output_path ${outputFilename}
    """
}

process addFeaturesToGraph {
    input:
    path graphFile

    output:
    path outputFilename

    script:
    outputFilename = "${graphFile.baseName.split('_')[0]}_w_features.gml"
    """
    #!/usr/bin/env python3

    import networkx as nx
    from Bio import pairwise2
    from Bio.Align import substitution_matrices

    substitution_matrix = substitution_matrices.load('BLOSUM62') 

    def add_alignmets_scores_to_edges(graph):
        for u, v, d in graph.edges(data=True):
            seq1 = graph.nodes[u]['seq']
            seq2 = graph.nodes[v]['seq']
            alignmentScore = pairwise2.align.globaldx(seq1, seq2, substitution_matrix, score_only=True)
            d['alignmentScore'] = alignmentScore
        
        return graph

    G = nx.read_gml('$graphFile')
    G_w_features = add_alignmets_scores_to_edges(G)

    nx.write_gml(G_w_features, '${outputFilename}')
    """
}

process createNeo4jDB {

    output:
    val dbName

    script:
    """
    #!/usr/bin/env python3

    from neo4j import GraphDatabase

    NEO4J_AUTH = ('${params.neo4jUser}','${params.neo4jPassword}')
    driver = GraphDatabase.driver('${params.neo4jURI}', auth=NEO4J_AUTH)

    with driver.session() as session:
        session.run("CREATE DATABASE ${params.neo4jDB} IF NOT EXISTS")

    dbName = '${params.neo4jDB}'
    """
}


process loadChimerasToDB {

    input:
    path chimerasFile
    val dbName

    output:
    tuple val(graphName), val(numberOfNodes), val(numberOfEdges)

    script:
    graphName = chimerasFile.baseName
    numberOfNodes = 0
    numberOfEdges = 0
    """
    #!/usr/bin/env python3

    import pandas as pd
    from graphdatascience import GraphDataScience

    NEO4J_AUTH = ('${params.neo4jUser}','${params.neo4jPassword}')
    gds = GraphDataScience('${params.neo4jURI}', auth=NEO4J_AUTH, database='${dbName}')

    graphName = '${chimerasFile.baseName}'

    chimeras_df = pd.read_csv('${chimerasFile}')
    chimeras_df.rename(columns={"Fisher\'s exact test p-value": 'fishersPValue'}, inplace=True)
    node_columns = '${chimerasNodesColumns.join(', ')}'.split(', ')
    edge_columns = '${chimerasEdgesColumns.join(', ')}'.split(', ')
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
    """

}

process addAlignmentsToGraph {
    input:
    tuple path(graphFile), path(alignmentsFile)

    output:
    path outputFilename

    script:
    outputFilename = "${graphFile.baseName.split('_')[0]}_w_alignments.gml"
    """
    #!/usr/bin/env python3

    import pandas as pd
    import networkx as nx

    G = nx.read_gml('${graphFile}')
    alignments_df = pd.read_csv('${alignmentsFile}', sep='\t')
    alignments_df['origin'] = 'alignment'

    node_properties = [ 'sseq', 'sstrand', 'origin']
    edge_properties = [ 'sstart', 'send', 'qstart', 'qend', 'evalue', 'bitscore', 'score', 'length', 'pident', 'nident', 'mismatch', 'gapopen', 'gaps', 'ppos', 'origin']

    def add_alignments_to_graph(graph, alignments_df, node_properties, edge_properties):
        
        for i, row in alignments_df.iterrows():

            node_id = f"{row['qseqid']}_{i}"
            graph.add_node(node_id)

            for att in node_properties:
                graph.nodes[node_id][att] = row[att]

            graph.add_edge(row['qseqid'], node_id)

            for att in edge_properties:
                graph.edges[row['qseqid'], node_id][att] = row[att]

        return graph

    G_w_alignments = add_alignments_to_graph(G, alignments_df, node_properties, edge_properties)

    nx.write_gml(G_w_alignments, '${outputFilename}')
    """
}

process createNodeAndEdgesDataframes {
    input:
    path graphFile

    output:
    tuple path(nodesOutputFilename), path(edgesOutputFilename)

    script:
    nodesOutputFilename = "${graphFile.baseName}_nodes.csv"
    edgesOutputFilename = "${graphFile.baseName}_edges.csv"
    """
    #!/usr/bin/env python3

    import pandas as pd
    import networkx as nx

    G = nx.read_gml('${graphFile}')

    nodes_df = pd.DataFrame([{'id': n, **d} for n, d in G.nodes(data=True)])
    edges_df = pd.DataFrame([{**d, 'from': u, 'to': v} for u, v, d in G.edges(data=True)])

    nodes_df.to_csv("${nodesOutputFilename}")
    edges_df.to_csv("${edgesOutputFilename}")

    """
}

process loadToDB {

    input:
    path graphFile

    output:
    tuple val(graphName), val(numberOfNodes), val(numberOfEdges)

    script:
    graphName = graphFile.baseName
    numberOfNodes = 0
    numberOfEdges = 0
    """
    #!/usr/bin/env python3

    import networkx as nx
    from py2neo import Graph, Node, Relationship

    graphName = '${graphFile.baseName}'

    graph = Graph(uri = '$params.neo4jURI', user = '$params.neo4jUser', password = '$params.neo4jPassword')

    graph_data = nx.read_gml('${graphFile}')

    numberOfNodes = 0
    for node in graph_data.nodes(data=True) if node['origin'] == 'chimera':
        node_properties = node[1]
        node_properties['id'] = node[0]
        node_properties['graph'] =  
        node = Node(['SEQENCE', 'CHIMERA'], **node_properties)
        graph.create(node)
        numberOfNodes += 1

    numberOfEdges = 0
    for edge in graph_data.edges(data=True):
        edge_properties = edge[2]
        edge_properties['from'] = edge[0]
        edge_properties['to'] = edge[1]
        edge_properties['graph'] = graphName
        edge = Relationship(Node('CHIMERA', id = edge[0]), '', Node('CHIMERA', id = edge[1]), **edge_properties)
        graph.create(edge)
        numberOfEdges += 1

    """    

}

workflow interactomeModeling_wf {

    take:
    chimeras_ch
    alignments_ch

    main:

    // db_ch = createNeo4jDB()
    db_ch = Channel.of(params.neo4jDB)
    loadResults_ch = loadChimerasToDB(chimeras_ch, db_ch)



    graphs_gml_ch = dataToChimeraGraph(chimeras_ch) | addFeaturesToGraph

    keyValueGraphGML_ch = graphs_gml_ch.map(it -> [getKeyFromFilePath(it), it])
    keyValueAlignmests_ch = alignments_ch.map(it -> [getKeyFromFilePath(it), it])

    graphAndAlignments_ch = keyValueGraphGML_ch.join(keyValueAlignmests_ch).map(it -> [it[1], it[2]])
    alignmentGraphsGML_ch = addAlignmentsToGraph( graphAndAlignments_ch)

    graphsEdgesAndNodesDF_ch = createNodeAndEdgesDataframes(alignmentGraphsGML_ch)

    emit:
    graphsGML_ch = alignmentGraphsGML_ch
    graphsDF_ch = graphsEdgesAndNodesDF_ch

}

