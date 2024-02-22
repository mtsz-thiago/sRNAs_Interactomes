#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.kmer_sz = 4


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

workflow interactomeModeling_wf {

    take:
    chimeras_ch
    alignments_ch

    main:
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

