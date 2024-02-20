#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.kmer_sz = 4

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
    path graphFile

    script:
    """
    #!/usr/bin/env python3

    import networkx as nx
    from Bio import pairwise2
    from Bio.Align import substitution_matrices

    substitution_matrix = substitution_matrices.load('BLOSUM62') 


    def add_alignmets_socres_to_edges(graph):
        for u, v, d in graph.edges(data=True):
            seq1 = graph.nodes[u]['seq']
            seq2 = graph.nodes[v]['seq']
            alignmentScore = pairwise2.align.globaldx(seq1, seq2, substitution_matrix, score_only=True)
            d['alignmentScore'] = alignmentScore
        
        return graph

    G = nx.read_gml('${graphFile}')
    G_w_features = add_alignmets_socres_to_edges(G)

    nx.write_gml(G_w_features, '$graphFile')
    """
}

workflow interactomeModeling_wf {

    take:
    chimeras_ch 

    main:
    graphs_gml_ch = dataToChimeraGraph(chimeras_ch) | addFeaturesToGraph

    emit:
    graphs_ch = graphs_gml_ch
}

