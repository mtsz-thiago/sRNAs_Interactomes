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

workflow interactomeModeling_wf {

    take:
    chimeras_ch 

    main:
    graphs_gml_ch = dataToChimeraGraph(chimeras_ch)

    emit:
    graphs_ch = graphs_gml_ch
}

