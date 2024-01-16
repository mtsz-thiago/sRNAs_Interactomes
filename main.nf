#!/usr/bin/env nextflow
 
params.queries_files = [
    "$baseDir/output/queries/EP_RNA1.fa", 
    "$baseDir/output/queries/EP_RNA2.fa", 
    "$baseDir/output/queries/ESP_RNA1.fa",
    "$baseDir/output/queries/ESP_RNA2.fa",
    "$baseDir/output/queries/SP_RNA1.fa",
    "$baseDir/output/queries/SP_RNA2.fa"]

params.output_dir = "$baseDir/output"
params.queries_files_chunk_sizes = 10

process align {
    container 'ncbi/blast'

    input:
    path query_file

    output:
    path "alignments_results.txt"

    script:
    """
    blastn -query ${query_file} -db nt -remote -out alignments_results.txt -outfmt 6
    """
}

workflow {
    queries_ch = Channel
                    .fromPath(params.queries_files)
                    .splitFasta(by: params.queries_files_chunk_sizes, file:true)
    
    aligments_ch = align(queries_ch)

    aligments_ch.collectFile(
        name: "all_results.txt", 
        keepHeader: true, storeDir: params.output_dir)
}