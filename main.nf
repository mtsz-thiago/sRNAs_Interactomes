#!/usr/bin/env nextflow
 
params.queries_files = [
    // "$baseDir/output/queries/EP_RNA1.fa", 
    // "$baseDir/output/queries/EP_RNA2.fa", 
    "$baseDir/output/queries/ESP_RNA1.fa",
    "$baseDir/output/queries/ESP_RNA2.fa",
    "$baseDir/output/queries/SP_RNA1.fa",
    "$baseDir/output/queries/SP_RNA2.fa"]

params.output_dir = "$baseDir/output"

process align {
    container 'ncbi/blast'

    input:
    file query_file

    output:
    file "alignments_results.txt"

    script:
    """
    blastn -query ${query_file} -db nt -remote -out alignments_results.txt -outfmt 6
    """
}

workflow {
    queries_ch = Channel.fromPath(params.queries_files)
    
    aligments_ch = align(queries_ch)

    aligments_ch.view()
}