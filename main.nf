#!/usr/bin/env nextflow
 
params.interactome_files = [
    "$baseDir/data/EP.csv", 
    "$baseDir/data/ESP.csv",
    "$baseDir/data/SP.csv"]
params.output_dir = "$baseDir/output"


workflow {
    interactome_data_ch = Channel.fromPath(params.interactome_files)
        | splitCsv( sep: ',', header: true )
        // | view
        | map( row -> ">${row['RNA1 name']}\n${row['RNA1 seq']}\n>${row[1]}\n${row[18]}\n")
        
    interactome_data_ch.collectFile(name: "query.fa", storeDir: params.output_dir)
}