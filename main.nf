#!/usr/bin/env nextflow
 
params.suplementary_data_urls = "$baseDir/data/sup_data_urls.csv"
params.output_dir = "$baseDir/output"

/*
 * Split a fasta file into multiple files
 */
process download_data_using_wget {
    
    input:
    tuple val(name), val(url)
    
    output:
    tuple val(name), file("sup_file")

    
    script:
    """
    echo "Downloading ${name} from ${url}"
    wget -O sup_file ${url}
    """
}
 
/*
 * Define the workflow
 */
workflow {
    suplementary_files_ch = Channel.fromPath(params.suplementary_data_urls) \
        | splitCsv(header: true, sep: ',') \
        | map { row-> tuple(row.name, row.url)} 
        | download_data_using_wget
        | collectFile(storeDir: params.output_dir )
}