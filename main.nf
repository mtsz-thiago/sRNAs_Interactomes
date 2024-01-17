#!/usr/bin/env nextflow

params.queries_files = [
    "$baseDir/output/queries/EP_RNA1.fa",
    "$baseDir/output/queries/EP_RNA2.fa",
    "$baseDir/output/queries/ESP_RNA1.fa",
    "$baseDir/output/queries/ESP_RNA2.fa",
    "$baseDir/output/queries/SP_RNA1.fa",
    "$baseDir/output/queries/SP_RNA2.fa"]

params.data_dir = "$baseDir/data"
params.output_dir = "$baseDir/output"
params.queries_files_chunk_sizes = 10

params.species_on_dataset = ["Salmonella enterica"]

process downloadSalmonellaGenome {
    // container 'biocontainers/ncbi-datasets-cli:15.12.0_cv23.1.0-4'

    input:
    val species_on_dataset

    output:
    path "salmonella_genome.gzip"

    script:
    """
    curl -o salmonella_genome.gzip https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_000006945.2/download?include_annotation_type=GENOME_FASTA
    """
}

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

    Channel.of(params.species_on_dataset)
        | downloadSalmonellaGenome
        | collectFile(
            name: "salmonella_genome.gzip",
            storeDir: params.data_dir
        )

    // queries_ch = Channel
    //                 .fromPath(params.queries_files)
    //                 .splitFasta(by: params.queries_files_chunk_sizes, file:true)

    // aligments_ch = align(queries_ch)

    // aligments_ch.collectFile(
    //     name: "all_results.txt",
    //     keepHeader: true,
    //     storeDir: params.output_dir)
}