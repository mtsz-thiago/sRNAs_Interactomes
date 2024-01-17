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
    path "salmonella_genome.fna.gz"

    script:
    """
    curl -o salmonella_genome.fna.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/945/GCF_000006945.2_ASM694v2/GCF_000006945.2_ASM694v2_genomic.fna.gz
    """
}

process makeBlastDB {
    container 'ncbi/blast'

    input:
    path genome_file

    output:
    path "SalmonellaDB"

    script:
    """
    mkdir -p 'SalmonellaDB'
    cd 'SalmonellaDB' && makeblastdb -in ../${genome_file} -dbtype nucl -out salmonella_genome_db -title 'SalmonellaDB'
    """
}

process unzipGenome {
    input:
    path genome_gz_file

    output:
    path "salmonella_genome.fna"

    script:
    """
    gzip -fd ${genome_gz_file}
    """
}

process alignLocally {
    container 'ncbi/blast'

    input:
    path query_file
    path db_file

    output:
    path "alignments_results.tsv"

    script:
    """
    blastn -query ${query_file} -db ${db_file}/salmonella_genome_db -out alignments_results.tsv -outfmt 6
    """
}

workflow {

    salmonella_gz_genome_ch = Channel.of(params.species_on_dataset)
        | downloadSalmonellaGenome
    
    salmonella_gz_genome_ch | collectFile(
            name: "salmonella_genome.fna.gz",
            storeDir: params.data_dir
        )
    
    salmonella_genome_ch = unzipGenome(salmonella_gz_genome_ch)

    salmonella_db_ch = makeBlastDB(salmonella_genome_ch)

    queries_ch = Channel
                    .fromPath(params.queries_files)
                    .splitFasta(by: params.queries_files_chunk_sizes, file:true)

    aligments_ch = alignLocally(queries_ch, salmonella_db_ch)
    
    aligments_ch.collectFile(
        name: "all_results.txt",
        keepHeader: true,
        storeDir: params.output_dir)
}