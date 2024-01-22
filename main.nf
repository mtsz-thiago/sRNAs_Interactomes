#!/usr/bin/env nextflow
params.output_dir = "$baseDir/output"
params.data_file = "$baseDir/data/Liu_sup5_data.xlsx"
params.cache_dir = "$baseDir/data"
params.queries_files_chunk_sizes = 10
params.salmonella_ref_genome_ftp_url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/945/GCF_000006945.2_ASM694v2/GCF_000006945.2_ASM694v2_genomic.fna.gz"

process createFastaFilesFromSupData {

    input:
    path data_file

    output:
    path "queries/*"

    script:
    """
    python3 $projectDir/src/sup_data_to_fasta.py -i ${data_file} -o queries
    """

}

process downloadSalmonellaGenome {

    publishDir params.cache_dir, mode: 'copy'

    output:
    path "salmonella_genome.fna.gz"

    script:
    """
    curl -o salmonella_genome.fna.gz ${params.salmonella_ref_genome_ftp_url}
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
    tuple path(query_file), path(db_file)

    output:
    tuple val(prefix), path("alignments_results.tsv")

    script:
    prefix = query_file.name.split("\\.")[0].split("_")[0]
    """
    blastn -query ${query_file} -db ${db_file}/salmonella_genome_db -out alignments_results.tsv -outfmt 6
    sed -i '1i Query_ID\tSubject_ID\tPIdentity\tAlignment_Length\tMismatches\tGap_Openings\tQuery_Start\tQuery_End\tSubject_Start\tSubject_End\tE_value\tBit_Score' alignments_results.tsv
    """
}

workflow {

    // download salmonella genome
    salmonella_gz_genome_ch = downloadSalmonellaGenome()
    
    salmonella_gz_genome_ch | collectFile(
            name: "salmonella_genome.fna.gz",
            storeDir: params.output_dir
        )

    salmonella_genome_ch = unzipGenome(salmonella_gz_genome_ch)

    // make blast db
    salmonella_db_ch = makeBlastDB(salmonella_genome_ch)

    // create fasta files from sup data
    queries_ch = createFastaFilesFromSupData(params.data_file) | flatten 
    
    // queries_ch.count().view(it -> "Number of queries: ${it}")
    queries_ch.collectFile(
        storeDir: "$params.output_dir/queries"
    )

    // load queries from fasta files
    splited_queries_ch = queries_ch
        .splitFasta(by: params.queries_files_chunk_sizes, file:true)

    splited_queries_ch.countFasta().view(c -> "Number of queries: ${c}")
    
    // combine queries and db before calling alignLocally
    query_db_ch = splited_queries_ch.combine(salmonella_db_ch)
    aligments_ch = alignLocally(query_db_ch)
    
    // store output 
    aligments_ch.map(it -> tuple("${it[0]}_alignments_results.tsv", it[1]))
                .collectFile(
                    keepHeader: true,
                    storeDir: params.output_dir
                )
}

