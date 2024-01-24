#!/usr/bin/env nextflow
params.output_dir = "$baseDir/output"
params.data_file = "$baseDir/data/Liu_sup5_data.xlsx"
params.cache_dir = "$baseDir/data"
params.queries_files_chunk_sizes = 10
params.salmonella_ref_genome_ftp_url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/945/GCF_000006945.2_ASM694v2/GCF_000006945.2_ASM694v2_genomic.fna.gz"

blast_word_sz_list = [7, 8, 9, 10, 11]

def getScenarioFromFileName(queryFilePath) {
    return queryFilePath.name.split("\\.")[0].split("_")[0]
}

def getScenarionFromAlignedFIleNmeChunk(alginedFileName) {
    return "${alginedFileName}_alignments_results.tsv"
}

def mapAlignedTuplesToGroupChunks(it) {
    return tuple(getScenarionFromAlignedFIleNmeChunk(it[0]), it[1])
}

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
    tuple path(query_file), path(db_file), val(blast_word_sz)

    output:
    tuple val(prefix), path("alignments_results.tsv")

    script:
    prefix = getScenarioFromFileName(query_file) + "w${blast_word_sz}"
    """
    blastn -query ${query_file} -word_size ${blast_word_sz} -db ${db_file}/salmonella_genome_db -out alignments_results.tsv -outfmt "6 qseqid qgi qacc qaccver qlen sseqid sallseqid sgi sallgi sacc saccver sallacc slen qstart qend sstart send qseq sseq evalue bitscore score length pident nident mismatch positive gapopen gaps ppos frames qframe sframe btop staxids sscinames scomnames sblastnames sskingdoms stitle salltitles sstrand qcovs qcovhsp"
    sed -i '1i qseqid\tqgi\tqacc\tqaccver\tqlen\tsseqid\tsallseqid\tsgi\tsallgi\tsacc\tsaccver\tsallacc\tslen\tqstart\tqend\tsstart\tsend\tqseq\tsseq\tevalue\tbitscore\tscore\tlength\tpident\tnident\tmismatch\tpositive\tgapopen\tgaps\tppos\tframes\tqframe\tsframe\tbtop\tstaxids\tsscinames\tscomnames\tsblastnames\tsskingdoms\tstitle\tsalltitles\tsstrand\tqcovs\tqcovhsp' alignments_results.tsv
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
    queries_ch = createFastaFilesFromSupData(params.data_file).flatten()
    
    // queries_ch.count().view(it -> "Number of queries: ${it}")
    queries_ch.collectFile(
        storeDir: "$params.output_dir/queries"
    )

    // load queries from fasta files
    splited_queries_ch = queries_ch
        .splitFasta(by: params.queries_files_chunk_sizes, file:true)

    splited_queries_ch.countFasta().view(c -> "Number of queries: ${c}")
    
    // add word sizes to blast input channel
    blast_word_sz_ch = channel.fromList(blast_word_sz_list)

    // combine queries and db before calling alignLocally
    blast_inputs_ch = splited_queries_ch
                        .combine(salmonella_db_ch)
                        .combine(blast_word_sz_ch)

    aligments_ch = alignLocally(blast_inputs_ch)
    
    // store output 
    aligments_ch.map( it -> mapAlignedTuplesToGroupChunks(it) )
                .collectFile(
                    keepHeader: true,
                    storeDir: params.output_dir
                )
}

