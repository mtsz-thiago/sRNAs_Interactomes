#!/usr/bin/env nextflow
params.output_dir = "$baseDir/output"
params.data_file = "$baseDir/data/Liu_sup5_data.xlsx"
params.cache_dir = "$baseDir/data"
params.queries_files_chunk_sizes = 10

params.salmonella_id = "GCF_000210855.2"

genomeFilename = "GCF_000210855.2_ASM21085v2_genomic.fna"
cdsFilename = "cds_from_genomic.fna"


// params.blastWordSZ_list = "4,7,11,15"
params.blast_word_sizes = "4,7,11,15"

def getScenarioFromFileName(queryFilePath) {
    return queryFilePath.name.split("\\.")[0].split("_")[0]
}

def getScenarionFromAlignedFIleNmeChunk(alginedFileName) {
    return "${alginedFileName}_alignments_results.tsv"
}

def mapAlignedTuplesToGroupChunks(it) {
    return tuple(getScenarionFromAlignedFIleNmeChunk(it[0]), it[1])
}

def getWordSizesFromStringParam(word_sizes_str) {
    return ((String)word_sizes_str).split(',').collect {it as Integer}
}

process downloadSalmonellaDataset {
    container 'biocontainers/ncbi-datasets-cli:15.12.0_cv23.1.0-4'
    publishDir params.cache_dir, mode: 'copy'

    input:
    val salmonellaRefSeqID

    output:
    path "salmonella_dataset.zip"

    script:
    """
    datasets download genome accession ${salmonellaRefSeqID} --filename salmonella_dataset.zip --include gff3,rna,cds,protein,genome,seq-report
    """

}

process createFastaFilesFromSupData {
    // publishDir params.cache_dir, mode: 'copy'

    input:
    path data_file

    output:
    tuple path("files/queries/*"), path("files/expected/*")

    script:
    """
    python3 $projectDir/src/sup_data_to_fasta.py -i ${data_file} -o files
    """

}

process makeBlastGenomeDB {
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

process extractGenomeDataFromZip {
    input:
    path genome_ziped_file

    output:
    path "ncbi_dataset/data/${params.salmonella_id}/*"

    script:
    """
    unzip ${genome_ziped_file}
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

    salmonellaZipedDataset_ch = downloadSalmonellaDataset( channel.of(params.salmonella_id) )
    salmonellaDataset_ch = extractGenomeDataFromZip(salmonellaZipedDataset_ch).flatten()
    
    salmonellaGenome_ch = salmonellaDataset_ch.filter { it -> it.getName() == genomeFilename }
    salmonellaCDS_ch = salmonellaDataset_ch.filter { it -> it.getName() == cdsFilename }

    salmonellaCDS_ch.count().view()

    // make blast db
    salmonella_db_ch = makeBlastDB(salmonellaCDS_ch)
    
    // create fasta files from sup data
    queries_and_expected_ch = createFastaFilesFromSupData(params.data_file).flatten()
    branched_queries_and_expected_ch = queries_and_expected_ch.branch { 
        queries:  it.getParent().getName() == "queries"
        expected:  it.getParent().getName() == "expected"
    }

    queries_ch = branched_queries_and_expected_ch.queries
    expected_results_ch = branched_queries_and_expected_ch.expected
    
    queries_ch.collectFile(
        storeDir: "$params.output_dir/queries"
    )

    expected_results_ch.collectFile(
        storeDir: "$params.output_dir/expected_alignments_results"
    )

    // load queries from fasta files
    splited_queries_ch = queries_ch
        .splitFasta(by: params.queries_files_chunk_sizes, file:true)
    
    // add word sizes to blast input channel
    blast_word_sz_list = getWordSizesFromStringParam(params.blast_word_sizes)
    blast_word_sz_ch = channel.from(blast_word_sz_list)

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

