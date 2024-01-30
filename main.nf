#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.output_dir = "$baseDir/output"
params.data_file = "$baseDir/data/Liu_sup5_data.xlsx"
params.cache_dir = "$baseDir/data"
params.queries_files_chunk_sizes = 10

params.salmonella_id = "GCF_000210855.2"

genomeFilename = "GCF_000210855.2_ASM21085v2_genomic.fna"
cdsFilename = "cds_from_genomic.fna"
params.blast_word_sizes = "4,7,11,15"

include { blast_wf as blastWFFullGenome } from "./module.blast" params(  queriesChunckSize: params.queries_files_chunk_sizes,
                                                    wordSizes_list: getWordSizesFromStringParam(params.blast_word_sizes))
include { blast_wf as blastWFCDS } from "./module.blast" params(  queriesChunckSize: params.queries_files_chunk_sizes,
                                                    wordSizes_list: getWordSizesFromStringParam(params.blast_word_sizes))

def getScenarioFromFileName(queryFilePath) {
    return queryFilePath.name.split("\\.")[0].split("_")[0]
}

// def getScenarionFromAlignedFIleNmeChunk(alginedFileName) {
//     return "${alginedFileName}_alignments_results.tsv"
// }

// def mapAlignedTuplesToGroupChunks(it) {
//     return tuple(getScenarionFromAlignedFIleNmeChunk(it[0]), it[1])
// }

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
    path("files/queries/*"), emit: queries
    path("files/expected/*"), emit: expected

    script:
    """
    python3 $projectDir/src/sup_data_to_fasta.py -i ${data_file} -o files
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


workflow {

    salmonellaZipedDataset_ch = downloadSalmonellaDataset( channel.of(params.salmonella_id) )
    salmonellaDataset_ch = extractGenomeDataFromZip(salmonellaZipedDataset_ch).flatten()
    
    salmonellaGenome_ch = salmonellaDataset_ch.filter { it -> it.getName() == genomeFilename }
    salmonellaCDS_ch = salmonellaDataset_ch.filter { it -> it.getName() == cdsFilename }

    // create fasta files from sup data
    queriesAndExpected_ch  = createFastaFilesFromSupData(params.data_file)
    queries_ch = queriesAndExpected_ch.queries.flatten()
    expected_results_ch = queriesAndExpected_ch.expected.flatten()
    
    queries_ch.collectFile(
        storeDir: "$params.output_dir/queries"
    )

    expected_results_ch.collectFile(
        storeDir: "$params.output_dir/expected_alignments_results"
    )
    
    // run blast against full genome
    blastFullGenomeResults_ch = blastWFFullGenome(queries_ch, salmonellaGenome_ch)
    blastFullGenomeResults_ch.aligmentsResults_ch.collectFile(
        storeDir: "$params.output_dir/full_genome_alignments"
    )

    // run blast against cds
    blastResultsCDS_ch = blastWFCDS(queries_ch, salmonellaCDS_ch)
    blastResultsCDS_ch.aligmentsResults_ch.collectFile(
        storeDir: "$params.output_dir/cds_alignments"
    )
}

