#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.output_dir = "$baseDir/output"
params.data_file = "$baseDir/data/Liu_sup5_data.xlsx"
params.cache_dir = "$baseDir/data"
params.queries_files_chunk_sizes = 20

params.salmonella_id = "GCF_000210855.2"

genomeFilename = "GCF_000210855.2_ASM21085v2_genomic.fna"
cdsFilename = "cds_from_genomic.fna"
params.blast_word_sizes = "7,11"

include { blast_wf as blastWFFullGenome } from "./module.blast" params(  queriesChunckSize: params.queries_files_chunk_sizes,
                                                    wordSizes_list: getWordSizesFromStringParam(params.blast_word_sizes))
include { blast_wf as blastWFCDS } from "./module.blast" params(  queriesChunckSize: params.queries_files_chunk_sizes,
                                                    wordSizes_list: getWordSizesFromStringParam(params.blast_word_sizes))

def getScenarioWordSizeKey(queryFilePath) {
    return queryFilePath.getName().split("\\.")[0].split("_")[0]
}

def getScenarionFromFilename(filename) {
    def matcher = (filename =~ /(.*)-(.*)$/)
    if (matcher.matches()) {
        return matcher[0][1]
    } else {
        return null
    }
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

process mergeResultsAndExpected {

    input:
    tuple val(scenarioAndWordSize), path('genome'), path('cds'), path('expected')

    output:
    path("${scenarioAndWordSize}_results.csv")

    //TODO this should be moved to a .py for clearity
    script:
    """
    #!/usr/bin/python3
    import pandas as pd

    full_genomeAlignments_df = pd.read_csv('genome', sep='\t')
    cds_alignments_df = pd.read_csv('cds', sep='\t')
    expected_results_df = pd.read_csv('expected', sep=',')

    merged_df = full_genomeAlignments_df.merge(cds_alignments_df, on='qseqid', how='outer', suffixes=('_genome','_cds'))
    merged_df = merged_df.merge(expected_results_df, left_on='qseqid', right_on='query_id', how='right', suffixes=('','_expected'))
    merged_df.to_csv('${scenarioAndWordSize}_results.csv', sep=',')
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
    expectedResults_ch = queriesAndExpected_ch.expected.flatten()
    
    queries_ch.collectFile(
        storeDir: "$params.output_dir/queries"
    )

    expectedResults_ch.collectFile(
        storeDir: "$params.output_dir/expected_alignments_results"
    )
    
    // run blast against full genome
    blastFullGenomeResults_ch = blastWFFullGenome(queries_ch, salmonellaGenome_ch)
    fullGenomeAlignments_ch = blastFullGenomeResults_ch.aligmentsResults_ch
    fullGenomeAlignments_ch.countLines().view(it -> "Number genomic alginments ${it}")
    fullGenomeAlignments_ch.collectFile(
        storeDir: "$params.output_dir/full_genome_alignments"
    )

    // run blast against cds
    blastResultsCDS_ch = blastWFCDS(queries_ch, salmonellaCDS_ch)
    cds_alignments_ch = blastResultsCDS_ch.aligmentsResults_ch
    cds_alignments_ch.countLines().view(it -> "Number CDS alginments ${it}")
    cds_alignments_ch.collectFile(
        storeDir: "$params.output_dir/cds_alignments"
    )

    // merge results
    keyFileGenomeAlignment_ch = fullGenomeAlignments_ch.map(it -> [getScenarioWordSizeKey(it), it])
    keyFileGenomeCDS_ch =  cds_alignments_ch.map(it -> [getScenarioWordSizeKey(it), it])
    
    keyFileExpectedResults_ch = expectedResults_ch.map(it -> [getScenarionFromFilename(it.getName()), it])

    filesToMerge = keyFileGenomeAlignment_ch
                        .join(keyFileGenomeCDS_ch)
                        .map(it -> [getScenarionFromFilename(it[0]), it[0], it[1], it[2]])
                        .combine(keyFileExpectedResults_ch, by:0)
                        .map(it -> [it[1],it[2],it[3],it[4]] )

    filesToMerge.count().view(it -> "Merging ${it} files")
    // filesToMerge.view()
    mergedResults_ch = mergeResultsAndExpected(filesToMerge)
    mergedResults_ch.collectFile(
        storeDir: "$params.output_dir"
    )
}

