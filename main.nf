#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.output_dir = "$baseDir/output"
params.data_file = "$baseDir/data/Liu_sup5_data.xlsx"
params.cache_dir = "$baseDir/data"
params.queries_files_chunk_sizes = 20

params.salmonella_id = "GCF_000210855.2"

genomeFilename = "GCF_000210855.2_ASM21085v2_genomic.fna"
cdsFilename = "cds_from_genomic.fna"
params.blast_word_sizes = "11"

include { blast_wf as blastWFFullGenome} from "./module/blast" params(  queriesChunckSize: params.queries_files_chunk_sizes,
                                                    wordSizes_list: getWordSizesFromStringParam(params.blast_word_sizes))
include { indexSubjectSequences as index1} from "./module/blast"
include { blast_wf as blastWFCDS } from "./module/blast" params(  queriesChunckSize: params.queries_files_chunk_sizes,
                                                    wordSizes_list: getWordSizesFromStringParam(params.blast_word_sizes))
include { indexSubjectSequences as index2} from "./module/blast"

include { interactomeModeling_wf as graphModelingWF } from "./module/graph" params(kmer_sz: 4)

def getWordSizesFromStringParam(word_sizes_str) {
    return ((String)word_sizes_str).split(',').collect {it as Integer}
}
def getScenarioWordSizeKey(queryFilePath) {
    return queryFilePath.getName().split("\\.")[0]
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

process fromXLSXtoCSV {
    // publishDir params.cache_dir, mode: 'copy'

    input:
    path dataFile

    output:
    path ('*.csv')

    script:
    """
    #!/usr/bin/python3
    import pandas as pd

    def load_xlsx(input_file):
        sup_dict = pd.read_excel(input_file, sheet_name=None)
        sup_dict = {k.replace('+', '_'): v for k, v in sup_dict.items() if k != 'Legend'}
        return sup_dict

    sup_dict = load_xlsx('${dataFile}')

    for k, df in sup_dict.items():
        df.to_csv(k + '.csv', index=False)
    """

}

process flattenChimeraData {

    input:
    path dataFile

    output:
    path outputFilename

    script:
    outputFilename = dataFile.baseName.replaceAll(".csv", "_flattened.csv")
    """
    #!/usr/bin/python3

    import pandas as pd
    import uuid

    def add_query_id_column(df):
        df['query_id'] = df.groupby('seq').name.transform(lambda x: uuid.uuid4().hex)
        return df

    def flatten_df(df):
        dfRNA1 = df[[c for c in df.columns if 'RNA2' not in c]]
        rna1columnsMap = {c:c.strip('RNA1').strip() for c in dfRNA1.columns if 'RNA1' in c}
        dfRNA1 = dfRNA1.rename(columns=rna1columnsMap)
        dfRNA1['origin'] = 'RNA1'
        dfRNA1['chimera_idx'] = dfRNA1.index
        
        dfRNA2 = df[[c for c in df.columns if 'RNA1' not in c]]
        rna2columnsMap = {c:c.strip('RNA2').strip() for c in dfRNA2.columns if 'RNA2' in c}
        dfRNA2 = dfRNA2.rename(columns=rna2columnsMap)
        dfRNA2['origin'] = 'RNA2'
        dfRNA2['chimera_idx'] = dfRNA2.index
        
        flattened_df = pd.concat([dfRNA1, dfRNA2], ignore_index=True)
        
        return flattened_df

    
    df = pd.read_csv('${dataFile}')
    flattened_df = flatten_df(df)
    with_id_df = add_query_id_column(flattened_df)

    with_id_df.to_csv('${outputFilename}', index=False)

    """

}

process createQueriesFromChimeraData {
    input:
    path dataFile

    output:
    path outputFilename

    script:
    outputFilename = dataFile.baseName.replaceAll(".csv", ".fasta")
    """
    #!/usr/bin/python3

    import pandas as pd
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    def map_to_fasta(df):
        bio_df = df.groupby('seq').apply(lambda x: map_record_to_SeqRecord(x), include_groups=False)
        return bio_df

    def map_record_to_SeqRecord(r):
        id = f"{r.query_id.iloc[0]}"
        # print(id)
        name = ''
        seq = Seq(r.name)
        description = f"{id} from {r.chimera_idx.values}"
        return SeqRecord(
            id=id,
            description=description,
            name=name, 
            seq=seq)

    df = pd.read_csv('${dataFile}')
    fasta_df = map_to_fasta(df)
    SeqIO.write(fasta_df, '${outputFilename}', "fasta")
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

process mergeChimeraAndAlignmentsResults {

    input:
    tuple path(chimeras), path(genome), path(cds)

    output:
    path outputFilename

    //TODO this should be moved to a .py for clearity
    script:
    outputFilename = "${genome}.csv"
    """
    #!/usr/bin/python3
    import pandas as pd

    full_genomeAlignments_df = pd.read_csv('${genome}', sep='\t')
    cds_alignments_df = pd.read_csv('${cds}', sep='\t')
    chimeras_results_df = pd.read_csv('${chimeras}', sep=',')

    merged_df = full_genomeAlignments_df.merge(cds_alignments_df, on='qseqid', how='outer', suffixes=('_genome','_cds'))
    merged_df = merged_df.merge(chimeras_results_df, left_on='qseqid', right_on='query_id', how='right', suffixes=('','_chimeras'))
    merged_df.to_csv('${outputFilename}', sep=',')
    """
}


workflow dataAnalysisWF {

    take:
    chimerasData_ch
    salmonellaGenome_ch
    salmonellaCDS_ch

    main:

    chimera_ch = flattenChimeraData(chimerasData_ch)
    chimera_ch.collectFile(
        storeDir: "$params.output_dir/chimeras"
    )

    queries_ch = createQueriesFromChimeraData(chimera_ch)
    queries_ch.collectFile(
        storeDir: "$params.output_dir/queries"
    )

    blastFullGenomeResults_ch = blastWFFullGenome(queries_ch, salmonellaGenome_ch)
    fullGenomeAlignments_ch = blastFullGenomeResults_ch.aligmentsResults_ch
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

    interactomeGraphs_ch = graphModelingWF(chimera_ch, fullGenomeAlignments_ch)

    //   keyFileGenomeCDS_ch =  cds_alignments_ch.map(it -> [getScenarioWordSizeKey(it), it])
    
    // keyFileExpectedResults_ch = expectedResults_ch.map(it -> [getScenarionFromFilename(it.getName()), it])

    //   filesToMerge = keyFileGenomeAlignment_ch
    //                     .join(keyFileGenomeCDS_ch)
    //                     .map(it -> [getScenarionFromFilename(it[0]), it[0], it[1], it[2]])
    //                     .combine(keyFileExpectedResults_ch, by:0)
    //                     .map(it -> [it[1],it[2],it[3],it[4]] )


    keyValueBlastFullGenomeResults_ch = blastFullGenomeResults_ch.map( it -> [ getScenarioWordSizeKey(it), it])
    keyValueBlastResultsCDS_ch = blastResultsCDS_ch.map( it -> [ getScenarioWordSizeKey(it), it])

    genomeAlingments_ch = keyValueBlastFullGenomeResults_ch.join(keyValueBlastResultsCDS_ch).map(it-> [it[1], it[2]])
    filesToMerge = chimera_ch
                        .combine(genomeAlingments_ch)

    
    filesToMerge.count().view(it -> "Merging ${it} files")
    merged_ch = mergeChimeraAndAlignmentsResults(filesToMerge)

    emit:
    alignments_ch = merged_ch 
    graphsGML_ch = interactomeGraphs_ch.graphsGML_ch
    graphsDF_ch = interactomeGraphs_ch.graphsDF_ch.flatten()
    

}

workflow {

    salmonellaZipedDataset_ch = downloadSalmonellaDataset( channel.of(params.salmonella_id) )
    salmonellaDataset_ch = extractGenomeDataFromZip(salmonellaZipedDataset_ch).flatten()
    
    salmonellaGenome_ch = salmonellaDataset_ch.filter { it -> it.getName() == genomeFilename } | index1
    salmonellaCDS_ch = salmonellaDataset_ch.filter { it -> it.getName() == cdsFilename } | index2

    chimerasData_ch = fromXLSXtoCSV(params.data_file).flatten()

    analysisInputs_ch = chimerasData_ch
                            .combine(salmonellaGenome_ch)
                            .combine(salmonellaCDS_ch)

    analysisInputs_ch.count().view(it -> "Number of analysis inputs ${it}")
    chiumeraInput_ch = analysisInputs_ch.map(it -> it[0])
    genomeInput_ch = analysisInputs_ch.map(it -> it[1])
    cdsInput_ch = analysisInputs_ch.map(it -> it[2])

    analysisResults_ch = dataAnalysisWF(chiumeraInput_ch, genomeInput_ch, cdsInput_ch)

    analysisResults_ch.alignments_ch.collectFile(
        storeDir: "$params.output_dir/alignments"
    )
    analysisResults_ch.graphsGML_ch.collectFile(
        storeDir: "$params.output_dir/graphs"
    )
    analysisResults_ch.graphsDF_ch.collectFile(
        storeDir: "$params.output_dir/graphs"
    )

}

