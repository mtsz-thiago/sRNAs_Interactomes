#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// GLOBAL PARAMS
genomeFilename = "GCF_000210855.2_ASM21085v2_genomic.fna"
cdsFilename = "cds_from_genomic.fna"
params.output_dir = "$baseDir/output"
params.data_file = "$baseDir/data/Liu_sup5_data.xlsx"
params.cache_dir = "$baseDir/data"
params.salmonella_id = "GCF_000210855.2"
// BLAST PARAMS
params.queriesChunckSize = 10
params.wordSizes_list = [11]
params.blastCullingLimit = 2
params.blastHSQSLimit = 50
params.minAlignmentCoverageThreshold = 0.9
params.minPidentThreshold = 0.9
// GRAPH PARAMS
params.kmer_sz = 4
params.neo4jURI = "bolt://neo4j:7687"
params.neo4jUser = "neo4j"
params.neo4jPassword = "Password"
params.neo4jDB = "neo4j"

include { blast_wf as blastWFFullGenome } from "./module/blast" params(
    queriesChunckSize: params.queriesChunckSize,
    wordSizes_list: params.wordSizes_list,
    blastCullingLimit: params.blastCullingLimit,
    blastHSQSLimit: params.blastHSQSLimit,
    minAlignmentCoverageThreshold: params.minAlignmentCoverageThreshold,
    minPidentThreshold: params.minPidentThreshold
)

include { blast_wf as blastWFCDS } from "./module/blast" params(
    queriesChunckSize: params.queriesChunckSize,
    wordSizes_list: params.wordSizes_list,
    blastCullingLimit: params.blastCullingLimit,
    blastHSQSLimit: params.blastHSQSLimit,
    minAlignmentCoverageThreshold: params.minAlignmentCoverageThreshold,
    minPidentThreshold: params.minPidentThreshold
)

include { interactomeModeling_wf as graphModelingWF } from "./module/graph" params(
    kmer_sz: params.kmer_sz
    neo4jURI: params.neo4jURI
    neo4jUser: params.neo4jUser
    neo4jPassword: params.neo4jPassword
    neo4jDB: params.neo4jDB
    )

def getScenarioWordSizeKey(queryFilePath) {
    return queryFilePath.getName().split("\\.")[0]
}

def getScenarionFromFilename(filename) {
    def matcher = (filename =~ /(.*)-(.*)$/)
    if (matcher.matches()) {
        return matcher[0][1]
    } else {
        return null
    }
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

process flattenDataFile {

    input:
    path dataFile

    output:
    path "*.csv"

    script:
    """
    #!/usr/bin/env python

    import pandas as pd
 
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
        
        w_query_id = flattened_df.assign(
            query_id = pd.factorize(flattened_df['seq'])[0]
        )

        return w_query_id
        
    sup_dict = pd.read_excel('${dataFile}', sheet_name=None)
    sup_dict = {k.replace('+', '_'): v for k, v in sup_dict.items() if k != 'Legend'}

    for k, raw in sup_dict.items():
        df = flatten_df(raw)
        df.to_csv(k + '.csv', index=False)

    """
}

process extractAlignmentQueries {
    input:
    path chimerasFile

    output:
    path outputFilename

    script:
    outputFilename = chimerasFile.baseName.replace(".csv", ".fasta")
    """
    #!/usr/bin/python3
    import pandas as pd
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    def map_to_fasta(df):
        return df.groupby('query_id').apply(map_record_to_SeqRecord).reset_index(drop=True)

    def map_record_to_SeqRecord(r):
        id = f"{r.name}"
        name = ''
        seq = Seq(r.iloc[0].seq)
        description = f"{id} from {r.chimera_idx.values}"
        return SeqRecord(
            id=id,
            description=description,
            name=name, 
            seq=seq)

    df = pd.read_csv('${chimerasFile}')
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

process mergeChimerasAndAlignments {

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
    rawData_ch = channel.of(params.data_file)
    chimeras_ch  = flattenDataFile(rawData_ch).flatten()
    chimeras_ch.collectFile(
        storeDir: "$params.output_dir/chimeras"
    )

    queries_ch = extractAlignmentQueries(chimeras_ch)
    queries_ch.collectFile(
        storeDir: "$params.output_dir/queries"
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

    graphModelingWF(chimeras_ch, fullGenomeAlignments_ch)
    
    // merge results
    keyFileGenomeAlignment_ch = fullGenomeAlignments_ch.map(it -> [getScenarioWordSizeKey(it), it])
    keyFileGenomeCDS_ch =  cds_alignments_ch.map(it -> [getScenarioWordSizeKey(it), it])
    
    keyFileChimeras_ch = chimeras_ch.map(it -> [getScenarionFromFilename(it.getName()), it])

    filesToMerge = keyFileGenomeAlignment_ch
                        .join(keyFileGenomeCDS_ch)
                        .map(it -> [getScenarionFromFilename(it[0]), it[0], it[1], it[2]])
                        .combine(keyFileChimeras_ch, by:0)
                        .map(it -> [it[1],it[2],it[3],it[4]] )

    filesToMerge.count().view(it -> "Merging ${it} files")
    // filesToMerge.view()
    mergedResults_ch = mergeChimerasAndAlignments(filesToMerge)
    mergedResults_ch.collectFile(
        storeDir: "$params.output_dir/alignments_results"
    )
}

