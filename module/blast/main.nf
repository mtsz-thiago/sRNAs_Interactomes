#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.queriesChunckSize = 10
params.wordSizes_list = [4,7,11,15]
params.blastCullingLimit = 2
params.blastHSQSLimit = 50
params.minAlignmentCoverage = 0.9

blastOuputColumns = ['qseqid', 'qgi', 'qacc', 'qaccver', 'qlen', 'sseqid', 'sallseqid', 'sgi', 'sallgi', 'sacc', 'saccver', 'sallacc', 'slen', 'qstart', 'qend', 'sstart', 'send', 'qseq', 'sseq', 'evalue', 'bitscore', 'score', 'length', 'pident', 'nident', 'mismatch', 'positive', 'gapopen', 'gaps', 'ppos', 'frames', 'qframe', 'sframe', 'btop', 'staxids', 'sscinames', 'scomnames', 'sblastnames', 'sskingdoms', 'stitle', 'salltitles', 'sstrand', 'qcovs', 'qcovhsp']

def getScenarioFromFileName(queryFilePath) {
    return queryFilePath.name.split("\\.")[0].split("_")[0]
}

process indexSubjectSequences {
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

process alignLocally {
    container 'ncbi/blast'

    input:
    tuple path(query_file), path(db_file), val(blast_word_sz)

    output:
    tuple val(prefix), path(outputFilename)

    script:
    prefix = getScenarioFromFileName(query_file) + "w${blast_word_sz}"
    queryBaseName = query_file.baseName
    outputFilename = "${queryBaseName}_alignments_results.tsv"
    """
    blastn -query ${query_file} -word_size ${blast_word_sz} -db ${db_file}/salmonella_genome_db -out ${outputFilename} -culling_limit ${params.blastCullingLimit} -max_hsps ${params.blastHSQSLimit} -outfmt "6 ${blastOuputColumns.join(' ')}"
    sed -i '1i ${blastOuputColumns.join('\t')}' ${outputFilename}
    """
}

process filterHIghCovarageAlignments {
    input:
    file alignmentsFile

    output:
    file highCovarageAlignmentsFile

    script:
    highCovarageAlignmentsFile = "${alignmentsFile.baseName}_highCovarage.tsv"
    """
    awk -F'\t' '{if (\$41 > ${params.minAlignmentCoverage}) print}' ${alignmentsFile} > ${highCovarageAlignmentsFile}
    """
}

workflow blast_wf {

    take:
    queries_ch 
    subjectGenomes_ch

    main:
    wordSizes_ch = channel.from(params.wordSizes_list)
    subjectDB_ch = subjectGenomes_ch | indexSubjectSequences

    blastInputs_ch = queries_ch
        .splitFasta(by: params.queriesChunckSize, file:true)
        .combine(subjectDB_ch)
        .combine(wordSizes_ch)
        
    alignments_ch = alignLocally(blastInputs_ch)
    groupedAlignments_ch = alignments_ch.collectFile(
                            item -> [item[0], item[1]],
                            keepHeader: true,
                        )

    highCovarageAlignments_ch = filterHIghCovarageAlignments(groupedAlignments_ch)

    emit:
    aligmentsResults_ch = highCovarageAlignments_ch
}

