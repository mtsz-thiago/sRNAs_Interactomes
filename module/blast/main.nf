#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.queriesChunckSize = 10
params.wordSizes_list = [4,7,11,15]

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
    tuple val(prefix), path("alignments_results.tsv")

    script:
    prefix = getScenarioFromFileName(query_file) + "w${blast_word_sz}"
    """
    blastn -query ${query_file} -word_size ${blast_word_sz} -db ${db_file}/salmonella_genome_db -out alignments_results.tsv -outfmt "6 qseqid qgi qacc qaccver qlen sseqid sallseqid sgi sallgi sacc saccver sallacc slen qstart qend sstart send qseq sseq evalue bitscore score length pident nident mismatch positive gapopen gaps ppos frames qframe sframe btop staxids sscinames scomnames sblastnames sskingdoms stitle salltitles sstrand qcovs qcovhsp"
    sed -i '1i qseqid\tqgi\tqacc\tqaccver\tqlen\tsseqid\tsallseqid\tsgi\tsallgi\tsacc\tsaccver\tsallacc\tslen\tqstart\tqend\tsstart\tsend\tqseq\tsseq\tevalue\tbitscore\tscore\tlength\tpident\tnident\tmismatch\tpositive\tgapopen\tgaps\tppos\tframes\tqframe\tsframe\tbtop\tstaxids\tsscinames\tscomnames\tsblastnames\tsskingdoms\tstitle\tsalltitles\tsstrand\tqcovs\tqcovhsp' alignments_results.tsv
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

    emit:
    aligmentsResults_ch = groupedAlignments_ch
}

