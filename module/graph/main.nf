#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.kmer_sz = 4
params.neo4jURI = "bolt://neo4j:7687"
params.neo4jUser = "neo4j"
params.neo4jPassword = "Password"
params.neo4jDB = "neo4j"

chimerasNodesColumns = ['name','Strand','from','to','type','seq','query_id']
chimerasEdgesColumns = ['ligation from','ligation to','Number of interactions','Odds Ratio','fishersPValue']

alignmentsNodesColumns = ['sseqid', 'sgi', 'sacc', 'saccver', 'slen', 'sseq', 'staxids', 'sscinames', 'scomnames', 'sblastnames', 'sskingdoms', 'stitle', 'sstrand']
alignmentsEdgesColumns = ['qseqid', 'qstart', 'qend', 'sstart', 'send', 'qseq', 'evalue', 'bitscore', 'score', 'length', 'pident', 'nident', 'mismatch', 'positive', 'gapopen', 'gaps', 'ppos', 'qframe', 'sframe', 'btop', 'qcovs', 'qcovhsp']

process loadChimerasToDB {

    input:
    path chimerasFile

    output:
    val graphName

    script:
    graphName = chimerasFile.baseName.split('-')[0]
    """
    #!/usr/bin/env python3

    import pandas as pd
    from neo4j import GraphDatabase

    NEO4J_AUTH = ('${params.neo4jUser}','${params.neo4jPassword}')
    driver = GraphDatabase.driver('${params.neo4jURI}', auth=NEO4J_AUTH, database='${params.neo4jDB}')

    graphName = '${chimerasFile.baseName}'

    chimeras_df = pd.read_csv('${chimerasFile}')
    chimeras_df.rename(columns={"Fisher\'s exact test p-value": 'fishersPValue'}, inplace=True)
    node_columns = '${chimerasNodesColumns.join(', ')}'.split(', ')
    edge_columns = '${chimerasEdgesColumns.join(', ')}'.split(', ')
    nodes_data = chimeras_df[node_columns]
    edges_data = chimeras_df[edge_columns]

    nodes = nodes_data.assign(
        nodeId=nodes_data.index,
        labels=lambda x: [["SEQUENCE", "CHIMERA", "${graphName}"]] * len(nodes_data)
    )

    def find_pair(i, r, chimeras_df):
        return chimeras_df.query(f"({r['chimera_idx']} == chimera_idx) and ({i} != index)").index[0]

    targetNode = [find_pair(i, r, chimeras_df) for i,r in chimeras_df.iterrows() ]
    edges = edges_data.assign(
        sourceNode=edges_data.index,
        targetNode=targetNode
    )
    
    with driver.session() as session:
        # Load nodes
        for index, row in nodes.iterrows():
            session.run(
                "CREATE (n:SEQUENCE:CHIMERA:${graphName} {nodeId: \$nodeId, name: \$name, Strand: \$Strand, from: \$from, to: \$to, type: \$type, seq: \$seq, query_id: \$query_id})",
                **row
            )

        # Load edges
        for index, row in edges.iterrows():
            session.run(
                "MATCH (a:SEQUENCE:CHIMERA {nodeId: \$sourceNode}), (b:SEQUENCE:CHIMERA {nodeId: \$targetNode}) "
                "CREATE (a)-[r:LIGATES {Number_of_interactions: \$Number_of_interactions, Odds_Ratio: \$Odds_Ratio, fishersPValue: \$fishersPValue}]->(b)",
                sourceNode=row['sourceNode'], targetNode=row['targetNode'], Number_of_interactions=row['Number of interactions'], Odds_Ratio=row['Odds Ratio'], fishersPValue=row['fishersPValue']
            )
    """

}

process loadAlignmentsToDB {

    input:
    tuple val(graphName), path(alignmentsFile)

    output:
    val(graphName)

    script:
    """
    #!/usr/bin/env python3

    import pandas as pd
    from neo4j import GraphDatabase

    NEO4J_AUTH = ('${params.neo4jUser}','${params.neo4jPassword}')
    driver = GraphDatabase.driver('${params.neo4jURI}', auth=NEO4J_AUTH, database='${params.neo4jDB}')

    graphName = '${alignmentsFile.baseName}'

    alignments_df = pd.read_csv('${alignmentsFile}', sep='\t')

    node_columns = '${alignmentsNodesColumns.join(', ')}'.split(', ')
    edge_columns = '${alignmentsEdgesColumns.join(', ')}'.split(', ')

    nodes_data = alignments_df[node_columns]
    edges_data = alignments_df[edge_columns]

    nodes = nodes_data.assign(
        nodeId=nodes_data.index,
        labels=lambda x: [["SEQUENCE", "ALIGNMENT", "${graphName}"]] * len(nodes_data)
    )

    edges = edges_data.assign(
        sourceNode=edges_data.index,
        targetNode=edges_data['qseqid']
    )

    with driver.session() as session:
        # Load nodes
        for index, row in nodes.iterrows():
            session.run(
                "CREATE (n:SEQUENCE:ALIGNMENT:${graphName} {nodeId: \$nodeId, sseqid: \$sseqid, sgi: \$sgi, sacc: \$sacc, saccver: \$saccver, slen: \$slen, sseq: \$sseq, staxids: \$staxids, sscinames: \$sscinames, scomnames: \$scomnames, sblastnames: \$sblastnames, sskingdoms: \$sskingdoms, stitle: \$stitle, sstrand: \$sstrand})",
                **row
            )

        # Load edges
        for index, row in edges.iterrows():
            session.run(
                "MATCH (a:SEQUENCE:ALIGNMENT:${graphName} {nodeId: \$sourceNode}), (b:SEQUENCE:CHIMERA:${graphName} {query_id: \$qseqid}) "
                "CREATE (a)-[r:ALIGNS {qseqid: \$qseqid, qstart: \$qstart, qend: \$qend, sstart: \$sstart, send: \$send, qseq: \$qseq, evalue: \$evalue, bitscore: \$bitscore, score: \$score, length: \$length, pident: \$pident, nident: \$nident, mismatch: \$mismatch, positive: \$positive, gapopen: \$gapopen, gaps: \$gaps, ppos: \$ppos, qframe: \$qframe, sframe: \$sframe, btop: \$btop, qcovs: \$qcovs, qcovhsp: \$qcovhsp}]->(b)",
                **row
            )

    """
}


workflow interactomeModeling_wf {

    take:
    chimeras_ch
    alignments_ch

    main:
    loadResults_ch = loadChimerasToDB(chimeras_ch)
    
    keyALignments_ch = alignments_ch.map(it-> [it.baseName.split('-')[0], it]).join(loadResults_ch)
    loadAlignmentsToDB(keyALignments_ch)
}

