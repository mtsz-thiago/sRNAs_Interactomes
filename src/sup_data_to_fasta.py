#!/usr/bin/env python

import argparse
import pandas as pd
import os
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
<<<<<<< HEAD
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Blast import NCBIWWW, NCBIXML
import uuid

def deduplicate_queries_by_seq(df):
    duplicated_seqs = df.seq.duplicated()
    return df[~duplicated_seqs]

=======
import uuid

>>>>>>> feature/starting_graphs_modeling
def flatten_df(df):
    
    dfRNA1 = df[[c for c in df.columns if 'RNA2' not in c]]
    rna1columnsMap = {c:c.strip('RNA1').strip() for c in dfRNA1.columns if 'RNA1' in c}
    dfRNA1 = dfRNA1.rename(columns=rna1columnsMap)
    dfRNA1['origin'] = 'RNA1'
<<<<<<< HEAD
    dfRNA2 = df[[c for c in df.columns if 'RNA1' not in c]]
    rna2columnsMap = {c:c.strip('RNA2').strip() for c in dfRNA2.columns if 'RNA2' in c}
    dfRNA2 = dfRNA2.rename(columns=rna2columnsMap)
    dfRNA1['origin'] = 'RNA1'
    
    flattened_df = pd.concat([dfRNA1, dfRNA2], ignore_index=True)
    
=======
    dfRNA1['chimera_idx'] = dfRNA1.index
    
    dfRNA2 = df[[c for c in df.columns if 'RNA1' not in c]]
    rna2columnsMap = {c:c.strip('RNA2').strip() for c in dfRNA2.columns if 'RNA2' in c}
    dfRNA2 = dfRNA2.rename(columns=rna2columnsMap)
    dfRNA2['origin'] = 'RNA2'
    dfRNA2['chimera_idx'] = dfRNA2.index
    
    flattened_df = pd.concat([dfRNA1, dfRNA2], ignore_index=True)
    
>>>>>>> feature/starting_graphs_modeling
    return flattened_df

def extract_queries_and_expected(df):
    flattened_df = flatten_df(df)
<<<<<<< HEAD
    deduplicated_df = deduplicate_queries_by_seq(flattened_df)
    with_query_id = add_query_id_column(deduplicated_df)
=======
    with_query_id = add_query_id_column(flattened_df)
>>>>>>> feature/starting_graphs_modeling
    bio_df = map_to_fasta(with_query_id)
    return bio_df, with_query_id

def map_raw_data_to_queries_and_expected(df_dict):
    queries_dict = {}
    expected_dict = {}
    for k, df in df_dict.items():
<<<<<<< HEAD
        bio_df, expected_df = extract_queries_and_expected(df)
        queries_dict[k] = bio_df
        expected_dict[k] = expected_df
    return queries_dict, expected_dict

def add_query_id_column(df):
    df['query_id'] = df.apply(lambda x: str(uuid.uuid1()), axis=1)
=======
        bio_df, expected = extract_queries_and_expected(df)
        queries_dict[k] = bio_df
        expected_dict[k] = expected
    return queries_dict, expected_dict      

def add_query_id_column(df):
    df['query_id'] = df.groupby('seq').name.transform(lambda x: uuid.uuid4().hex)
>>>>>>> feature/starting_graphs_modeling
    return df

def load_xlsx(input_file):
    sup_dict = pd.read_excel(input_file, sheet_name=None)
    sup_dict = {k.replace('+', '_'): v for k, v in sup_dict.items() if k != 'Legend'}
    return sup_dict

def map_to_fasta(df):
<<<<<<< HEAD
    return df.apply(axis=1, func=map_record_to_SeqRecord)

def map_record_to_SeqRecord(r):
    strand_code = 1 if r[f"Strand"] == "+" else -1
    id = r['query_id']
    name = id
    seq = Seq(r[f"seq"])
    from_pos = r[f"from"]
    to_pos = r[f"to"]
    type = r[f"type"]
    description = f"{id} origin: {r.origin} expected {r.name}"
=======
    bio_df = df.groupby('seq').apply(lambda x: map_record_to_SeqRecord(x), include_groups=False)
    return bio_df

def map_record_to_SeqRecord(r):
    id = f"{r.query_id.iloc[0]}"
    # print(id)
    name = ''
    seq = Seq(r.name)
    description = f"{id} from {r.chimera_idx.values}"
>>>>>>> feature/starting_graphs_modeling
    return SeqRecord(
        id=id,
        description=description,
        name=name, 
        seq=seq)

def write_csvs(seq_records, expected, output):
    
    if not os.path.exists(f"{output}/queries"): 
        os.makedirs(f"{output}/queries")
    if not os.path.exists(f"{output}/expected"): 
        os.makedirs(f"{output}/expected")
    
    for k, df in seq_records.items():
        SeqIO.write(df, f"{output}/queries/{k}.fa", "fasta")
        
    for k, df in expected.items():
        df.to_csv(f"{output}/expected/{k}.csv", index=False)

def main(input_file, output_dir):    
    if not os.path.exists(output_dir): 
        os.makedirs(output_dir)
    
    sup_data_dict = load_xlsx(input_file)
    
    queries_dict, expected_dict = map_raw_data_to_queries_and_expected(sup_data_dict)

    write_csvs(queries_dict, expected_dict, output_dir)
<<<<<<< HEAD
   
=======
>>>>>>> feature/starting_graphs_modeling

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert supplementary data to FASTA format.")
    parser.add_argument("-i", "--input", help="Input .xlsx file with Liu, et al. supplementary data")
    parser.add_argument("-o", "--output_dir", default=os.getcwd(), help="output_dir directory")
    args = parser.parse_args()
    main(args.input, args.output_dir)



