#!/usr/bin/env python

import argparse
import pandas as pd
import os
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Blast import NCBIWWW, NCBIXML
import uuid

def deduplicate_queries_by_seq(df):
    duplicated_seqs = df.seq.duplicated()
    return df[~duplicated_seqs]

def from_pair_to_single_seq_record_df(df):
    
    RNA1_related_columns = [k for k in df.columns if "RNA1" in k]
    RNA2_related_columns = [k for k in df.columns if "RNA2" in k]
    not_RNA1_or_RNA2_related = [k for k in df.columns if k not in RNA1_related_columns and k not in RNA2_related_columns]
    common_columns_suffixes = [k.split(" ")[1] for k in RNA1_related_columns]

    rna1_df = df[not_RNA1_or_RNA2_related + RNA1_related_columns].copy()
    rna2_df = df[not_RNA1_or_RNA2_related + RNA2_related_columns].copy()

    rna1_suffixes_dict = {k:v for k,v in zip(RNA1_related_columns, common_columns_suffixes)}
    rna1_df.rename(columns=rna1_suffixes_dict, inplace=True)
    rna1_df["origin"] = "RNA1"

    rna2_suffixes_dict = {k:v for k,v in zip(RNA2_related_columns, common_columns_suffixes)}
    rna2_df.rename(columns=rna2_suffixes_dict, inplace=True)
    rna2_df["origin"] = "RNA2"
    
    queries_df = pd.concat([rna1_df, rna2_df], ignore_index=True)
    large_seq_queries_df = deduplicate_queries_by_seq(queries_df)
    
    return large_seq_queries_df

def transform_raw_dfs_to_queries(raw_dfs_dict):
    return {k: from_pair_to_single_seq_record_df(v) for k,v in raw_dfs_dict.items()}

def load_xlsx(input_file):
    sup_dict = pd.read_excel(input_file, sheet_name=None)
    sup_dict = {k.replace('+', '_'): v for k, v in sup_dict.items() if k != 'Legend'}
    return sup_dict

def map_record_to_SeqRecord(r):
    strand_code = 1 if r[f"Strand"] == "+" else -1
    id = str(uuid.uuid4())
    name = id + "_" + r[f"origin"]
    seq = Seq(r[f"seq"])
    from_pos = r[f"from"]
    to_pos = r[f"to"]
    type = r[f"type"]
    description = f"{id} origin: {r.origin} expected {r.name}"
    return SeqRecord(
        id=id,
        description=description,
        features=[
                SeqFeature(FeatureLocation(from_pos, to_pos, strand=strand_code),type=type)
            ],
        name=name, 
        seq=seq)

def map_SeqRecords_to_fasta(seq_records):
    bio_dfs = {}
    expected_dfs = {}
    for k, df in seq_records.items():
        
        df1 = df[df.origin == 'RNA1']
        df2 = df[df.origin == 'RNA2']
        
        bio_dfs[k+"_RNA1"] = df1.apply(axis=1, func=lambda x: map_record_to_SeqRecord(x))
        bio_dfs[k+"_RNA2"] = df2.apply(axis=1, func=lambda x: map_record_to_SeqRecord(x))
        
        
        expected_dfs[k+"_RNA1"] = pd.concat([df1.name, bio_dfs[k+"_RNA1"].map(lambda i: i.id)], axis=1)
        expected_dfs[k+"_RNA1"].rename(columns={"name": "expected", 0: "query_id"}, inplace=True)
        expected_dfs[k+"_RNA2"] = pd.concat([df1.name, bio_dfs[k+"_RNA2"].map(lambda i: i.id)], axis=1)
        expected_dfs[k+"_RNA2"].rename(columns={"name": "expected", 0: "query_id"}, inplace=True)
        
    return bio_dfs, expected_dfs

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
    
    queries_dict = transform_raw_dfs_to_queries(sup_data_dict)

    seg_records_dict, expected_dfs = map_SeqRecords_to_fasta(queries_dict)
    write_csvs(seg_records_dict, expected_dfs, output_dir)
   

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert supplementary data to FASTA format.")
    parser.add_argument("-i", "--input", help="Input .xlsx file with Liu, et al. supplementary data")
    parser.add_argument("-o", "--output_dir", default=os.getcwd(), help="output_dir directory")
    args = parser.parse_args()
    main(args.input, args.output_dir)



