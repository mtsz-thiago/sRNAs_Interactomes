import argparse
import pandas as pd
import os
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Blast import NCBIWWW, NCBIXML

def load_inputs(input_files):    
    data_dfs = { Path(f).stem:pd.read_csv(f) for f in input_files}
    return data_dfs

def map_record_to_SeqRecord(r, col_name="RNA1"):
    strand_code = 1 if r[f"{col_name} Strand"] == "+" else -1
    id = r[f"{col_name} name"]
    name = r[f"{col_name} name"]
    seq = Seq(r[f"{col_name} seq"])
    from_pos = r[f"{col_name} from"]
    to_pos = r[f"{col_name} to"]
    type = r[f"{col_name} type"]
    description = '>' + id + ' ' + ' '
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
    for k, df in seq_records.items():
        bio_dfs[k+"_RNA1"] = df.apply(axis=1, func=lambda x: map_record_to_SeqRecord(x, "RNA1"))
        bio_dfs[k+"_RNA2"] = df.apply(axis=1, func=lambda x: map_record_to_SeqRecord(x, "RNA2"))
    return bio_dfs

def write_fasta(seq_records, output):
    for k, df in seq_records.items():
        SeqIO.write(df, f"{output}/{k}.fasta", "fasta")

def main(input_files, output):
    
    if not os.path.exists(output): 
        os.makedirs(output)
            
    input_files_dict = load_inputs(input_files)
    seg_records_dict = map_SeqRecords_to_fasta(input_files_dict)
    write_fasta(seg_records_dict, output)
   

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert supplementary data to FASTA format.")
    parser.add_argument("-i", "--inputs", nargs="+", help="Input CSV files")
    parser.add_argument("-o", "--output", help="Output directory")
    args = parser.parse_args()
    main(args.inputs, args.output)



