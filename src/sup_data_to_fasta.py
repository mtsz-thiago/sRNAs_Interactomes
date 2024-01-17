import argparse
import pandas as pd
import os
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Blast import NCBIWWW, NCBIXML


def load_xlsx(input_file):
    sup_dict = pd.read_excel(input_file, sheet_name=None)
    del sup_dict['Legend']
    return sup_dict

def map_record_to_SeqRecord(r, col_name="RNA1"):
    strand_code = 1 if r[f"{col_name} Strand"] == "+" else -1
    id = r[f"{col_name} name"]
    name = r[f"{col_name} name"]
    seq = Seq(r[f"{col_name} seq"])
    from_pos = r[f"{col_name} from"]
    to_pos = r[f"{col_name} to"]
    type = r[f"{col_name} type"]
    description = id + ' ' + ' '
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
        SeqIO.write(df, f"{output}/{k}.fa", "fasta")

def main(input_file, output_dir):    
    if not os.path.exists(output_dir): 
        os.makedirs(output_dir)
            
    sup_data_dict = load_xlsx(input_file)
    seg_records_dict = map_SeqRecords_to_fasta(sup_data_dict)
    write_fasta(seg_records_dict, output_dir)
   

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert supplementary data to FASTA format.")
    parser.add_argument("-i", "--input", help="Input .xlsx file with Liu, et al. supplementary data")
    parser.add_argument("-o", "--output_dir", default=os.getcwd(), help="output_dir directory")
    args = parser.parse_args()
    main(args.input, args.output_dir)



