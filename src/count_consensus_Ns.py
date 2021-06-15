
#!/usr/bin/env python3

## Timothy Thurman


## Usage ~/path_to_script/count_consensus_Ns.py -i SOMETHING.fa -o OUTPUT.csv

# import re # regular expressions
import argparse # to parse arguments
import pyfastx
import pandas as pd

# Do argument parsing
parser = argparse.ArgumentParser(description="""Count base composition, particularly Ns, in consensus seq""")
parser.add_argument("-i", "--input", dest="inFile", help="Consensus seqeunces, in multi-sample .fasta")
parser.add_argument("-o", "--output", dest="outFile", help="Text file for results output")
args = parser.parse_args()



res_df = pd.DataFrame({"sample" : [None], "A" : [None], "C" : [None], "G" : [None], "N" : [None], "T" : [None], "total" : [None]})
for seq in pyfastx.Fasta(args.inFile):
        name_dict = {"sample" :seq.name}
        length_dict = {"total" : len(seq)}
        comp = seq.composition
        comp.update(name_dict)
        comp.update(length_dict)
        one_res_df = pd.DataFrame.from_records([comp])
        res_df = res_df.append(one_res_df)
        
res_df = res_df.iloc[1: , :]

res_df['perc_N'] = res_df['N']/res_df['total']
res_df['n_less5k'] = res_df['N'] < 5000


res_df.to_csv(args.outFile, index=False) 