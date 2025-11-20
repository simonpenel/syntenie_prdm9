import pandas as pd
import argparse
from collections import Counter
parser = argparse.ArgumentParser()

parser.add_argument('-i', '--id2spec', type=str, required=True, help='ID2SPEC input file')
parser.add_argument('-f', '--fam2id', type=str, required=True, help='FAM2ID input file')
parser.add_argument('-o', '--fam2idspec', type=str, required=True, help='FAM2IDSPEC output file')

args = parser.parse_args()

df_id2spec = pd.read_csv(args.id2spec, sep=';', header=0,names=["SeqID","Species"])
#print(df_id2spec)
df_fam2id = pd.read_csv(args.fam2id, sep=';', header=0,names=["Family","SeqID"])
#print(df_fam2id)
df_fam2idspec = df_fam2id.join(df_id2spec.set_index('SeqID'),on='SeqID', how="inner",lsuffix='_caller', rsuffix='_other')
#df_fam2idspec.to_csv(args.fam2idspec,sep=';' , index=False , na_rep="NA")
fam_nbseq = Counter(df_fam2id["Family"])
df_fam_nbseq = pd.DataFrame.from_dict(fam_nbseq,orient='index').reset_index()
df_fam_nbseq = df_fam_nbseq.set_axis(['Family', 'Nb sequences'], axis=1)
#df_fam_nbseq.to_csv("lol.csv",sep=';' , index=False , na_rep="NA")
df_fam2spec = df_fam2idspec.drop("SeqID", axis=1).drop_duplicates()
fam_nbspec = Counter(df_fam2spec["Family"])
df_fam_nbspec = pd.DataFrame.from_dict(fam_nbspec,orient='index').reset_index()
df_fam_nbspec = df_fam_nbspec.set_axis(['Family', 'Nb species'], axis=1)
#df_fam_nbspec.to_csv("lal.csv",sep=';' , index=False , na_rep="NA")
df_fam_nbseqspec = df_fam_nbseq.join(df_fam_nbspec.set_index('Family'),on='Family', how="inner",lsuffix='_caller', rsuffix='_other')
df_fam_nbseqspec.to_csv(args.fam2idspec,sep=';' , index=False , na_rep="NA")