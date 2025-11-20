import pandas as pd
import argparse
from collections import Counter
parser = argparse.ArgumentParser()

parser.add_argument('-n', '--nb_species', type=int, required=True, help='nb of species')
parser.add_argument('-s', '--fam2nbseqspec', type=str, required=True, help='FAM2NBSEQNBSPEC inputfile')
parser.add_argument('-f', '--fam2id', type=str, required=True, help='FAM2ID input file')
parser.add_argument('-o', '--monofam2id', type=str, required=True, help=' monogenic FAM2IDSPEC output file')

args = parser.parse_args()

df_fam_nbseqspec = pd.read_csv(args.fam2nbseqspec, sep=';', header=0)
print(df_fam_nbseqspec)
df_fam2id = pd.read_csv(args.fam2id, sep=';', header=0,names=["Family","SeqID"])
print(df_fam2id)
nb_species = args.nb_species
print("Nb of species = "+ str(nb_species))
# pas plus de 3xplus de seuqnec que d'especes
# au moins 1/3 des especes
df_monogenic_fam_nbseqspec  = df_fam_nbseqspec[df_fam_nbseqspec ['Nb sequences'] <= df_fam_nbseqspec['Nb species'] * 3 ]
df_broad_monogenic_fam_nbseqspec  = df_monogenic_fam_nbseqspec[df_monogenic_fam_nbseqspec['Nb sequences'] >= int( nb_species / 3) ]

print(df_broad_monogenic_fam_nbseqspec)
df_fam2id_monogenic = df_fam2id.join(df_broad_monogenic_fam_nbseqspec.set_index('Family'),on='Family', how="inner",lsuffix='_caller', rsuffix='_other')
print(df_fam2id_monogenic)
df_fam2id_monogenic.to_csv(args.monofam2id,sep=';' , index=False)

# #df_fam2idspec.to_csv(args.fam2idspec,sep=';' , index=False , na_rep="NA")
# fam_nbseq = Counter(df_fam2id["Family"])
# df_fam_nbseq = pd.DataFrame.from_dict(fam_nbseq,orient='index').reset_index()
# df_fam_nbseq = df_fam_nbseq.set_axis(['Family', 'Nb sequences'], axis=1)
# #df_fam_nbseq.to_csv("lol.csv",sep=';' , index=False , na_rep="NA")
# df_fam2spec = df_fam2idspec.drop("SeqID", axis=1).drop_duplicates()
# fam_nbspec = Counter(df_fam2spec["Family"])
# df_fam_nbspec = pd.DataFrame.from_dict(fam_nbspec,orient='index').reset_index()
# df_fam_nbspec = df_fam_nbspec.set_axis(['Family', 'Nb species'], axis=1)
# #df_fam_nbspec.to_csv("lal.csv",sep=';' , index=False , na_rep="NA")
# df_fam_nbseqspec = df_fam_nbseq.join(df_fam_nbspec.set_index('Family'),on='Family', how="inner",lsuffix='_caller', rsuffix='_other')
# df_fam_nbseqspec.to_csv(args.fam2idspec,sep=';' , index=False , na_rep="NA")