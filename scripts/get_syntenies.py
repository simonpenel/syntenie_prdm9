import argparse
import time
import os
import sys
import pandas as pd
from Bio import SeqIO
 
parser = argparse.ArgumentParser()


parser.add_argument('-i', '--input', type=str, required=True, help='whole_summary.csv from protein analysis')
parser.add_argument('-p', '--positions', type=str,required=True, help='path to gff file parsed results on gene position')
parser.add_argument('-o', '--output', type=str, required=True, help='output file path')
parser.add_argument('-c', '--cluster', type=str, required=True, help='input cluster path')
parser.add_argument('-s', '--stats', type=str, required=True, help='input cluster stats')

args = parser.parse_args()

nb_voisins = 50

def fun_gene_2_prot(dico, gene:str):
    if gene == "NO DATA":
        return("NO DATA")
    if gene in dico :
        proteins = dico[gene]
    else :
        print(gen)
        sys.exit("Gene non trouve")
        #proteins = ["NO PROT"]
    return(proteins)

def fun_prot_2_fam(dico, proteins):
    if proteins == "NO DATA":
        return("NO DATA")
    families = []
    for protein in proteins:
        if protein in dico :
            family  = dico[protein]
        else :
            family = ["NO FAM", 0, 0]
        families.append(family)
    return(families)

def fun_unique_families(families):
    if families == "NO DATA":
        return("NO DATA")
    unique_families = []
    for family in families:
        if not family  in unique_families :
            unique_families.append(family)
    return(unique_families)    


def fun_largest_families(families):
    if families == "NO DATA":
        return("NO DATA")
    largest_family = families[0]
    for family in families:
        if family[2] > largest_family[2]:
            largest_family  = family
        if (family[2] == largest_family[2])  and (family[1] < largest_family[1]):
            largest_family  = family           
    return(largest_family[0])    


def processCluster(fcluster:str):
    #df_cluster = pd.read_csv(fcluster, sep=';',header=None,names=["Family","SeqID"])
    df_cluster = pd.read_csv(fcluster, sep=';',header=0)
    dico = df_cluster.set_index('SeqID').T.to_dict('list')
    return(dico)

def processStats(fstats:str):
    #df_cluster = pd.read_csv(fcluster, sep=';',header=None,names=["Family","SeqID"])
    df_stats = pd.read_csv(fstats, sep=';',header=0)
    dico = df_stats.set_index('Family').T.to_dict('list')
    return(dico)

def processPositions(gff:str,
    dico_prot_2_contig,
    dico_prot_2_gene,
    dico_contig_2_synteny,
    dico_contig_2_indexes,
    dico_gene_2_proteins
    ):
    # dico_prot_2_contig = {}  # prot name -> contig ( on ne gadre que le 1er contig)
    # dico_prot_2_gene = {}    # prot_name -> gene_id
    # dico_contig_2_synteny = {} # contig -> arrayf of genes
    # dico_contig_2_indexes = {} # contig -> dico gene-> index in array
    with open(gff, 'r') as reader:
        lines = reader.readlines()[1:]
        for line in lines:
            info  = line.rstrip('\n').split(";")
            assembly = info[0]
            prot_name = info[1]
            contig = info[2]
            gene_name = info[3]
            gene_id = info[4]
            dico_prot_2_gene[prot_name] = gene_id
            if not gene_id in dico_gene_2_proteins:
                dico_gene_2_proteins[gene_id] = []
                dico_gene_2_proteins[gene_id].append(prot_name)
            else :
                dico_gene_2_proteins[gene_id].append(prot_name)
            if prot_name in dico_prot_2_contig :
                print("Skiping " + prot_name + " in "+contig+" because  it is already in "+dico_prot_2_contig[prot_name])
            else :
                dico_prot_2_contig[prot_name] = contig
            if not contig in dico_contig_2_synteny:
                dico_contig_2_synteny[contig] = []
                dico_contig_2_synteny[contig].append(gene_id)
            else :
                if not gene_id in dico_contig_2_synteny[contig]:
                    dico_contig_2_synteny[contig].append(gene_id)
    for contig in dico_contig_2_synteny.keys() :
        index_gene = 0
        dico_index = {}
        for gene in dico_contig_2_synteny[contig]:
            dico_index[gene] = index_gene
            index_gene +=1
        dico_contig_2_indexes[contig] = dico_index


dico_prot_2_contig = {}  # prot name -> contig ( on ne gadre que le 1er contig)
dico_prot_2_gene = {}    # prot_name -> gene_id
dico_contig_2_synteny = {} # contig -> arrayf of genes
dico_contig_2_indexes = {} # contig -> dico gene-> index in array
dico_gene_2_proteins = {} #  gene-> array of proteins

processPositions(args.positions,dico_prot_2_contig,dico_prot_2_gene,dico_contig_2_synteny,dico_contig_2_indexes,dico_gene_2_proteins)                            

dico_cluster = processCluster(args.cluster)
dico_stats = processStats(args.stats)


logfile=open(args.output+".log", 'w')
outfile=open(args.output, 'w')
outfile.write("SeqID")
voisin = -nb_voisins + 1
while(voisin < 0):
    outfile.write(";"+str(voisin))
    voisin +=1
voisin = 1
while(voisin < nb_voisins):
    outfile.write(";"+str(voisin))
    voisin +=1
outfile.write("\n")

outfilecheck=open(args.output+".check", 'w')
outfilecheck.write("SeqID")
voisin = -nb_voisins + 1
while(voisin < 0):
    outfilecheck.write(";"+str(voisin))
    voisin +=1
voisin = 1
while(voisin < nb_voisins):
    outfilecheck.write(";"+str(voisin))
    voisin +=1
outfilecheck.write("\n")



logfile.write("PROCESSING PROTEINS\n")
logfile.write("===================\n")
with open(args.input, 'r') as reader:
    lines = reader.readlines()[1:]
    for line in lines:
        splitline = line.split(';')
        prot_name = splitline[1]
        logfile.write("Processing "+ prot_name+"\n")
        print("Processing "+ prot_name)
        contig = dico_prot_2_contig[prot_name]
        print("Contig = "+contig)
        gene_id = dico_prot_2_gene[prot_name]
        print("GeneID = "+gene_id)
        synteny = dico_contig_2_synteny[contig]
        indexes = dico_contig_2_indexes[contig]

        index = indexes[gene_id]
        print("Index " + gene_id + " = "+str(index))
        if synteny[index] != gene_id:
            sys.exit("error in index")
        left_genes = []
        for voisin in  range(index - nb_voisins, index -1):
            print(str(voisin))
            if voisin >= 0 :
                print("Gene id "+ synteny[voisin] )
                left_genes.append(synteny[voisin])
            else :
                left_genes.append("NO DATA")
        print(left_genes)


        proteins = list(map(lambda x: fun_gene_2_prot(dico_gene_2_proteins,x), left_genes))
        print(proteins)
        families = list(map(lambda x: fun_prot_2_fam(dico_cluster,x), proteins))
        print(families)
        unique_families = list(map(fun_unique_families, families))
        print(unique_families)
        largest_families = list(map(fun_largest_families, unique_families))
        print(largest_families)

        # get contig
        




#df.to_csv(args.output, sep=';')



