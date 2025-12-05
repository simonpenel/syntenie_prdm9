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
            family = ["NO FAMILY", 0, 0]
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

def genes_2_families(genes,dico_gene_2_proteins,dico_cluster):
    proteins = list(map(lambda x: fun_gene_2_prot(dico_gene_2_proteins,x), genes))
    families = list(map(lambda x: fun_prot_2_fam(dico_cluster,x), proteins))
    unique_families = list(map(fun_unique_families, families))
    largest_families = list(map(fun_largest_families, unique_families))
    return(largest_families)
    

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
    dico_gene_2_proteins,
    dico_geneid_2genename
    ):
    with open(gff, 'r') as reader:
        lines = reader.readlines()[1:]
        for line in lines:
            info  = line.rstrip('\n').split(";")
            assembly = info[0]
            prot_name = info[1]
            contig = info[2]
            gene_name = info[3]
            gene_id = info[4]
            if not gene_id in dico_geneid_2genename:
                dico_geneid_2genename[gene_id] = gene_name
            else :
                if dico_geneid_2genename[gene_id] != gene_name :
                    sys.exit("Error in gene names")
            dico_prot_2_gene[prot_name] = gene_id
            if not gene_id in dico_gene_2_proteins:
                dico_gene_2_proteins[gene_id] = []
                if not prot_name in dico_gene_2_proteins[gene_id] :
                    dico_gene_2_proteins[gene_id].append(prot_name)
            else :
                if not prot_name in dico_gene_2_proteins[gene_id] :
                    dico_gene_2_proteins[gene_id].append(prot_name)
            if prot_name in dico_prot_2_contig :
                print("Skiping " + prot_name + " in "+contig+" because  it is already in "+dico_prot_2_contig[prot_name],file=sys.stderr)
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
            #print("DEBUG INDEX "+contig+" gene = "+gene+" name = "+dico_geneid_2genename[gene] +" index = "+str(index_gene))
            dico_index[gene] = index_gene
            index_gene +=1
        dico_contig_2_indexes[contig] = dico_index


dico_prot_2_contig = {}  # prot name -> contig ( on ne gadre que le 1er contig)
dico_prot_2_gene = {}    # prot_name -> gene_id
dico_contig_2_synteny = {} # contig -> arrayf of genes
dico_contig_2_indexes = {} # contig -> dico gene-> index in array
dico_gene_2_proteins = {} #  gene-> array of proteins
dico_geneid_2genename = {} #  gene id > gene name

print("Processing "+args.positions)
processPositions(
    args.positions,
    dico_prot_2_contig,
    dico_prot_2_gene,
    dico_contig_2_synteny,
    dico_contig_2_indexes,
    dico_gene_2_proteins,
    dico_geneid_2genename)        
print("Done.")                    

print("Processing "+args.cluster)
dico_cluster = processCluster(args.cluster)
print("Done.") 
print("Processing "+args.stats)
dico_stats = processStats(args.stats)
print("Done.") 

logfile=open(args.output+".log", 'w')
outfile=open(args.output, 'w')
outfile.write("SeqID")
voisin = -nb_voisins 
while(voisin < 0):
    outfile.write(";"+str(voisin))
    voisin +=1
voisin = 1
while(voisin <= nb_voisins):
    outfile.write(";"+str(voisin))
    voisin +=1
outfile.write("\n")


with open(args.input, 'r') as reader:
    lines = reader.readlines()[1:]
    for line in lines:
        splitline = line.split(';')
        prot_name = splitline[1]
        print("Processing "+ prot_name)
        logfile.write("Processing "+ prot_name+" : ")
        contig = dico_prot_2_contig[prot_name]
        logfile.write("Contig = "+contig+" ")
        gene_id = dico_prot_2_gene[prot_name]
        logfile.write("GeneID = "+gene_id+" ")
        gene_name = dico_geneid_2genename[gene_id]
        logfile.write("Gene name = "+gene_name+" ")
        synteny = dico_contig_2_synteny[contig]
        indexes = dico_contig_2_indexes[contig]

        index = indexes[gene_id]
        logfile.write("Index " + gene_id + " = "+str(index)+"\n")
        if synteny[index] != gene_id:
            sys.exit("error in index")
        left_genes = []
        print("Search for left genes")
        for voisin in  range(index - nb_voisins, index ):
            if voisin >= 0 :
                left_genes.append(synteny[voisin])
            else :
                left_genes.append("NO DATA")

        right_genes = []   
        print("Search for right genes")         
        for voisin in  range(index + 1, index + nb_voisins + 1):
            if voisin < len(synteny):
                right_genes.append(synteny[voisin])
            else :
                right_genes.append("NO DATA")

        logfile.write("Left genes of ["+prot_name+"] = ")
        for gid in left_genes:
            if gid == "NO DATA":
                logfile.write(gid +" ")
            else:
                logfile.write(gid +":"+dico_geneid_2genename[gid]+" ")
        logfile.write("\n")

        logfile.write("Right genes of ["+prot_name+"] = ")
        for gid in right_genes:
            if gid == "NO DATA":
                logfile.write(gid +" ")
            else:
                logfile.write(gid +":"+dico_geneid_2genename[gid]+" ")
        logfile.write("\n")

        left_families = genes_2_families(left_genes,dico_gene_2_proteins,dico_cluster)
        right_families = genes_2_families(right_genes,dico_gene_2_proteins,dico_cluster)

        logfile.write("Left families of ["+prot_name+"] = ")
        for fam in left_families:
            logfile.write(fam+" ")
        logfile.write("\n")

        logfile.write("Right families of ["+prot_name+"] = ")
        for fam in right_families:
            logfile.write(fam+" ")
        logfile.write("\n")

        outfile.write(prot_name)
        for family in left_families :
            outfile.write(";"+str(family))
        for family in right_families :
                outfile.write(";"+str(family))    
        outfile.write("\n")




#df.to_csv(args.output, sep=';')



