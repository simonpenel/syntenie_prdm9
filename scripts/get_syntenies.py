import argparse
import time
import os
import sys
import pandas as pd
from Bio import SeqIO
 
parser = argparse.ArgumentParser()


parser.add_argument('-i', '--input', type=str, required=True, help='whole_summary.csv from protein analysis')
parser.add_argument('-g', '--gff', type=str,required=True, help='path to gff file, REQUIRED if type=prot')
parser.add_argument('-o', '--output', type=str, required=True, help='output file path')
parser.add_argument('-c', '--cluster', type=str, required=True, help='input cluster path')
parser.add_argument('-s', '--stats', type=str, required=True, help='input cluster stats')

args = parser.parse_args()

nb_voisins = 50

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


## Function to read gff files and get chromosome,
# start and end for locus.

def preProcessGff(gff:str):
    with open(gff, 'r') as reader:
        print("Preprocessing gff... (This may take several minutes)")        
        print("Reading gff...")
        for line in reader:
            if line.startswith('#'):
                continue
            if line.split('\t')[2] == 'gene' or line.split('\t')[2] == 'pseudogene':
            
                splitline = line.split('\t')
                chrom     = splitline[0]
                strand    = splitline[6]
                genepos   = [splitline[3], splitline[4]]
                gene_info   = splitline[8]
                pseudo = 0
                if line.split('\t')[2] == 'pseudogene':
                    pseudo = 1
                #identifiant = gene_info.split(';')[0].split("=")
                #if identifiant[0] == "ID" :
                #        gene_id = identifiant[1]
                #        dico_gene[gene_id] = [chrom,strand,genepos]   
                flag = 0                     
                xrefs = gene_info.split(';')[1].split(",")
                for xref in xrefs:
                    ref=xref.split(':')
                    if ref[0] == "Dbxref=GeneID" or ref[0] == "GeneID":
                        gene_id = ref[1]
                        dico_gene[gene_id] = [chrom,strand,genepos,pseudo]
                        if strand == "+":
                            if not chrom in dico_chromo_direct:
                                dico_chromo_direct[chrom]  = []
                                dico_chromo_direct[chrom].append([gene_id,genepos])
                            else :
                                dico_chromo_direct[chrom].append([gene_id,genepos])
                        else:
                            if not chrom in dico_chromo_reverse:
                                dico_chromo_reverse[chrom]  = []
                                dico_chromo_reverse[chrom].append([gene_id,genepos])
                            else :
                                dico_chromo_reverse[chrom].append([gene_id,genepos])

                        flag = 1
                if flag == 0 :
                    xrefs = gene_info.split(';')
                    for xref in xrefs:
                        ref=xref.split('=')
                        if ref[0] == "locus_tag" :
                            gene_id = ref[1].rstrip()
                            dico_gene[gene_id] = [chrom,strand,genepos,pseudo]
                            if strand == "+":
                                if not chrom in dico_chromo_direct:
                                    dico_chromo_direct[chrom]  = []
                                    dico_chromo_direct[chrom].append([gene_id,genepos])
                                else :
                                    dico_chromo_direct[chrom].append([gene_id,genepos])
                            else:
                                if not chrom in dico_chromo_reverse:
                                    dico_chromo_reverse[chrom]  = []
                                    dico_chromo_reverse[chrom].append([gene_id,genepos])
                                else :
                                    dico_chromo_reverse[chrom].append([gene_id,genepos])

                            flag = 1
                if flag == 0 :
                    print("debug2 "+gene_info)
                    
def processGff(gff:str):
    with open(gff, 'r') as reader:
        print("Processing gff... (This may take several minutes)")        
        print("Reading gff...")
        for line in reader:
            if line.startswith('#'):
                continue
            if line.split('\t')[2] == 'CDS':
                cds_gene_id = "none"
                prot_name = "none"
                splitline = line.split('\t')
                cds_info   = splitline[8]
                flag_gene = 0
                flag_prot = 0
                
                xrefs = cds_info.split(';')[2].split(",")
                for xref in xrefs:
                    ref=xref.split(':')
                    if ref[0] == "Dbxref=GeneID" or ref[0] == "GeneID" :
                        cds_gene_id = ref[1]
                        flag_gene = 1
                    if ref[0] == "GenBank"  or ref[0]== "Dbxref=GenBank":
                        prot_name = ref[1]
                        flag_prot = 1  
                    if ref[0] == "Genbank"  or ref[0] == "Dbxref=Genbank":
                        prot_name = ref[1]
                        flag_prot = 1  
                    if ref[0] == "NCBI_GP"  or ref[0] == "Dbxref=NCBI_GP":
                        prot_name = ref[1]
                        flag_prot = 1     
                if flag_gene == 0 :
                    xrefs = cds_info.split(';')
                    for xref in xrefs:
                        ref=xref.split('=')
                        if ref[0] == "locus_tag" :
                            cds_gene_id = ref[1]
                            flag_gene = 1
                            
                if flag_gene == 0 :
                    print("Warning : No gene")    
                if flag_prot == 0 :
                    print("Warning : No protein")                      
                if cds_gene_id == "none":
                    print("WARNING: Unable to get gene id : ")
                    print(cds_info)
                if prot_name == "none":
                    print("WARNING: Unable to get protein name : ")
                    print(cds_info)
                if prot_name != "none" and cds_gene_id != "none":
                    if not prot_name in dico_prot:
                        if  cds_gene_id in dico_gene:
                            dico_prot[prot_name] = dico_gene[cds_gene_id]
                            if not cds_gene_id in dico_gene_2_prot:
                                dico_gene_2_prot[cds_gene_id] = []
                                dico_gene_2_prot[cds_gene_id].append(prot_name)
                            else:
                                dico_gene_2_prot[cds_gene_id].append(prot_name)
                        else :
                            print("ERROR "+cds_gene_id + "(" + cds_info+ ")")    
dico_cluster = processCluster(args.cluster)
dico_stats = processStats(args.stats)
dico_chromo_direct = {}
dico_chromo_reverse = {}
dico_gene = {}                                                              
preProcessGff(args.gff)
dico_prot = {}     
dico_gene_2_prot = {}
processGff(args.gff)
dico_prot_position = {}
dico_chromo_direct_proteins = {}
dico_chromo_reverse_proteins = {}

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

logfile.write("DIRECT STRANDS:\n")
logfile.write("===============\n")
for chromo in dico_chromo_direct:
    logfile.write(chromo+": ")
    infos =  dico_chromo_direct[chromo]
    sorted_infos = sorted(infos, key=lambda d:int(d[1][0]))
    position = 0
    proteins = []

    for info in sorted_infos :
        gene = info[0]
        prot = dico_gene[gene]
        logfile.write("Gene [" + gene+"] ("+str(position)+")-> ")
        if gene in dico_gene_2_prot:
            prot_names = dico_gene_2_prot[gene] 
            proteins.append(prot_names)
            for prot_name in prot_names:
                logfile.write(prot_name+", ")
                dico_prot_position[prot_name] = [chromo, "+",position]
            position += 1

        else :
            logfile.write(" no protein ")
        logfile.write("; ")
    if position != len(proteins):
        print("last position  = "+str(position))
        print("nb of proteins  = "+str(len(proteins)))
        sys.exit("something is wrong")

    dico_chromo_direct_proteins[chromo]=  proteins  
    logfile.write("\n")

logfile.write("REVERSE STRANDS:\n")
logfile.write("===============\n")
for chromo in dico_chromo_reverse:
    logfile.write(chromo+":\n")
    infos =  dico_chromo_reverse[chromo]
    logfile.write("Avant\n")
    logfile.write(str(infos))
    logfile.write("\n")
    sorted_infos = sorted(infos, key=lambda d:int(d[1][0]))
    logfile.write("Apres\n")
    logfile.write(str(sorted_infos))
    logfile.write("\n")
    position = 0
    proteins = []

    for info in sorted_infos :
        gene = info[0]
        prot = dico_gene[gene]
        logfile.write("Gene [" + gene+"] ("+str(position)+")-> ")
        if gene in dico_gene_2_prot:
            prot_names = dico_gene_2_prot[gene] 
            proteins.append(prot_names)
            for prot_name in prot_names:
                logfile.write(prot_name+", ")
                dico_prot_position[prot_name] = [chromo, "-",position]
            position += 1

        else :
            logfile.write(" no protein ")
        logfile.write("; ")
    if position != len(proteins):
        print("last position  = "+str(position))
        print("nb of proteins  = "+str(len(proteins)))
        sys.exit("something is wrong")

    dico_chromo_reverse_proteins[chromo]=  proteins  
    logfile.write("\n")

df = pd.DataFrame() 

logfile.write("PROCESSING PROTEINS\n")
logfile.write("===================\n")
with open(args.input, 'r') as reader:
    lines = reader.readlines()[1:]
    for line in lines:
        splitline = line.split(';')
        prot_name = splitline[1]
        logfile.write("Processing "+ prot_name+"\n")
        if prot_name in dico_prot_position :
            left_families = []
            right_families = []
            logfile.write(str(dico_prot_position[prot_name]))
            logfile.write("\n")
            info_synt = dico_prot_position[prot_name]
            chromo = info_synt[0]
            strand = info_synt[1]
            pos = info_synt[2]
            if strand == "+":
                synteny = dico_chromo_direct_proteins[chromo]
                # check 
                flag_check = False
                for check in synteny[pos]:
                    if check == prot_name:
                        flag_check = True
                if not flag_check:
                    sys.exit("check error")
                for voisin in  range(pos - nb_voisins, pos -1):
                    logfile.write("Pos "+str(voisin) +" ")
                    if voisin > 0:
                        logfile.write(str(synteny[voisin]))
                        associated_families = []
                        logfile.write("Pos "+str(voisin) +" ")
                        logfile.write(str(synteny[voisin]))
                        for protein_voisine in synteny[voisin]:
                            print("Proteine voisine "+protein_voisine)
                            if protein_voisine in dico_cluster:
                                associated_family = dico_cluster[protein_voisine][0]
                                print("Famille =  "+dico_cluster[protein_voisine][0])
                            else:
                                print("Pas de famille")   
                                associated_family = "NO FAMILY"
                            if not associated_family in associated_families :
                                associated_families.append(associated_family)
                    else :
                        associated_families = ["NO DATA"]
                    left_families.append(associated_families)
                    logfile.write(str(associated_families))
                    logfile.write("\n")
                for voisin in  range(pos + 1, pos + nb_voisins):
                    logfile.write("Pos "+str(voisin) +" ")
                    if voisin < len(synteny):
                        logfile.write(str(synteny[voisin]))
                        associated_families = []
                        logfile.write("Pos "+str(voisin) +" ")
                        logfile.write(str(synteny[voisin]))
                        for protein_voisine in synteny[voisin]:
                            print("Proteine voisine "+protein_voisine)
                            if protein_voisine in dico_cluster:
                                associated_family = dico_cluster[protein_voisine][0]
                                print("Famille =  "+dico_cluster[protein_voisine][0])
                            else:
                                print("Pas de famille")   
                                associated_family = "NO FAMILY"
                            if not associated_family in associated_families :
                                associated_families.append(associated_family)
                    else :
                        associated_families = ["NO DATA"]
                    right_families.append(associated_families)
                    logfile.write(str(associated_families))
                    logfile.write("\n")                    

            if strand == "-":
                synteny = dico_chromo_reverse_proteins[chromo]
                # check 
                flag_check = False
                for check in synteny[pos]:
                    if check == prot_name:
                        flag_check = True
                if not flag_check:
                    sys.exit("check error")
                for voisin in  range(pos - nb_voisins, pos -1):
                    logfile.write("Pos "+str(voisin) +" ")
                    if voisin > 0:
                        logfile.write(str(synteny[voisin]))
                        associated_families = []
                        logfile.write("Pos "+str(voisin) +" ")
                        logfile.write(str(synteny[voisin]))
                        for protein_voisine in synteny[voisin]:
                            print("Proteine voisine "+protein_voisine)
                            if protein_voisine in dico_cluster:
                                associated_family = dico_cluster[protein_voisine][0]
                                print("Famille =  "+dico_cluster[protein_voisine][0])
                            else:
                                print("Pas de famille")   
                                associated_family = "NO FAMILY"
                            if not associated_family in associated_families :
                                associated_families.append(associated_family)
                    else :
                        associated_families = ["NO DATA"]
                    left_families.append(associated_families)
                    logfile.write(str(associated_families))
                    logfile.write("\n")
                for voisin in  range(pos + 1, pos + nb_voisins):
                    logfile.write("Pos "+str(voisin) +" ")
                    if voisin < len(synteny):
                        logfile.write(str(synteny[voisin]))
                        associated_families = []
                        logfile.write("Pos "+str(voisin) +" ")
                        logfile.write(str(synteny[voisin]))
                        for protein_voisine in synteny[voisin]:
                            print("Proteine voisine "+protein_voisine)
                            if protein_voisine in dico_cluster:
                                associated_family = dico_cluster[protein_voisine][0]
                                print("Famille =  "+dico_cluster[protein_voisine][0])
                            else:
                                print("Pas de famille")   
                                associated_family = "NO FAMILY"
                            if not associated_family in associated_families :
                                associated_families.append(associated_family)
                    else :
                        associated_families = ["NO DATA"]
                    right_families.append(associated_families)
                    logfile.write(str(associated_families))
                    logfile.write("\n")                  
            
            outfile.write(prot_name)
            
            for family in left_families :
                print("\nleft families")
                nbspec_max = 0
                nbseq_max = 0
                selected_family = family[0]
                for fam in family:
                    print("debug "+fam)
                    if fam != "NO FAMILY" and fam != "NO DATA" :
                        if fam in dico_stats:
                            stats = dico_stats[fam]
                            print(stats)
                            nbseq = stats[0]
                            nbspec = stats[1]
                            if nbspec > nbspec_max :
                                nbspec_max = nbspec
                                nbseq_max = nbseq
                                selected_family = fam
                            elif nbspec == nbspec_max and nbseq > nbseq_max :
                                selected_family = fam
                                nbseq_max = nbseq
                    print("Selected = "+ selected_family)
                outfile.write(";"+selected_family)
            for family in right_families :
                print("\nright families")
                nbspec_max = 0
                nbseq_max = 0
                selected_family = family[0]
                for fam in family:
                    print("debug "+fam)
                    if fam != "NO FAMILY" and fam != "NO DATA" :
                        if fam in dico_stats:
                            stats = dico_stats[fam]
                            print(stats)
                            nbseq = stats[0]
                            nbspec = stats[1]
                            if nbspec > nbspec_max :
                                nbspec_max = nbspec
                                nbseq_max = nbseq
                                selected_family = fam
                            elif nbspec == nbspec_max and nbseq > nbseq_max :
                                selected_family = fam
                                nbseq_max = nbseq
                    print("Selected = "+ selected_family)
                outfile.write(";"+selected_family) 
            outfile.write("\n")

            outfilecheck.write(prot_name)
            for family in left_families :
                outfilecheck.write(";"+str(family))
            for family in right_families :
                outfilecheck.write(";"+str(family))    
            outfilecheck.write("\n")

        else :
            logfile.write("not found\n")



#df.to_csv(args.output, sep=';')



