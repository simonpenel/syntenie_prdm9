# version 1.1
# ============


import argparse
import time
import os
import sys
import pandas as pd
import re
import numpy as np
import subprocess
#from BCBio import GFF
from Bio import SeqIO
from Bio import codonalign
from Bio.Seq import Seq
from Bio.codonalign import CodonSeq
from Bio.SeqRecord import SeqRecord
parser = argparse.ArgumentParser()

parser.add_argument('-a', '--accession', type=str,required=True, help='accession')
parser.add_argument('-g', '--gff', type=str,required=True, help='path to gff file')
parser.add_argument('-o', '--output', type=str, required=True, help='output file path')

args = parser.parse_args()


#-------------------------------------------------------------------------
# Fonction processExonsGff
# fill 3 dictionnaries :
# dico_exons_pos: (contig, mrna) -> liste des (start,end,strand) des exons
# dico_cds_pos:   (contig, mrna) -> liste des (start,end,strand) des cds
# dico_prot:      protein name   -> liste des (contig,mrna)
#
# Notes on GFF3 Fields:
#
# Fields must be tab-separated. Also, all but the final field in each feature line must contain a value; "empty" columns should be denoted with a '.'
#
#     seqname - name of the chromosome or scaffold; chromosome names can be given with or without the 'chr' prefix. Important note: the seqname must be one used within Ensembl, i.e. a standard chromosome name or an Ensembl identifier such as a scaffold ID, without any additional content such as species or assembly. See the example GFF output below.
#     source - name of the program that generated this feature, or the data source (database or project name)
#     feature - feature type name, e.g. Gene, Variation, Similarity
#     start - Start position* of the feature, with sequence numbering starting at 1.
#     end - End position* of the feature, with sequence numbering starting at 1.
#     score - A floating point value.
#     strand - defined as + (forward) or - (reverse).
#     frame - One of '0', '1' or '2'. '0' indicates that the first base of the feature is the first base of a codon, '1' that the second base is the first base of a codon, and so on..
#     attribute - A semicolon-separated list of tag-value pairs, providing additional information about each feature.
#
# *- Both, the start and end position are included. For example, setting start-end to 1-2 describes two bases, the first and second base in the sequence.
#
#-------------------------------------------------------------------------
def processExonsGff(gff:str):
    with open(gff, 'r') as reader:
        print("Processing gff... (This may take several minutes)")
        print("Reading gff...")
        for line in reader:
            if line.startswith('#'):
                continue
            split_line = line.split('\t')
             # Processing gene
            if split_line[2] == 'gene':
                #print(split_line)
                gene_info = split_line[8]
                gene_infos = gene_info.split(';')
                gene_id = "none"
                gene_name = "none"
                for info in gene_infos:
                    #print("debug "+info)
                    if info.split("=")[0] == "Name":
                        #print("Name "+ info.split("=")[1] )
                        gene_name = info.split("=")[1]
                    if info.split("=")[0] == "ID":
                        #print("Name "+ info.split("=")[1] )
                        gene_id = info.split("=")[1]  
                        gene_id = re.sub(r'gene-', '', gene_id)

                    xrefs = info.split(",")

                    for xref in xrefs:
                        #print("test "+xref)
                        ref=xref.split(':')
                        #print("test2 "+ref[0])

                        if ref[0] == "Dbxref=GeneID" or ref[0] == "GeneID" :
                            gene_id = ref[1]
                            #print("gene id ")
                            #print(gene_id)
                #print("GENE NAME "+gene_name + " : "+ gene_id)
                if gene_id == "none":
                    print(line)
                    sys.exit(1)

                dico_gene[gene_id] = gene_name

                # xrefs = gene_info.split(';')[2].split(",")
                # for xref in xrefs:
                #     ref=xref.split(':')
                #     if ref[0] == "Dbxref=GeneID" or ref[0] == "GeneID" :
                #         gene_id = ref[1]
                #         print("gene id ")
                #         print(gene_id)
                #exit(1)        
            # Processing exons
            if split_line[2] == 'exon':
                contig=split_line[0]
                start=split_line[3]
                end=split_line[4]
                strand=split_line[6]
                exon_info = split_line[8]
                exon_infos = exon_info.split(';')
                parent = "none"
                for info in exon_infos:
                    infos = info.split("=")
                    if infos[0] == "Parent":
                        parent = infos[1]
                if (contig,parent) in dico_exons_pos:
                    dico_exons_pos[(contig,parent)].append([start,end,strand])
                else :
                    dico_exons_pos[(contig,parent)] = []
                    dico_exons_pos[(contig,parent)].append([start,end,strand])

            #processing cds
            if split_line[2] == 'CDS':
                #print("debug CDS")
                #print(split_line)
                contig=split_line[0]
                start=split_line[3]
                end=split_line[4]
                strand=split_line[6]
                frame=split_line[7]
                cds_gene_id = "none"
                prot_name = "none"
                cds_info   = split_line[8]
                cds_infos = cds_info.split(';')
                parent = "none"
                for info in cds_infos:
                    infos = info.split("=")
                    if infos[0] == "Parent":
                        parent = infos[1]
                flag_gene = 0
                flag_prot = 0
                xrefs = cds_info.split(';')[2].split(",")
                for xref in xrefs:
                    ref=xref.split(':')
                    if ref[0] == "Dbxref=GeneID" or ref[0] == "GeneID" :
                        cds_gene_id = ref[1]
                        #print("cds gene id ")
                        #print(cds_gene_id)
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
                if cds_gene_id == "none":
                    fwarn.write("WARNING: Unable to get gene id : ")
                    fwarn.write(cds_info)
                if prot_name == "none":
                    fwarn.write("WARNING: Unable to get protein name : ")
                    fwarn.write(cds_info)
                if prot_name != "none" and cds_gene_id != "none":
                    parent = dico_gene[cds_gene_id]
                    if prot_name in dico_prot:
                        dico_prot[prot_name].append([contig,parent])
                    else :
                        dico_prot[prot_name] = []
                        dico_prot[prot_name].append([contig,parent])
                    if (contig,parent) in dico_cds_pos:
                        dico_cds_pos[(contig,parent)].append([start,end,strand,frame])
                    else :
                        dico_cds_pos[(contig,parent)] = []
                        dico_cds_pos[(contig,parent)].append([start,end,strand,frame])



# output file
f = open(args.output, "w")
# output file
ferr = open(args.output+".ERROR", "w")
# log file
fwarn = open(args.output+".warnings", "w")
f.write("Accession;SeqID;Contig;Gene;Start;End;Strand\n")




# initialise dictionnaries
dico_exons_pos = {}     # Dico ["nom Contig","nom mrna"] => liste [start exon, end exon, strand]
dico_cds_pos = {}		# Dico ["nom Contig","nom mrna"] => liste [start cds, end cds, strand]
dico_prot = {}  		# Dico "nom protein" => liste ["nom Contig","nom mrna"]
dico_gene = {}          # Dico Gene id gene name

processExonsGff(args.gff)


for seqid in dico_prot.keys():
    #print("Test "+seqid)
    #contig_mrnas = dico_prot[seqid]
    contig_mrna_dico  = {}
     # get the uniques couples (contig, mrna) associated to the protein
    for contig_mrna in dico_prot[seqid]:
        if not (contig_mrna[0],contig_mrna[1]) in contig_mrna_dico:
            contig_mrna_dico[(contig_mrna[0],contig_mrna[1])] = "ok"
    contig_mrnas = contig_mrna_dico.keys()
    #print(contig_mrnas)
    for contig_mrna in contig_mrnas:
        #print(contig_mrna)
        contig = contig_mrna[0]
        mrna = contig_mrna[1]
        cds_pos = dico_cds_pos[(contig,mrna)]
        #print(len(cds_pos))
        strand = cds_pos[0][2]
        if strand == "+":
            first = cds_pos[0][0]
            last = cds_pos[len(cds_pos)-1][1]
        elif strand == "-":
            last = cds_pos[0][1]
            first = cds_pos[len(cds_pos)-1][0]
        else:
            sys.exit(1)
        f.write(args.accession+";"+seqid+";"+contig+";"+mrna+";"+str(first)+";"+str(last)+";"+strand+"\n")
        #print(cds_pos[0])
        #print(cds_pos[len(cds_pos)-1])
        #print(str(first)+":"+str(last))

