## Main snakemake for genewise analysis on genomic data
## Date : Decembre 2024
## Authors :


# Import python modules
# ---------------------
import os
import json

# Function to load JSON files
# ---------------------------
def load_json(file_path):
    with open(file_path, "r") as file:
        return json.load(file)

# Assign environment variables
# ----------------------------
globals().update(load_json("../environment_path.json"))


# Configuration
# -------------
configfile: "analyse.json"
configfile: "assemblies.json"


# List of assemblies
# ------------------
ACCESSNB = config["assembly_list"]

# Name of the resources directory. 
# --------------------------------
PROTEIN_RESOURCES_DIR_NAME  = config["protein_resources_dir_name"]

# -----------------------------------------------
# all : inputs define the to be files generated . 
# -----------------------------------------------       
            
rule all:
    input:
        expand("syntenie/{accession}/synteny.csv",accession=ACCESSNB), 
        
#rule get_proteins_chromosome_info:
#    input:
#        summary = pathGTDriftData
#        + "genome_assembly/{accession}/analyses/"+PROTEIN_RESOURCES_DIR_NAME + "whole_summary.csv",
#        gff = pathGTDriftData + "genome_assembly/{accession}/annotation/genomic.gff", 
#    output:
#        "syntenie/{accession}/whole_summary_with_chromosomes.csv",
#    shell:
#        """
#        python3 ../scripts/get_syntenies.py -i {input.summary}   -g {input.gff} -o {output}
#        """

rule get_syntenies:
    input:
        famfile_csv="concatenated_proteomes.cluster.sorted.monogenic.csv",
        fam_stat_csv="concatenated_proteomes_fam_nb_seq_species.csv",
        summary = pathGTDriftData
        + "genome_assembly/{accession}/analyses/"+PROTEIN_RESOURCES_DIR_NAME + "whole_summary.csv",
        #gff = pathGTDriftData + "genome_assembly/{accession}/annotation/genomic.gff", 
        positions = "analyse_gff/analyse.{accession}.csv"
    output:
        "syntenie/{accession}/synteny.csv",
    shell:
        """
        python3 ../scripts/get_syntenies.py  -i {input.summary}  -p {input.positions} -s {input.fam_stat_csv} -c {input.famfile_csv} -o {output}
        """        