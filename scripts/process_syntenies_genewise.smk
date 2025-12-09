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

# Name of the GENEWISE results directory. 
# --------------------------------
GENEWISE_RESULTS_DIR_NAME  = config["genewise_results_dir_name"]

# -----------------------------------------------
# all : inputs define the to be files generated . 
# -----------------------------------------------       
            
rule all:
    input:
        expand("syntenie_genewise/{accession}/synteny.csv",accession=ACCESSNB), 
        


rule get_syntenies:
    input:
        famfile_csv="concatenated_proteomes.cluster.sorted.monogenic.csv",
        #coordinates="/beegfs/banque/peneldb/gtdrift_template/pipeline/scripts/analyses/genewise_analysis/results/{accession}/Step3_genewise/genewise_prediction.txt",  
        #coordinates="/beegfs/banque/peneldb/gtdrift_template/pipeline/scripts/analyses/genewise_analysis/results/{accession}/Step3_genewise/gw.concat",
        fam_stat_csv="concatenated_proteomes_fam_nb_seq_species.csv",
        summary = pathGTDriftData
        + "genome_assembly/{accession}/analyses/"+GENEWISE_RESULTS_DIR_NAME + "whole_summary_genewise.csv",
        #gff = pathGTDriftData + "genome_assembly/{accession}/annotation/genomic.gff", 
        parsed = "analyse_gff/analyse.{accession}.csv"
    output:
        "syntenie_genewise/{accession}/synteny.csv",
    shell:
        """
        python3 ../scripts/get_syntenies_genewise.py  -i {input.summary}  -p {input.parsed} -s {input.fam_stat_csv} -c {input.famfile_csv} -o {output}
        """        