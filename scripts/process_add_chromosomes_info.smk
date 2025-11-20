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


# Name of global results directory.
# The directory is located in pathGTDriftGlobalResults
# ----------------------------------------------------
GLOBAL_RESULTS = config["analyse_dir_name"]

# Name of genome specific results directory.
# The directory is located in genome_assembly/{accession}/analyses/
# -----------------------------------------------------------------
GENOME_RESULTS = config["analyse_dir_name"]




    

        
# -----------------------------------------------
# all : inputs define the to be files generated . 
# -----------------------------------------------       
            
rule all:
    input:
        protein_with_chromo=expand(pathGTDriftData
        + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS
        +"whole_summary_with_chromosomes.csv",accession=ACCESSNB), 
        results=pathGTDriftGlobalResults
        + GLOBAL_RESULTS + "results_with_chromosomes.csv"    
        
rule get_proteins_chromosome_info:
    input:
        summary = pathGTDriftData
        + "genome_assembly/{accession}/analyses/"+PROTEIN_RESOURCES_DIR_NAME + "whole_summary.csv",
        protein = pathGTDriftData + "genome_assembly/{accession}/annotation/protein.faa", 
        gff = pathGTDriftData + "genome_assembly/{accession}/annotation/genomic.gff", 
    output:
        pathGTDriftData
        + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS +"whole_summary_with_chromosomes.csv",
    shell:
        """
        python3 ../utils/python/get_chromosome_positions.py -i {input.summary}  -p {input.protein} -g {input.gff} -o {output}
        """

rule concatenate_chromosome_info:
    input:
         # Whole analyse summary
         # ---------------------
         whole_summary=expand(
            pathGTDriftData
            + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS 
            + "whole_summary_with_chromosomes.csv",
            accession=ACCESSNB)
    output:
        # Concatenation of assemblies results
        # -----------------------------------
        pathGTDriftGlobalResults
        + GLOBAL_RESULTS + "results_with_chromosomes.csv"
    script:
        "../utils/python/merge.py"
