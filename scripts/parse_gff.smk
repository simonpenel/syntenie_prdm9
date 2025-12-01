## snakemake for analysis the dna sequences of the SET domains 

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




# function to get the name of the reference alignment for a domain
# ----------------------------------------------------------------
def get_domain(wildcards):
    domain = wildcards.domain
    return domain

# function to get the name of the reference alignment for a domain
# ----------------------------------------------------------------
def get_reference(wildcards):
    fname = DOMAIN_REFERENCES.get(wildcards.domain, "")
    return fname



# Rules
# -----

# -----------------------------------------------
# all : inputs define the to be files generated .
# -----------------------------------------------
rule all:
    input:
        analyse = expand("analyse_gff/analyse.{accession}.csv",accession=ACCESSNB)    
# -------------------------------------------------------
# get_SET_dna
# analyse of the dna sequences of the zincfinger proteins
# -------------------------------------------------------
rule analyse_gff:
    """
    analyse gff
    """
    input:
        # gff file of the assembly
        gff = pathGTDriftData + "genome_assembly/{accession}/annotation/genomic.gff"
    output:
        results = "analyse_gff/analyse.{accession}.csv"
    shell:
        """
        python3 ../scripts/explore_gff.py  -a {wildcards.accession} -g {input.gff} -o {output.results}
        """
