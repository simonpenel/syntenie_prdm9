
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
nb_spec = len(ACCESSNB)
print("Nb of species  = "+str(nb_spec))

# Rules
# -----

# -----------------------------------------------
# all : inputs define the to be files generated . 
# -----------------------------------------------
rule all:
    """
    sortie blast_compare
    """
    input:"concatenated_proteomes.blastout"

rule diamond_blast:
    """
    Compare proteins:  all vs all
    """    
    input:
        diamond_db="concatenated_proteomes.dmnd",
        concat="concatenated_proteomes.faa"
    output:
        diamond_out="concatenated_proteomes_diamond_fast.blastout"
    shell:
        """
        diamond blastp --faster  -d  concatenated_proteomes -q {input.concat} -o {output.diamond_out} -p 10 
        """

