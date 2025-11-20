
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
    Generate the candidates in fasta format and a summary of candidates in genomes.
    """
    #input:"concatenated_proteomes.faa"
    #input:"concatenated_proteomes.blastout"
    #input:"concatenated_proteomes.cluster"
    input:"concatenated_proteomes.cluster.sorted","concatenated_proteomes_id2species.csv",
        famfile_csv="concatenated_proteomes.cluster.sorted.csv",
        fam_nb_seq_spec="concatenated_proteomes_fam_nb_seq_species.csv",
        famfile_monogenic_csv="concatenated_proteomes.cluster.sorted.monogenic.csv",
    #input:"concatenated_proteomes.id2species"

rule concatenate_proteome:
    """
    Concatenate proteomes
    """    
    input:
        proteomes=expand(pathGTDriftData + "genome_assembly/{accession}/annotation/protein.faa",accession=ACCESSNB)
    output:
        concat="concatenated_proteomes.faa"
    shell:
        "cat {input.proteomes} > {output.concat}" 


rule diamond_db:
    """
    Build diamond db
    """    
    input:
        concat="concatenated_proteomes.faa"
    output:
        diamond_db="concatenated_proteomes.dmnd"
    shell:
        """
        diamond makedb --in {input.concat} -d  concatenated_proteomes -p 10
        """ 


rule diamond_blast:
    """
    Compare proteins:  all vs all
    """    
    input:
        diamond_db="concatenated_proteomes.dmnd",
        concat="concatenated_proteomes.faa"
    output:
        diamond_out="concatenated_proteomes.blastout"
    shell:
        """
        diamond blastp  -d  concatenated_proteomes -q {input.concat} -o {output.diamond_out} -p 20 --outfmt 6
        """

rule silix:
    """
    Clustering
    """    
    input:
        concat="concatenated_proteomes.faa",
        diamond="concatenated_proteomes.blastout"
    output:
        cluster="concatenated_proteomes.cluster"
    shell:
        """
        silix {input.concat}   {input.diamond} -f FAM > {output.cluster} || true
        """ 

rule sort_cluster:
    """
    Sorting
    """    
    input:
        cluster="concatenated_proteomes.cluster"
    output:
        sorted="concatenated_proteomes.cluster.sorted"
    shell:
        """
       sort {input.cluster}  > {output.sorted}
        """ 

rule get_species_for_id:
    """
    Species names
    """    
    input:
        concat="concatenated_proteomes.faa",
    output:
        id2spec="concatenated_proteomes_id2species.csv"
    shell:
        """
        grep ">"  {input.concat}  |sed -f ../scripts/sedfile  > {output.id2spec} 

        """
rule get_csv:
    """
    csv files
    """    
    input:
        famfile="concatenated_proteomes.cluster.sorted"
    output:
        famfile_csv="concatenated_proteomes.cluster.sorted.csv"
    shell:
        """
        awk ' {{OFS=";"}}; {{print $1,$2}}' {input.famfile} > {output.famfile_csv}
        """

rule fam_id_species:
    input:
        famfile_csv="concatenated_proteomes.cluster.sorted.csv",
        id2spec="concatenated_proteomes_id2species.csv"
    output:
        fam_nb_seq_spec="concatenated_proteomes_fam_nb_seq_species.csv"
    shell:
        """
        python3 ../scripts/merge_fam_species.py  -i {input.id2spec}  -f {input.famfile_csv} -o {output.fam_nb_seq_spec}
        """

rule select_monogenic:
    input:
        famfile_csv="concatenated_proteomes.cluster.sorted.csv",
        fam_nb_seq_spec="concatenated_proteomes_fam_nb_seq_species.csv"
    output:
        famfile_monogenic_csv="concatenated_proteomes.cluster.sorted.monogenic.csv",
    shell:
        """
        python3 ../scripts/get_monogenic.py -n {nb_spec} -s {input.fam_nb_seq_spec}  -f {input.famfile_csv} -o {output.famfile_monogenic_csv}
        """




