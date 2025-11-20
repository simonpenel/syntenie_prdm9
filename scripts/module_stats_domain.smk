## Module for domain analysis on protein data
## Date : Decembre 2024
## Authors : L. Duret, A. Raimbault, S. Penel 
## Purpose: Search for domain sequences in the proteome.

# Rules
# -----

# ------------------------------------------------------
# hmm_search
# search the  HMM profile of the domain in the proteome.
# ------------------------------------------------------
rule hmm_search:
    """
    Proteome search using the HMMs.
    """
    input:
        # the hmm profile of the domain
        model=get_reference_file,
        # the proteome
        protein=pathGTDriftData + "genome_assembly/{accession}/annotation/protein.faa",
    output:
        # table of of per-sequence hits
        table=pathGTDriftData
        + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS + "hmm_search/tbl/{domain}",
        # table of per-domain hits
        domains=pathGTDriftData
        + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS + "hmm_search/domtbl/{domain}_domains",
    shell:
        "{RUNCMD} hmmsearch -E 1E-3 --domE 1E-3 --tblout {output.table} --domtblout {output.domains} --noali {input.model} {input.protein}"


# ------------------------------------------
# formating_hmm_sequence_hit
# write per-sequence hits in tabular format.
# ------------------------------------------
rule formating_hmm_sequence_hit:
    """
    Formating per-sequence hits from hmmsearch.
    """
    input:
        # per-sequence hits from hmmsearch
        pathGTDriftData
        + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS + "hmm_search/tbl/{domain}",
    output:
        # per-sequence hits in tabular format
        pathGTDriftData
        + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS + "hmm_search/tbl/{domain}_tabulated",
    script:
        "../utils/python/hmmsearch_parser.py"

# ----------------------------------------
# formating_hmm_domain_hit
# write per-domain hits in tabular format.
# ----------------------------------------
rule formating_hmm_domain_hit:
    """
    Formating per-domain hits from hmmsearch.
    """
    input:
        # per-domain hits from hmmsearch
        per_domain=pathGTDriftData
        + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS + "hmm_search/domtbl/{domain}_domains",
    output:
        # per-domain hits in tabular format in which overlapping zinc
        # finger domains are merged to create one big domain with 
        # multiple repetitions.
        domain_summary=pathGTDriftData
        + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS + "hmm_search/domtbl/{domain}_domains_summary",
    script:
        "../utils/python/domain_parser.py"

# ------------------------------------------------------------
# Function sending the accession number
# ------------------------------------------------------------
def accession_nb(wildcards):
    return wildcards.accession



# ------------------------------------------------------------
# summarize_hmm_results
# write results of hmm search  for each assembly and domain.
# Output format example:
# ;SeqID;SET Query;SET E-value;SET Score;Nb SET domains;SET domain start;SET domain end;Taxid
# ------------------------------------------------------------
rule summarize_hmm_results:
    """
    Creation of a summary table of hmm_search results.
    """
    input:
        # organisms_data file
        organisms_file=pathGTDriftData + "organisms_data",
        # path of all per-sequence hits in tabular format 
        domain_per_sequence_tabulated=pathGTDriftData + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS + "hmm_search/tbl/{domain}_tabulated",
        # path of all per-domain hits in tabular format with overlapping zinc finger domains                     
        domain_per_domain_summary=pathGTDriftData + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS + "hmm_search/domtbl/{domain}_domains_summary",        
    output:
        # domain protein statistics for each assembly.
        pathGTDriftData
        + "genome_assembly/{accession}/analyses/" + GENOME_RESULTS +"summary_hmmsearch_{accession}_{domain}.csv",
    params:    
        accession=accession_nb,  
    script:
        "../utils/python/table_domain_builder_single.py"

