#!/usr/bin/env python


from os import path

project_dir = path.dirname(path.abspath(__file__))

# psmpa2 default files directory
default_psmpa2_dir = path.join(project_dir, "default_files", "psmpa2")

default_psmpa2_blast_database = path.join(default_psmpa2_dir, "blast_db", "rna")

default_psmpa2_database = {
    "mean_float": path.join(default_psmpa2_dir, 'psmpa2_database_mean_float.tsv.gz'),
    "mean_int": path.join(default_psmpa2_dir, 'psmpa2_database_mean_int.tsv.gz'),
    "median_float": path.join(default_psmpa2_dir, 'psmpa2_database_median_float.tsv.gz'),
    "median_int": path.join(default_psmpa2_dir, 'psmpa2_database_median_int.tsv.gz')
}

default_psmpa2_database_copy_number = path.join(default_psmpa2_dir, 'psmpa2_database_16S_count.tsv.gz')

default_psmpa2_database_lineage = path.join(default_psmpa2_dir, 'psmpa2_database_lineage.tsv.gz')

# psmpa1 default files directory
default_psmpa1_dir = path.join(project_dir, "default_files", "psmpa1")

default_psmpa1_ref_dir = path.join(default_psmpa1_dir, "pro_ref")

default_psmpa1_fasta = path.join(default_psmpa1_ref_dir, "pro_ref.fna")

default_psmpa1_tree = path.join(default_psmpa1_ref_dir, "pro_ref.tre")

default_psmpa1_hmm = path.join(default_psmpa1_ref_dir, "pro_ref.hmm")

default_psmpa1_model = path.join(default_psmpa1_ref_dir, "pro_ref.model")

default_psmpa1_raxml_info = path.join(default_psmpa1_ref_dir, "pro_ref.raxml_info")

# Inititalize default trait table files for hsp.
default_psmpa1_tables = {
    "16S": path.join(default_psmpa1_dir, "16S.tsv.gz"),
    "BGC": path.join(default_psmpa1_dir, "bgc.tsv.gz")
}

# psmpa-fungi default files directory
default_psmpa_fungi_dir = path.join(project_dir, "default_files", "psmpa_fungi")

default_psmpa_fungi_blast_database = path.join(default_psmpa_fungi_dir, "blast_db", "ref18S")

default_psmpa_fungi_database = {"default": path.join(default_psmpa_fungi_dir, "psmpa_fungi_database_default.tsv.gz")}
