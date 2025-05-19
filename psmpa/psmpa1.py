#!/usr/bin/env python


import sys
from os import path
import pandas as pd
import shutil
import os
from psmpa.default import default_psmpa1_tables, default_psmpa1_ref_dir
from psmpa.place_seqs import identify_ref_files, place_seqs_pipeline
from psmpa.util import (make_output_dir, check_files_exist, read_fasta,
                        get_query_seq_id, read_seqabun, sample_bgc_calculate, create_output_folder, TemporaryDirectory)
from psmpa.wrap_hsp import castor_hsp_workflow


def psmpa1_pipeline(study_fasta,
                    input_table,
                    output_folder,
                    processes,
                    placement_tool,
                    min_align,
                    edge_exponent,
                    chunk_size,
                    calc_nsti,
                    hsp_method,
                    remove_intermediate,
                    force,
                    verbose):
    """Function that contains wrapper commands for full psmpa1 pipeline.

    Args:
        Descriptions of all of these input arguments/options are given in the psmpa1 script.

    Returns:
        Output folder with all result files.
    """

    # check and create the output folder
    if force:
        try:
            create_output_folder(output_folder)
        except:
            shutil.rmtree(output_folder)
            os.makedirs(output_folder)
    else:
        create_output_folder(output_folder)

    out_tree = path.join(output_folder, "out.tre")

    # Check that all input files exist.
    ref_msa, tree, hmm, model = identify_ref_files(default_psmpa1_ref_dir, placement_tool)
    files2check = [study_fasta, ref_msa, tree, hmm, model]
    if input_table:
        files2check.append(input_table)

    # This will throw an error if any input files are not found.
    check_files_exist(files2check)

    # Check that sequence names in FASTA overlap with input table.
    if input_table:
        check_overlapping_seqs(study_fasta, input_table, verbose)

    if verbose:
        print("Placing sequences onto reference tree", file=sys.stderr)

    # Define folders for intermediate files (unless --remove_intermediate set).
    # Run place_seqs.
    if remove_intermediate:
        with TemporaryDirectory() as temp_dir:
            place_seqs_pipeline(study_fasta=study_fasta,
                                placement_tool=placement_tool,
                                out_tree=out_tree,
                                threads=processes,
                                out_dir=temp_dir,
                                min_align=min_align,
                                chunk_size=chunk_size,
                                verbose=verbose)

    else:
        intermediate_dir = path.join(output_folder, "intermediate")
        make_output_dir(intermediate_dir)
        # place_seqs_intermediate = path.join(intermediate_dir, "place_seqs")
        place_seqs_pipeline(study_fasta=study_fasta,
                            placement_tool=placement_tool,
                            out_tree=out_tree,
                            threads=processes,
                            out_dir=intermediate_dir,
                            min_align=min_align,
                            chunk_size=chunk_size,
                            verbose=verbose)

    if verbose:
        print("Finished placing sequences on output tree: " + out_tree,
              file=sys.stderr)

    # Check that input files exist for hsp.
    check_files_exist([out_tree])

    hsp_out = {}
    for func, trait_table in default_psmpa1_tables.items():
        if not input_table and (func == '16S'):
            continue
        check_files_exist([trait_table])
        hsp_outfile = 'psmpa1_' + func + '_result.tsv.gz'
        # Run hsp for each function database.
        hsp_table, ci_table = castor_hsp_workflow(tree_path=out_tree,
                                                  trait_table_path=trait_table,
                                                  hsp_method=hsp_method,
                                                  edge_exponent=edge_exponent,
                                                  chunk_size=chunk_size,
                                                  calc_nsti=calc_nsti,
                                                  calc_ci=False,
                                                  check_input=False,
                                                  num_proc=processes,
                                                  ran_seed=42,
                                                  verbose=verbose)
        if func == '16S':
            hsp_table = pd.merge(get_query_seq_id(study_fasta), hsp_table, how='left',
                                 left_on=['qseqid'], right_on=['qseqid']).set_index(['qseqid'])
        if func == 'BGC':
            # hsp_table = hsp_table.loc[:, (hsp_table != 0).any(axis=0)]
            hsp_table = pd.merge(get_query_seq_id(study_fasta), hsp_table, how='left',
                                 left_on=['qseqid'], right_on=['qseqid']).set_index(['qseqid'])
            hsp_table = calc_bgc_sum(hsp_table)

        hsp_out[hsp_outfile] = hsp_table
    for fn, df in hsp_out.items():
        fp = path.join(output_folder, fn)
        df.to_csv(fp, sep='\t', index=True, compression='gzip')

    if input_table:
        bgc_table = hsp_out['psmpa1_BGC_result.tsv.gz'].fillna(0)
        copy_number = pd.merge(get_query_seq_id(study_fasta), hsp_out['psmpa1_16S_result.tsv.gz'], how='left',
                               left_on=['qseqid'], right_on=['qseqid']).set_index(['qseqid'])
        # print(copy_number.head())
        feature_table = read_seqabun(input_table)
        # print(feature_table.head())
        feature_table_amended = feature_table.div(copy_number['16S_rRNA_Count'], axis=0).fillna(0)
        # print(feature_table_amended.head())
        sample_result = sample_bgc_calculate(feature_table_amended, bgc_table)
        sample_result.to_csv(path.join(output_folder, 'psmpa1_sample_result.tsv.gz'), sep='\t', compression='gzip')


def check_overlapping_seqs(in_seq, in_tab, verbose):
    """
    Check that ASV ids overlap between the input FASTA and sequence
    abundance table. Will throw an error if none overlap and will otherwise
    print number of overlapping ids to STDERR. Also throw warning if input
    ASV table contains a column called taxonomy.

    Args:
        in_seq: Study fasta file.
        in_tab: Feature table file.
        verbose: Print out details as commands are run.

    Returns:
        Error or None.
    """

    FASTA_ASVs = set(read_fasta(in_seq).keys())

    in_table = read_seqabun(in_tab)

    table_ASVs = set(in_table.index.values)

    num_ASV_overlap = len(table_ASVs.intersection(FASTA_ASVs))

    if 'taxonomy' in in_table.columns:
        print("Warning - column named \"taxonomy\" in abundance table - if "
              "this corresponds to taxonomic labels this should be removed "
              "before running this pipeline.", file=sys.stderr)

    # Throw error if 0 ASVs overlap between the two files.
    if num_ASV_overlap == 0:
        sys.exit("Stopping - no ASV ids overlap between input FASTA and "
                 "sequence abundance table")

    # Otherwise print to STDERR how many ASVs overlap between the two files
    # if verbose set.
    if verbose:
        print(str(num_ASV_overlap) + " of " + str(len(table_ASVs)) +
              " sequence ids overlap between input table and FASTA.\n",
              file=sys.stderr)


def calc_bgc_sum(bgc_table):
    """Calculate the sum of BGC for each query sequence.

    Args:
        bgc_table: A dataframe of predicted BGC distribution without sum column.

    Returns:
        A dataframe inserted a sum column after "qseqid"
    """

    bgc_table["sum"] = bgc_table.apply(lambda x: x.sum(), axis=1)
    sum_col = bgc_table.pop('sum')
    bgc_table.insert(0, 'sum', sum_col)
    return bgc_table
