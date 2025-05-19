import os
import shutil
from os import path
import numpy as np
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
from psmpa.default import default_psmpa_fungi_blast_database, default_psmpa_fungi_database
import pandas as pd
from psmpa.util import create_output_folder


def blast_runner(query, out, perc_identity=0, num_threads=4):
    """Run blast program in terminal.

    Use blast to match reference sequence with the highest identity for each query sequence.

    Args:
        query: Sequence file in fasta format.
        out: Blast format6 output.
        perc_identity: Blast identity threshold, 0-100.
        num_threads: The thread used when blast is run.

    Returns:
        stdout, stderr.
    """

    blast_cline = NcbiblastnCommandline(query=query, out=out, num_threads=num_threads, perc_identity=perc_identity,
                                        max_hsps=1, max_target_seqs=1, db=default_psmpa_fungi_blast_database, outfmt=6)
    stdout, stderr = blast_cline()
    return stdout, stderr


def amend_feature_table(feature_table, copy_number, blast_qseq_sseq_pair):
    """Correction of feature table using 16S copy number.

    Args:
        feature_table: A dataframe of feature table.
        copy_number: Copy number in blast_db.
        blast_qseq_sseq_pair: A dataframe of query sequence and target sequence pairs from blast result.

    Returns:
        A dataframe of amended feature table with the first column as index.
    """

    copy_number_database = pd.read_csv(copy_number, sep='\t', index_col=0, compression='gzip')
    qseq_copy_number = pd.merge(blast_qseq_sseq_pair, copy_number_database, how='left',
                                left_on=['sseqid'], right_index=True).drop(columns=['sseqid'])
    feature_table_amended = feature_table.div(qseq_copy_number['16S_rRNA_Count'], axis=0).fillna(0)
    return feature_table_amended, qseq_copy_number


class Matcher():
    """Match the blast result to the reference blast_db.
    The job of Matcher is to merge the dataframe that we need.
    First, we need to match the query sequence id to the blast results containing "sseqid" and "pident"
    because some query sequence ids will not appear in the blast result if blast doesn't match.
    Second, we need to match the lineage information containing "phylum", "genus", "species"
    to the previous merged dataframe.
    Last, we need to match the BGCs distribution to the previous merged dataframe.
    Attributes:
        qseqid: A list containing all the query sequence id.
        blast_result: A dataframe of blast result.
        psmpa2_database: A dataframe of psmpa2 blast_db.
    """

    def __init__(self, query, blast_out, method):
        """Inits Match Class with query seq id.
        Args:
            query: Fasta sequence file that input the program.
            blast_out: Blast result output file.
            method: Load the data processed by the specified method.
        """
        self.qseqid = self.__get_query_seq_id(query)
        self.blast_result = self.__read_blast_result(blast_out)
        self.psmpa2_database = self.__load_database(method)
        # self.lineage = self.__load_lineage(default_psmpa_fungi_database_lineage)

    def __get_query_seq_id(self, query):
        """Read the input fasta file and get the sequence id.
        Args:
            query: Input fasta file.
        Returns:
            A series of query sequence ids.
        """
        qseqid_ls = []
        for seq_record in SeqIO.parse(query, "fasta"):
            qseqid_ls.append(seq_record.id)
        qseqid = pd.Series(qseqid_ls, name='qseqid')
        return qseqid

    def __read_blast_result(self, blast_out):
        """Read the blast result output.
        Args:
            blast_out: blast output file.
        Returns:
            A dataframe containing blast result with the index as "qseqid".
        """
        col_names = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
                     'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
        blast_result = pd.read_csv(blast_out, sep='\t', header=None, index_col=0, names=col_names)
        return blast_result

    def __load_lineage(self, lineage_file):
        """Load the lineage information.
        Args:
            lineage_file: lineage tsv file containing bacterial information in blast_db.
        Returns:
            A dataframe of lineage of the bacteria in blast_db and the index is "id"
        """
        lineage = pd.read_csv(lineage_file, sep='\t', index_col=0, compression='gzip')
        return lineage

    def __load_database(self, method):
        """Load the psmpa2 blast_db.
        Load the psmpa2 blast_db containing the BGCs information and Set the first column as index.
        Args:
            blast_db: Direction of the target psmpa2 blast_db will be loaded.
        Returns:
            A dataframe of BGCs distribution with the index as "id"
        """
        # 分块，每一块是一个chunk，之后将chunk进行拼接
        database_chunk = pd.read_csv(default_psmpa_fungi_database[method], sep='\t', iterator=True,
                                     compression='gzip', index_col=0)
        loop = True
        chunk_size = 10000
        chunks = []
        while loop:
            try:
                chunk = database_chunk.get_chunk(chunk_size)
                chunks.append(chunk)
            except StopIteration:
                loop = False
        database = pd.concat(chunks)
        return database

    def qseqid_match_blastresult(self):
        """Match the query sequence id to the blast results.
        Returns:
            A dataframe containing "qseqid", "sseqid", "pident" columns.
        """
        need_col_names = ["sseqid", "pident"]
        qseqid_blastresult = pd.merge(self.qseqid, self.blast_result[need_col_names], how='left',
                                      left_on="qseqid", right_index=True).set_index(["qseqid"])
        return qseqid_blastresult

    def qseqid_blastresult_match_bgc(self):
        """match the bgc information to the previous merged dataframe.
        Returns:
            A dataframe containing "qseqid", "sseqid", "pident" columns and BGCs types.
        """
        qseqid_blastresult = self.qseqid_match_blastresult()
        qseqid_blastresult_bgc = pd.merge(qseqid_blastresult, self.psmpa2_database, how='left',
                                          left_on=['sseqid'], right_index=True)
        return qseqid_blastresult_bgc

    def qseqid_blastresult_match_lineage(self):
        """match the lineage information to the previous merged dataframe.
        Returns:
            A dataframe containing "qseqid", "sseqid", "pident", "lineage" columns.
        """
        qseqid_blastresult = self.qseqid_match_blastresult()
        qseqid_blastresult_lineage = pd.merge(qseqid_blastresult, self.lineage, how='left',
                                              left_on=['sseqid'], right_index=True)
        return qseqid_blastresult_lineage

    def qseqid_blastresult_lineage_match_bgc(self):
        """Match the BGCs distribution to the previous merged dataframe.
        Returns:
            A dataframe containing "qseqid", "sseqid", "pident", "lineage" columns and BGCs types.
        """
        qseqid_blastresult_lineage = self.qseqid_blastresult_match_lineage()
        qseqid_blastresult_lineage_bgc = pd.merge(qseqid_blastresult_lineage, self.psmpa2_database, how='left',
                                                  left_on=['sseqid'], right_index=True)
        return qseqid_blastresult_lineage_bgc


def psmpa_pipeline(study_fasta,
                   input_table,
                   output_folder,
                   blast_thread,
                   force,
                   method,
                   threshold,
                   verbose):
    """Function that contains wrapper commands for full psmpa2 pipeline.

    Descriptions of all of these input arguments/options are given in the
    psmpa2 script.

    Args:
          study_fasta: Sequence file in fasta format.
          input_table: Feature table in biom format.
          output_folder: Path to  output files.
          blast_thread: The thread used when blast is run.
          force: Overwrite the output folder if it exists.
          threshold: Set a blast threshold.
          verbose: Print out details as commands are running.

    Returns:
        A folder containing all result files.
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

    # verbose set
    if verbose:
        print('The output folder has been created.')

    # initialize variable
    df_to_save = {}  # {file_name: dataframe}

    # run blast
    blast_out = path.join(output_folder, 'blast_result.tsv')
    if threshold or blast_thread:
        blast_runner(query=study_fasta, out=blast_out, perc_identity=threshold, num_threads=blast_thread)
    else:
        blast_runner(query=study_fasta, out=blast_out)

    # verbose set
    if verbose:
        print("The BLAST program has been completed and the blast_result.tsv file has been generated.")

    # verbose set
    if verbose:
        print("Start analyzing blast results and predict BGC.")

    # match the corresponding information
    match = Matcher(query=study_fasta, blast_out=blast_out, method=method)
    psmpa_fungi_result = match.qseqid_blastresult_match_bgc().fillna(0)
    # psmpa_fungi_result = psmpa_fungi_result.loc[:, (psmpa_fungi_result != 0).any(axis=0)]
    bool_cond = psmpa_fungi_result['sseqid'] == 0
    psmpa_fungi_result[bool_cond] = psmpa_fungi_result[bool_cond].replace({0, '0'}, np.nan)
    df_to_save['psmpa_fungi_BGC_result.tsv.gz'] = psmpa_fungi_result

    # sample analysis
    # if input_table:
    #     # read feature table biom file
    #     feature_table = read_seqabun(input_table)
    #     # correction of feature table with 16S copy number
    #     feature_table_amended, qseq_copy_number = amend_feature_table(feature_table,
    #                                                                   default_psmpa_fungi_database_copy_number,
    #                                                                   match.qseqid_match_blastresult().drop(
    #                                                                       columns=['pident']))
    #     df_to_save['psmpa2_16S_result.tsv.gz'] = qseq_copy_number
    #     # select the data columns involved in the operation
    #     psmpa_fungi_result = psmpa_fungi_result.drop(columns=['sseqid', 'pident', 'lineage']).fillna(0)
    #     # calculate the total number of BGCs for each sample
    #     psmpa2_sample_result = sample_bgc_calculate(feature_table_amended, psmpa_fungi_result)
    #     df_to_save['psmpa2_sample_result.tsv.gz'] = psmpa2_sample_result

    # save all results
    for fn, df in df_to_save.items():
        fp = path.join(output_folder, fn)
        df.to_csv(fp, sep='\t', index=True, compression='gzip')

    # verbose set
    if verbose:
        query_seq_len = len(match.qseqid)
        print(f'{query_seq_len} sequences have been analyzed and the results have been saved.')
