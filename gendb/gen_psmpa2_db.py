#!/usr/bin/env python3

import argparse
import subprocess
import pandas as pd


class GenPSMPA2DB:
    """
    A class to build BLAST databases from 16S sequences and compute mean/median intensities
    merged with 16S copy number data, as well as handle lineage and copy data outputs.
    """
    def __init__(
        self,
        fasta_path,
        basedata_path,
        copydata_path,
        blastdb_prefix,
        mean_float_out_path,
        mean_int_out_path,
        median_float_out_path,
        median_int_out_path,
        lineage_out_path,
        copydata_out_path
    ):
        self.fasta_path = fasta_path
        self.basedata_path = basedata_path
        self.copydata_path = copydata_path
        self.blastdb_prefix = blastdb_prefix
        self.mean_float_out_path = mean_float_out_path
        self.mean_int_out_path = mean_int_out_path
        self.median_float_out_path = median_float_out_path
        self.median_int_out_path = median_int_out_path
        self.lineage_out_path = lineage_out_path
        self.copydata_out_path = copydata_out_path

    def build_blast_db(self):
        """
        Build a BLAST nucleotide database from the given FASTA file.
        """
        cmd = [
            'makeblastdb',
            '-in', self.fasta_path,
            '-dbtype', 'nucl',
            '-out', self.blastdb_prefix
        ]
        subprocess.run(cmd, check=True)
        print(f"[GenDB] BLAST database created with prefix: {self.blastdb_prefix}")

    def calc_bgc_statistics(self):
        """
        Compute mean and median BGC abundance from base data.
        Outputs both float and integer TSVs, with first column renamed to 'id'.
        """
        # Load base data
        df = pd.read_csv(self.basedata_path, sep='\t')
        cols = [
            '16S_ID', 'sum', 'PKSI', 'PKSother', 'NRPS', 'RiPPs',
            'Saccharides', 'Terpene', 'PKS-NRP_Hybrids', 'Others'
        ]
        # Check required columns
        missing = set(cols) - set(df.columns)
        if missing:
            raise KeyError(f"Missing required columns: {missing}")

        # Group by 16S_ID
        grouped = df.groupby('16S_ID')[cols[1:]]

        # Compute statistics
        mean_float = grouped.mean().reset_index()
        median_float = grouped.median().reset_index()

        mean_int = mean_float.copy()
        mean_int[cols[1:]] = mean_int[cols[1:]].round().astype(int)

        median_int = median_float.copy()
        median_int[cols[1:]] = median_int[cols[1:]].round().astype(int)

        # Rename first column
        for df_out in (mean_float, median_float, mean_int, median_int):
            df_out.rename(columns={'16S_ID': 'id'}, inplace=True)

        # Save outputs
        mean_float.to_csv(self.mean_float_out_path, sep='\t', index=False, compression='gzip')
        print(f"[GenDB] BGC mean (float) saved to: {self.mean_float_out_path}")

        mean_int.to_csv(self.mean_int_out_path, sep='\t', index=False, compression='gzip')
        print(f"[GenDB] BGC mean (int) saved to: {self.mean_int_out_path}")

        median_float.to_csv(self.median_float_out_path, sep='\t', index=False, compression='gzip')
        print(f"[GenDB] BGC median (float) saved to: {self.median_float_out_path}")

        median_int.to_csv(self.median_int_out_path, sep='\t', index=False, compression='gzip')
        print(f"[GenDB] BGC median (int) saved to: {self.median_int_out_path}")

    def export_lineage(self):
        """
        Read basedata table, extract lineage per 16S_ID, and write to file.
        Lineage format: p__phylum; g__genus; s__species
        """
        df = pd.read_csv(self.basedata_path, sep='\t')
        required = ['16S_ID', 'phylum', 'genus', 'species']
        missing = set(required) - set(df.columns)
        if missing:
            raise KeyError(f"Missing required columns for lineage: {missing}")

        records = []
        for id_val, group in df.groupby('16S_ID'):
            # Phylum: most frequent
            phylum_vals = group['phylum'].dropna().astype(str)
            phylum = phylum_vals.mode().iloc[0] if not phylum_vals.empty else ''

            # Genus: check consistency
            genus_vals = group['genus'].dropna().astype(str)
            unique_genus = genus_vals.unique()
            if len(unique_genus) > 1:
                print(f"16S_ID {id_val} group genus inconsistent: {', '.join(unique_genus)}")
            genus = genus_vals.mode().iloc[0] if not genus_vals.empty else ''

            # Species: take longest if multiple
            species_vals = group['species'].dropna().astype(str)
            unique_species = species_vals.unique()
            if len(unique_species) > 1:
                species = max(unique_species, key=len)
            else:
                species = unique_species[0] if len(unique_species) == 1 else ''

            lineage = f"{phylum}; {genus}; {species}"
            records.append({'id': id_val, 'lineage': lineage})

        df_lineage = pd.DataFrame(records)
        df_lineage.to_csv(self.lineage_out_path, sep='\t', index=False)
        print(f"[GenDB] Lineage data written to: {self.lineage_out_path}")


    def export_copydata(self):
        """
        Placeholder: process and export copy number data to file.
        """
        # TODO: implement copydata processing
        print(f"[GenDB] Copy number data would be written to: {self.copydata_out_path}")

    def run(self):
        """
        Execute the full pipeline.
        """
        self.build_blast_db()
        self.calc_bgc_statistics()
        self.export_lineage()
        self.export_copydata()


def main():
    parser = argparse.ArgumentParser(
        description="Build 16S BLAST DB and compute statistics with GenDB class."
    )
    parser.add_argument('--input16S', default="raw_data/16S.fasta",
                        help="16S FASTA file path.")
    parser.add_argument('--basedata', default="raw_data/basedata.tsv",
                        help="Base data TSV path.")
    parser.add_argument('--copydata', default="raw_data/16Scount.tsv",
                        help="16S copy number TSV path.")
    parser.add_argument('--blastdb_out', default="psmpa2/blast_db/rna",
                        help="BLAST DB prefix.")
    parser.add_argument('--mean_float_out', default="psmpa2/psmpa2_database_mean_float.tsv.gz",
                        help="Output gz TSV for mean float.")
    parser.add_argument('--mean_int_out', default="psmpa2/psmpa2_database_mean_int.tsv.gz",
                        help="Output gz TSV for mean int.")
    parser.add_argument('--median_float_out', default="psmpa2/psmpa2_database_median_float.tsv.gz",
                        help="Output gz TSV for median float.")
    parser.add_argument('--median_int_out', default="psmpa2/psmpa2_database_median_int.tsv.gz",
                        help="Output gz TSV for median int.")
    parser.add_argument('--lineage_out', default="psmpa2/psmpa2_database_lineage.tsv.gz",
                        help="Output path for lineage data.")
    parser.add_argument('--copydata_out', default="psmpa2/psmpa2_database_16S_count.tsv.gz",
                        help="Output path for processed copy data.")

    args = parser.parse_args()

    gendb = GenPSMPA2DB(
        fasta_path=args.input16S,
        basedata_path=args.basedata,
        copydata_path=args.copydata,
        blastdb_prefix=args.blastdb_out,
        mean_float_out_path=args.mean_float_out,
        mean_int_out_path=args.mean_int_out,
        median_float_out_path=args.median_float_out,
        median_int_out_path=args.median_int_out,
        lineage_out_path=args.lineage_out,
        copydata_out_path=args.copydata_out
    )
    gendb.run()


if __name__ == "__main__":
    main()