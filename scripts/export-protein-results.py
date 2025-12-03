#!/usr/bin/env python3
import argparse

from needle.blast import Results, group_matches
from needle.hits import hmm_clean, hmm_find_proteins
from needle.io import export_protein_hits

def main():
    parser = argparse.ArgumentParser(description="Export protein matches from BLAST TSV.")
    parser.add_argument("query_fasta", help="Path to query protein FASTA")
    parser.add_argument("target_fasta", help="Path to target genome FASTA")
    parser.add_argument("results_tsv", help="BLAST results TSV (NCBI headers)")
    parser.add_argument("hmm_dir", help="Directory of HMM profiles")
    parser.add_argument("genome_accession", help="Genome accession")
    parser.add_argument("output_dir", help="Path of output directory")
    args = parser.parse_args()

    res = Results(args.results_tsv, query_fasta_path=args.query_fasta, target_fasta_path=args.target_fasta)
    protein_matches = group_matches(res)

    # use HMM to find more fragments
    protein_matches = hmm_find_proteins(protein_matches, res, args.hmm_dir)

    protein_matches = [m for m in protein_matches if m.can_collate()]
    cleaned_protein_matches = hmm_clean(protein_matches, args.hmm_dir)

    export_protein_hits(
        args.genome_accession,
        cleaned_protein_matches,
        args.output_dir+"/proteins.tsv",
        args.output_dir+"/matches.tsv",
        args.output_dir+"/faa"
    )


if __name__ == "__main__":
    main()
