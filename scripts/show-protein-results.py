#!/usr/bin/env python3
import argparse

from needle.blast import Results, group_matches
from needle.hits import hmm_clean

def main():
    parser = argparse.ArgumentParser(description="Show protein matches from BLAST TSV.")
    parser.add_argument("query_fasta", help="Path to query protein FASTA")
    parser.add_argument("target_fasta", help="Path to target genome FASTA")
    parser.add_argument("results_tsv", help="BLAST results TSV (NCBI headers)")
    parser.add_argument("hmm_dir", help="Directory of HMM profiles")
    args = parser.parse_args()

    res = Results(args.results_tsv, query_fasta_path=args.query_fasta, target_fasta_path=args.target_fasta)
    protein_matches = group_matches(res)
    cleaned_protein_matches = hmm_clean(protein_matches, args.hmm_dir)

    for idx, pm in enumerate(cleaned_protein_matches, start=1):
        print(f"== ProteinMatch {idx} ==")
        print(f"Target: {pm.target_id}  Query: {pm.matches[0].query_accession if pm.matches else '-'}")
        print(f"Query range: {pm.query_start}-{pm.query_end}")
        print(f"Target range (5'->3' on inferred strand): {pm.target_start}-{pm.target_end}")
        print("Collated:")
        print(pm.collated_protein_sequence)
        print()


if __name__ == "__main__":
    main()
