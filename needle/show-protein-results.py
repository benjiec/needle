#!/usr/bin/env python3
import argparse

from blast import Results, group_matches


def main():
    parser = argparse.ArgumentParser(description="Show protein matches from BLAST TSV.")
    parser.add_argument("query_fasta", help="Path to query protein FASTA")
    parser.add_argument("target_fasta", help="Path to target genome FASTA")
    parser.add_argument("results_tsv", help="BLAST results TSV (NCBI headers)")
    args = parser.parse_args()

    res = Results(args.results_tsv, query_fasta_path=args.query_fasta, target_fasta_path=args.target_fasta)
    protein_matches = group_matches(res)

    for idx, pm in enumerate(protein_matches, start=1):
        print(f"== ProteinMatch {idx} ==")
        print(f"Target: {pm.target_id}  Query: {pm.matches[0].query_accession if pm.matches else '-'}")
        print(f"Query range: {pm.query_start}-{pm.query_end}")
        print(f"Target range (5'->3' on inferred strand): {pm.target_start}-{pm.target_end}")
        print(f"Covers start..end: {pm.covers_start_to_end}  Likely complete: {pm.likely_complete}  Overlap: {pm.query_overlap}")
        print()
        print(pm.pprint_target_protein_sequence())
        print()
        print("Collated:")
        print(pm.target_protein_sequence(res))
        print()


if __name__ == "__main__":
    main()


