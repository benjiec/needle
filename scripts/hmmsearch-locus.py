#!/usr/bin/env python3

import os
import argparse
from pathlib import Path

from needle.match import extract_subsequence_strand_sensitive, read_fasta_as_dict
from needle.hits import hmmsearch_to_dna_coords, compute_three_frame_translations

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("hmm_file")
    ap.add_argument("ncbi_download_dir")
    ap.add_argument("target_accession")
    ap.add_argument("target_start", type=int)
    ap.add_argument("target_end", type=int)
    args = ap.parse_args()

    hmm_file = os.path.abspath(args.hmm_file)
    genome_accession = Path(args.ncbi_download_dir).name

    fna_file = None
    for pattern in ["*.fna", "*.fasta"]:
        matches = list(Path(args.ncbi_download_dir).glob(pattern))
        if matches:
            fna_file = str(matches[0])
    assert fna_file

    genomic_fasta = read_fasta_as_dict(fna_file)
    assert args.target_accession in genomic_fasta

    translations = compute_three_frame_translations(genomic_fasta[args.target_accession], args.target_start, args.target_end)
    hmm_matches = hmmsearch_to_dna_coords(hmm_file, translations)
    for match in hmm_matches:
        print(match)


if __name__ == "__main__":
    main()
