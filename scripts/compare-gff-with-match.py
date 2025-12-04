#!/usr/bin/env python3
import os
import sys
import csv
import argparse
import subprocess
import tempfile
from pathlib import Path
from typing import List, Dict, Set

from needle.match import extract_subsequence_strand_sensitive, read_fasta_as_dict
from needle.hmm import parse_hmmsearch_domtbl
from needle.gff import parse_gff_to_hits
from needle.io import export_protein_hits


def run_cmd(cmd: List[str]) -> None:
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Command failed: {' '.join(cmd)}", file=sys.stderr)
        raise


def write_hmm_tsv(tsv_path: str, rows: List[Dict]) -> None:
    os.makedirs(os.path.dirname(tsv_path), exist_ok=True)
    if not rows:
        # Derive header from parser's typical keys
        header = ["target_name", "evalue", "score", "hmm_from", "hmm_to", "target_from", "target_to"]
        with open(tsv_path, "w") as f:
            f.write("\t".join(header) + "\n")
        return
    # Use keys from first row to preserve parser schema
    header = list(rows[0].keys())
    with open(tsv_path, "w") as f:
        f.write("\t".join(header) + "\n")
        for r in rows:
            f.write("\t".join(str(r.get(k, "")) for k in header) + "\n")


def main():
    ap = argparse.ArgumentParser(description="Join hmmscan domtblout with GFF-derived proteins.")
    ap.add_argument("hmm_file", help="HMM database file (pressed already)")
    ap.add_argument("ncbi_download_dir")
    ap.add_argument("needle_match_file")
    ap.add_argument("outdir", help="Output directory")
    args = ap.parse_args()

    hmm_file = os.path.abspath(args.hmm_file)
    hmm_model = Path(args.hmm_file).name.replace(".hmm","")

    genome_accession = Path(args.ncbi_download_dir).name
    faa_file = os.path.join(args.ncbi_download_dir, "protein.faa")
    gff_file = os.path.join(args.ncbi_download_dir, "genomic.gff")
    fna_file = None
    for pattern in ["*.fna", "*.fasta"]:
        matches = list(Path(args.ncbi_download_dir).glob(pattern))
        if matches:
            fna_file = str(matches[0])
    assert fna_file

    outdir = os.path.abspath(args.outdir)
    os.makedirs(outdir, exist_ok=True)

    # Run hmmscan producing domtblout
    with tempfile.TemporaryDirectory() as tmpdir:
        domtbl_path = os.path.join(tmpdir, "out.domtbl")
        out_path = os.path.join(tmpdir, "out.txt")
        cmd = ["hmmscan", "--domtblout", domtbl_path, "-o", out_path, hmm_file, faa_file]
        run_cmd(cmd)

        # Parse domtbl into rows for TSV, using hmmSEARCH not hmmSCAN convention
        domtbl_rows = parse_hmmsearch_domtbl(domtbl_path)

        # with hmmSCAN, the query and target coordinates are reversed from hmmSEARCH
        # hmmSCAN: scan proteins against HMM, so proteins are query, HMM is target
        # hmmSEARCH: search HMM against sequences, so HMM is query, sequence is target

        for row in domtbl_rows:
            hmm_from = row["query_from"] 
            hmm_to = row["query_to"] 
            protein_from = row["target_from"] 
            protein_to = row["target_to"] 
            row["query_from"] = protein_from
            row["query_to"] = protein_to
            row["target_from"] = hmm_from
            row["target_to"] = hmm_to

        domtbl_rows = [row for row in domtbl_rows if row["evalue"] < 1e-3]
        write_hmm_tsv(os.path.join(outdir, "hmmscan_output.tsv"), domtbl_rows)

    print("received", len(domtbl_rows), "matches from hmmscan")
    queries_from_hmmscan = list(set([m["query_name"] for m in domtbl_rows]))

    # Parse GFF -> ProteinHit list
    protein_hits = parse_gff_to_hits(gff_file, protein_id_attr="protein_id")
    print("parsed protein hits from GFF, total", len(protein_hits))

    # Filter to proteins that showed up in hmmscan output (by protein_id == query name)
    filtered_hits = [pm for pm in protein_hits if pm.query_accession in queries_from_hmmscan]
    print("filtered by hmmscan results to", len(filtered_hits))
    target_accessions = list(set([pm.target_accession for pm in filtered_hits]))

    # Load genomic fasta and set match target sequence
    genomic_fasta = read_fasta_as_dict(fna_file)
    for pm in filtered_hits:
        for m in pm.matches:
            m.target_sequence = extract_subsequence_strand_sensitive(genomic_fasta[m.target_accession], m.target_start, m.target_end)

    needle_rows = []
    with open(args.needle_match_file, newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
           if row["protein_hit_id"].startswith(hmm_model) and \
              row["target_accession"] in target_accessions:
               needle_rows.append(row)
    print("needle results matching", hmm_model, "on same target accessions", len(needle_rows))

    for gff_hit in filtered_hits:
        print(gff_hit.protein_hit_id, "found by hmmscan, on", gff_hit.target_accession, gff_hit.target_start, gff_hit.target_end)
        hmm_rows = [row for row in domtbl_rows if row["query_name"] == gff_hit.query_accession]
        assert hmm_rows
        for row in hmm_rows:
            print("    hmm model", row["target_from"], row["target_to"], "matches", gff_hit.query_accession, row["query_from"], row["query_to"], row["evalue"])
        needle_rows_on_target = [row for row in needle_rows if row["target_accession"] == gff_hit.target_accession]
        for row in needle_rows_on_target:
            print("    needle dna match", row["target_start"], row["target_end"], "covering hmm model", row["query_start"], row["query_end"]) 
        if len(needle_rows_on_target) == 0:
            print("    DID NOT FIND NEEDLE RESULT")

    if len(needle_rows):
        print("exporting needle rows")
        needle_tsv = os.path.join(outdir, "needle_matches.tsv")
        with open(needle_tsv, "w") as f:
            writer = csv.DictWriter(f, delimiter='\t', fieldnames=needle_rows[0].keys())
            writer.writeheader()
            for row in needle_rows:
               writer.writerow(row)

    print("exporting gff results")
    proteins_tsv = os.path.join(outdir, "gff_proteins.tsv")
    matches_tsv = os.path.join(outdir, "gff_matches.tsv")
    export_protein_hits(
        genome_accession=genome_accession,
        protein_hits=filtered_hits,
        proteins_tsv_path=proteins_tsv,
        nucmatches_tsv_path=matches_tsv,
        proteins_fasta_dir=None
    )


if __name__ == "__main__":
    main()
