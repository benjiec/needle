#!/usr/bin/env python3
import os
import sys
import csv
import argparse
import subprocess
import itertools
import tempfile
from typing import List, Dict, Set

from needle.match import extract_subsequence_strand_sensitive, read_fasta_as_dict
from needle.hmm import parse_hmmsearch_domtbl
from needle.gff import parse_gff_to_hits
from needle.io import export_protein_hits
from defaults import DefaultPath


def read_ko_names(path: str) -> Dict[str, str]:
    name_by_accession = {}
    with open(path, "r") as f:
        for raw_line in f:
            if not raw_line:
                continue
            line = raw_line.rstrip("\n")
            if not line:
                continue
            header_content = line.strip()
            accession = header_content.split("\t")[0]
            name_by_accession[accession] = " ".join(header_content.split("\t")[1:])
    return name_by_accession


def read_fasta_sequence_names(path: str) -> Dict[str, str]:
    name_by_accession = {}
    with open(path, "r") as f:
        for raw_line in f:
            if not raw_line:
                continue
            line = raw_line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                header_content = line[1:].strip()
                accession = header_content.split(None, 1)[0]
                name_by_accession[accession] = header_content
    return name_by_accession


def run_cmd(cmd: List[str]) -> None:
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Command failed: {' '.join(cmd)}", file=sys.stderr)
        raise


def print_comparison(curated_protein_name, status, hmm_rows, gff_hit, needle_rows):
    print("")
    print(curated_protein_name)
    print(status)

    print("    found by hmmscan, on", gff_hit.target_accession, gff_hit.target_start, gff_hit.target_end)
    for row in hmm_rows:
        print("    hmm model", row["target_from"], row["target_to"], "matches", gff_hit.query_accession, row["query_from"], row["query_to"], row["evalue"])

    if not needle_rows:
        print("    DID NOT FIND NEEDLE RESULT")
    else:
        for row in needle_rows:
            print("    needle dna match", row["protein_hit_id"], row["target_start"], row["target_end"], "covering hmm model", row["query_start"], row["query_end"]) 
        for m in gff_hit.matches:
            print("    gff match", m.target_start, m.target_end)


def main():
    ap = argparse.ArgumentParser(description="Join hmmscan domtblout with GFF-derived proteins.")
    ap.add_argument("hmm_model", help="HMM database file (pressed already)")
    ap.add_argument("genome_accession")
    ap.add_argument("needle_match_file")
    args = ap.parse_args()

    hmm_file = DefaultPath.kegg_hmm(args.hmm_model)
    hmm_model = args.hmm_model
    hmm_name = read_ko_names("data/ko.tsv")[hmm_model]

    genome_accession = args.genome_accession
    ncbi_download_dir = DefaultPath.ncbi_genome_dir(genome_accession)

    faa_file = DefaultPath.ncbi_genome_faa(genome_accession)
    gff_file = DefaultPath.ncbi_genome_gff(genome_accession)
    fna_file = DefaultPath.ncbi_genome_fna(genome_accession)
    assert fna_file

    # Run hmmscan producing domtblout
    with tempfile.TemporaryDirectory() as tmpdir:
        domtbl_path = os.path.join(tmpdir, "out.domtbl")
        out_path = os.path.join(tmpdir, "out.txt")
        cmd = ["hmmscan", "--domtblout", domtbl_path, "-o", out_path, hmm_file, faa_file]
        run_cmd(cmd)

        # Parse domtbl into rows for TSV, using hmmSEARCH not hmmSCAN convention
        hmmscan_rows = parse_hmmsearch_domtbl(domtbl_path)

        # with hmmSCAN, the query and target coordinates are reversed from hmmSEARCH
        # hmmSCAN: scan proteins against HMM, so proteins are query, HMM is target
        # hmmSEARCH: search HMM against sequences, so HMM is query, sequence is target

        for row in hmmscan_rows:
            hmm_from = row["query_from"] 
            hmm_to = row["query_to"] 
            protein_from = row["target_from"] 
            protein_to = row["target_to"] 
            row["query_from"] = protein_from
            row["query_to"] = protein_to
            row["target_from"] = hmm_from
            row["target_to"] = hmm_to

        hmmscan_rows = [row for row in hmmscan_rows if row["evalue"] < 1e-3]

    print(hmm_model, hmm_name)
    print("received", len(hmmscan_rows), "matches from hmmscan")
    queries_from_hmmscan = list(set([m["query_name"] for m in hmmscan_rows]))

    # Parse GFF -> ProteinHit list
    protein_hits = parse_gff_to_hits(gff_file, protein_id_attr="protein_id")
    print("parsed protein hits from GFF, total", len(protein_hits))

    protein_seq_names = read_fasta_sequence_names(faa_file)

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
               row["query_start"] = int(row["query_start"])
               row["query_end"] = int(row["query_end"])
               row["target_start"] = int(row["target_start"])
               row["target_end"] = int(row["target_end"])
               needle_rows.append(row)
    print("needle results matching", hmm_model, "on same target accessions", len(needle_rows))

    def sorter(gff_hit):
        hmm_rows = [row for row in hmmscan_rows if row["query_name"] == gff_hit.query_accession]
        return max([row["target_to"] for row in hmm_rows]) - \
               min([row["target_from"] for row in hmm_rows])
    filtered_hits = sorted(filtered_hits, key=sorter, reverse=True)

    for gff_hit in filtered_hits:
        hmm_aa_len = 0
        hmm_rows = [row for row in hmmscan_rows if row["query_name"] == gff_hit.query_accession]
        for row in hmm_rows:
            hmm_aa_len += (row["target_to"] - row["target_from"] + 1)

        gff_dna_len = 0
        for m in gff_hit.matches:
            gff_dna_len += (abs(m.target_start - m.target_end) + 1)
        gff_target_left = min([min(m.target_start, m.target_end) for m in gff_hit.matches])
        gff_target_right = max([max(m.target_start, m.target_end) for m in gff_hit.matches])

        curated_protein_name = protein_seq_names[gff_hit.query_accession]

        needle_rows_on_target = [row for row in needle_rows if row["target_accession"] == gff_hit.target_accession]

        if len(needle_rows_on_target) == 0:
            status = "NOT FOUND"
            print_comparison(curated_protein_name, status, hmm_rows, gff_hit, None)

        else:
            keyf = lambda row: row["protein_hit_id"]
            needle_rows_on_target = sorted(needle_rows_on_target, key=keyf)

            for needle_protein_hit_id, needle_rows_for_hit in itertools.groupby(needle_rows_on_target, keyf):
                needle_rows_for_hit = list(needle_rows_for_hit)

                needle_target_left = min([min(row["target_start"], row["target_end"]) for row in needle_rows_for_hit])
                needle_target_right = max([max(row["target_start"], row["target_end"]) for row in needle_rows_for_hit])

                if needle_target_left > gff_target_right or \
                   needle_target_right < gff_target_left:
                    status = "NOT SAME LOCUS"
                    # print_comparison(curated_protein_name, status, hmm_rows, gff_hit, needle_rows_for_hit)

		    # not actually reporting this, because we are also not
		    # reporting other needle matches that were not discovered
		    # by hmmscan

                else:
                    needle_aa_len = 0
                    needle_dna_len = 0
                    for row in needle_rows_for_hit:
                        needle_aa_len += (row["query_end"] - row["query_start"] + 1)
                        needle_dna_len += (abs(row["target_start"] - row["target_end"]) + 1)

                    status = "FOUND"
                    substatus = []
                    if float(abs(needle_aa_len - hmm_aa_len) / hmm_aa_len) > 0.1:
                        substatus.append("LEN DIFF HMMSCAN")
                        if len(needle_rows_on_target) != len(gff_hit.matches):
                            substatus.append("MISSING EXON")
                    if len(substatus):
                        status += " ("+", ".join(substatus)+")"
                    print_comparison(curated_protein_name, status, hmm_rows, gff_hit, needle_rows_for_hit)


if __name__ == "__main__":
    main()
