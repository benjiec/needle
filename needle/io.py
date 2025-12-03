import os
import errno
from .blast import ProteinMatch
from .hits import hmmsearch_score
from typing import Dict, List, Optional


def write_protein_row(f, genome_accession: str, pm: ProteinMatch, evalue: Optional[float], score: Optional[float]) -> None:
    pid = pm.protein_hit_id
    row = [
        pid,
        pm.query_accession,
        pm.target_accession,
        genome_accession,
        "" if evalue is None else f"{evalue}",
        "" if score is None else f"{score}",
        pm.collated_protein_sequence,
    ]
    f.write("\t".join(row) + "\n")


def write_nucmatch_rows(f, pm: ProteinMatch) -> None:
    pid = pm.protein_hit_id
    for m in pm.matches:
        row = [
            pid,
            m.target_accession,
            str(m.target_start),
            str(m.target_end),
            str(m.query_start),
            str(m.query_end),
            str(m.target_sequence_translated())
        ]
        f.write("\t".join(row) + "\n")


def write_fasta_record(f, pm: ProteinMatch) -> None:
    pid = pm.protein_hit_id
    seq = pm.collated_protein_sequence
    f.write(f">{pid}\n")
    f.write(seq + "\n")


def _create_dirs_for_file(path: str) -> None:
    directory = os.path.dirname(path)
    if directory:
        try:
            os.makedirs(directory, exist_ok=True)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise

def create_tsv_if_missing_with_header(tsv_path: str, header_columns: List[str]) -> None:
    if not os.path.exists(tsv_path):
        _create_dirs_for_file(tsv_path)
        with open(tsv_path, "w") as f:
            f.write("\t".join(header_columns) + "\n")


def assert_tsv_header(tsv_path: str, header_columns: List[str]) -> None:
    expected = "\t".join(header_columns)
    with open(tsv_path, "r") as f:
        first_line = f.readline().rstrip("\n")
    assert first_line == expected, f"Unexpected TSV header in {tsv_path!r}: {first_line!r} != {expected!r}"


def export_protein_hits(
    genome_accession: str,
    protein_matches: List[ProteinMatch],
    proteins_tsv_path: str,
    nucmatches_tsv_path: str,
    proteins_fasta_dir: str,
) -> None:
    filtered = [pm for pm in protein_matches if pm.can_produce_single_sequence()]
    protein_header = [
        "protein_hit_id",
        "query_accession",
        "target_accession",
        "target_genome_accession",
        "hmmsearch_evalue",
        "hmmsearch_score",
        "collated_protein_sequence",
    ]
    nucmatch_header = [
        "protein_hit_id",
        "target_accession",
        "target_start",
        "target_end",
        "query_start",
        "query_end",
        "target_sequence"
    ]

    create_tsv_if_missing_with_header(proteins_tsv_path, protein_header)
    create_tsv_if_missing_with_header(nucmatches_tsv_path, nucmatch_header)
    assert_tsv_header(proteins_tsv_path, protein_header)
    assert_tsv_header(nucmatches_tsv_path, nucmatch_header)

    _create_dirs_for_file(os.path.join(proteins_fasta_dir, "dummy"))
    if not os.path.exists(proteins_fasta_dir):
        os.makedirs(proteins_fasta_dir, exist_ok=True)

    with open(proteins_tsv_path, "a") as f_prot, open(nucmatches_tsv_path, "a") as f_nuc:
        for pm in filtered:
            evalue: Optional[float] = None
            score: Optional[float] = None
            if pm.hmm_file:
                score, evalue = hmmsearch_score(pm.hmm_file, pm.collated_protein_sequence)
            else:
                print("no hmm file, skip hmmsearch")
            write_protein_row(f_prot, genome_accession, pm, evalue, score)
            write_nucmatch_rows(f_nuc, pm)

        # Write protein FASTA records grouped by query_accession into per-query files
        by_query: Dict[str, List[ProteinMatch]] = {}
        for pm in filtered:
            qa = pm.query_accession
            by_query.setdefault(qa, []).append(pm)
        for query_accession, pms in by_query.items():
            faa_path = os.path.join(proteins_fasta_dir, f"{query_accession}.faa")
            _create_dirs_for_file(faa_path)
            with open(faa_path, "a") as f_faa:
                for pm in pms:
                    write_fasta_record(f_faa, pm)
