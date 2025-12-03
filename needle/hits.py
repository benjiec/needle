import os
import subprocess
import tempfile
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple
import csv
import io
import errno

from .blast import NucMatch, ProteinMatch, order_matches_for_junctions


@dataclass
class Candidate:
    assigned_overlap_to_left: Optional[int]  # for overlaps; None for gaps
    window_seq: str
    stitched: str
    left_trimmed: int
    right_kept: str


def generate_transition_candidates(
    left_aa: str,
    right_aa: str,
    overlap_len: int,
    gap_len: int,
    overlap_flanking_len: int = 20,
) -> List[Candidate]:

    candidates: List[Candidate] = []

    left_len = len(left_aa)
    right_len = len(right_aa)
    assert left_len >= overlap_len
    assert right_len >= overlap_len
    assert overlap_len == 0 or gap_len == 0

    left_window_start = max(left_len-overlap_len-overlap_flanking_len, 0)
    right_window_end_plus1 = min(overlap_len+overlap_flanking_len, right_len)

    for k in range(0, overlap_len + 1):
        # Assign k of the overlap to the left, and (overlap_len - k) to the right.

        left_trim = (overlap_len - k)
        left_prefix = left_aa[: left_len - left_trim]
        left_window = left_aa[left_window_start: left_len - overlap_len + k]
        right_suffix = right_aa[k: right_len]
        right_window = right_aa[k: right_window_end_plus1]

        if overlap_len == 0 and gap_len:
            gap = "X" * gap_len
            assigned_overlap_to_left = None
        else:
            gap = ""
            assigned_overlap_to_left = k

        stitched = left_prefix + gap + right_suffix
        window_seq = left_window + gap + right_window

        candidates.append(
            Candidate(assigned_overlap_to_left=assigned_overlap_to_left,
                      window_seq=window_seq,
                      stitched=stitched, left_trimmed=left_trim, right_kept=gap+right_suffix))

    return candidates


def run_command(cmd: str):
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def parse_hmmsearch_score_for_cand(domtbl_path):
    best_idx = None
    best_score = float("-inf")
    best_evalue = None

    expected_headers = "# target name        accession   tlen query name           accession   qlen   E-value  score"
    expected_eval_idx = 6
    expected_score_idx = 7
    assert expected_headers.replace("# target name", "target_name").replace("query name", "query_name").split()[expected_score_idx] == "score"
    assert expected_headers.replace("# target name", "target_name").replace("query name", "query_name").split()[expected_eval_idx] == "E-value"

    has_headers = False

    with open(domtbl_path, "r") as domf:
        for line in domf:
            if line.startswith(expected_headers):
                has_headers = True
            if not line or line.startswith("#") or has_headers is False:
                continue
            parts = line.strip().split()
            score = float(parts[expected_score_idx])
            name = parts[0]
            if name.startswith("cand_"):
                idx = int(name.split("_")[1])
                if score > best_score:
                    best_score = score
                    best_evalue = float(parts[expected_eval_idx])
                    best_idx = idx

    return best_idx, best_score, best_evalue

def hmmsearch(hmm_file_name, sequences, score=True):
    with tempfile.TemporaryDirectory() as tmpdir:
        fasta_path = os.path.join(tmpdir, "cands.faa")
        domtbl_path = os.path.join(tmpdir, "out.domtbl")
        with open(fasta_path, "w") as f:
            for i, cand in enumerate(sequences):
                f.write(f">cand_{i}\n{cand}\n")
        cmd = ["hmmsearch", "--domtblout", domtbl_path, hmm_file_name, fasta_path]
        run_command(cmd)
        if score is False:
            with open(domtbl_path, "r") as domf:
                return "".join(domf.readlines())
        return parse_hmmsearch_score_for_cand(domtbl_path)

def score_and_select_best_transition(
    candidates: List[Candidate],
    hmm_file_name: str,
) -> Candidate:

    if len(candidates) == 1:
        return candidates[0]
    sequences = [c.window_seq for c in candidates]
    best_idx, _1, _2 = hmmsearch(hmm_file_name, sequences)
    if best_idx is None:
        print("WARNING: cannot determine best candidate using hmmsearch, default to first transition")
        return candidates[0]
    return candidates[best_idx]

def get_hmmsearch_score_eval(hmm_file: str, protein_seq: str) -> Tuple[Optional[float], Optional[float]]:
    if not hmm_file or not protein_seq:
        return None, None
    sequences = [protein_seq]
    _, score, eval = hmmsearch(hmm_file, sequences)
    return score, eval


def aa_by_match(matches: List[NucMatch]) -> Dict[int, str]:
    mapping: Dict[int, str] = {}
    for m in matches:
        aa_full = m.target_sequence_translated()
        mapping[id(m)] = aa_full
    return mapping


def stitch_cleaned_sequence(
    ordered_pairs: List[Tuple[NucMatch, NucMatch, int, int]],
    best_candidates_by_pair_index: Dict[int, Candidate],
    aa_by_match: Dict[int, str],
) -> str:

    result = ""
    for idx, (left, _right, _overlap, _gap) in enumerate(ordered_pairs):
        cand = best_candidates_by_pair_index[idx]
        if idx == 0:
            result += cand.stitched
        else:
            result = result[: len(result)-cand.left_trimmed]
            result += cand.right_kept
    return result


def _clone(m: NucMatch) -> NucMatch:
    return NucMatch(
        query_accession=m.query_accession,
        target_accession=m.target_accession,
        query_start=m.query_start,
        query_end=m.query_end,
        target_start=m.target_start,
        target_end=m.target_end,
        e_value=m.e_value,
        identity=m.identity,
        matched_sequence=m.matched_sequence,
        query_sequence=m.query_sequence,
        target_sequence=m.target_sequence,
    )

def _trim_dna_front(dna: Optional[str], aa_count: int) -> Optional[str]:
    if dna is None or aa_count <= 0:
        return dna
    bases = 3 * aa_count
    if bases >= len(dna):
        return ""
    return dna[bases:]

def _trim_dna_back(dna: Optional[str], aa_count: int) -> Optional[str]:
    if dna is None or aa_count <= 0:
        return dna
    bases = 3 * aa_count
    if bases >= len(dna):
        return ""
    return dna[: len(dna) - bases]

def adjust_target_coordinates(left: NucMatch, right: NucMatch, cand: Candidate) -> Tuple[NucMatch, NucMatch]:
    nl = _clone(left)
    nr = _clone(right)

    # Gap: leave as-is
    if cand.assigned_overlap_to_left is None:
        return nl, nr

    nl.query_end -= cand.left_trimmed
    nl.target_end -= 3 * cand.left_trimmed
    nl.target_sequence = _trim_dna_back(nl.target_sequence, cand.left_trimmed)

    nr.query_start += cand.assigned_overlap_to_left
    nr.target_start += 3 * cand.assigned_overlap_to_left
    nr.target_sequence = _trim_dna_front(nr.target_sequence, cand.assigned_overlap_to_left)

    return nl, nr


def hmm_clean_protein(
    protein_match: ProteinMatch,
    hmm_file_name: str,
    overlap_flanking_len: int = 20,
) -> ProteinMatch:

    if len(protein_match.matches) < 2:
        new_protein_match = ProteinMatch(
            target_id=protein_match.target_id,
            matches=protein_match.matches,
            query_start=protein_match.query_start,
            query_end=protein_match.query_end,
            target_start=protein_match.target_start,
            target_end=protein_match.target_end,
            hmm_file=hmm_file_name
        )
        return new_protein_match

    # Compute AA per match and junction candidates
    aa_map = aa_by_match(protein_match.matches)
    pairs = order_matches_for_junctions(protein_match.matches)

    selected: Dict[int, Candidate] = {}
    for idx, (left, right, overlap_len, gap_len) in enumerate(pairs):
        cands = generate_transition_candidates(
            aa_map[id(left)], aa_map[id(right)], overlap_len, gap_len, overlap_flanking_len
        )
        best = cands[0] if len(cands) <= 1 else score_and_select_best_transition(cands, hmm_file_name)
        selected[idx] = best

    # Stitch the final AA from original matches and chosen splits
    cleaned_aa = stitch_cleaned_sequence(pairs, selected, aa_map)

    # Create new NucMatch objects
    new_matches: List[NucMatch] = []
    assert protein_match.matches
    current_left = pairs[0][0]
    for idx, (_left, right, _1, _2) in enumerate(pairs):
        selected_candidate = selected[idx]
        new_left, new_right = adjust_target_coordinates(current_left, right, selected_candidate)
        new_matches.append(new_left)
        current_left = new_right
    new_matches.append(current_left)

    cleaned_pm = ProteinMatch(
        target_id=protein_match.target_id,
        matches=new_matches,
        query_start=min(m.query_start for m in new_matches),
        query_end=max(m.query_end for m in new_matches),
        target_start=protein_match.target_start,
        target_end=protein_match.target_end,
        hmm_file=hmm_file_name
    )
    # Validate cleaned sequence matches the newly collated sequence from adjusted matches
    assert cleaned_pm.collated_protein_sequence == cleaned_aa, (
        f"Cleaned ProteinMatch collated sequence mismatch: "
        f"{cleaned_pm.collated_protein_sequence} != {cleaned_aa}"
    )
    return cleaned_pm


def hmm_clean(protein_matches: List[ProteinMatch], hmm_dir: str, overlap_flanking_len: int = 20) -> List[ProteinMatch]:
    cleaned: List[ProteinMatch] = []
    for pm in protein_matches:
        q_acc = pm.matches[0].query_accession if pm.matches else None
        if not q_acc:
            cleaned.append(pm)
            continue
        hmm_path = os.path.join(hmm_dir, f"{q_acc}.hmm")
        if os.path.exists(hmm_path):
            cleaned.append(hmm_clean_protein(pm, hmm_path, overlap_flanking_len))
        else:
            cleaned.append(pm)
    return cleaned


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
                score, evalue = get_hmmsearch_score_eval(pm.hmm_file, pm.collated_protein_sequence)
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
