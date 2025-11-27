import os
import subprocess
import tempfile
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

from .blast import NucMatch, ProteinMatch


@dataclass
class Candidate:
    split_k: Optional[int]  # for overlaps; None for gaps
    window_seq: str
    stitched_local: str


def order_matches_for_junctions(matches: List[NucMatch]) -> List[Tuple[NucMatch, NucMatch, int, int]]:
    if not matches:
        return []
    ordered = sorted(matches, key=lambda m: (m.query_start, m.query_end))
    pairs: List[Tuple[NucMatch, NucMatch, int, int]] = []
    for i in range(len(ordered) - 1):
        left = ordered[i]
        right = ordered[i + 1]
        overlap_len = max(0, left.query_end - right.query_start + 1)
        gap_len = max(0, right.query_start - left.query_end - 1)
        pairs.append((left, right, overlap_len, gap_len))
    return pairs


def _tail(seq: str, n: int) -> str:
    return seq[-n:] if n > 0 else ""


def _head(seq: str, n: int) -> str:
    return seq[:n] if n > 0 else ""


def generate_transition_candidates(
    left_aa: str,
    right_aa: str,
    overlap_len: int,
    gap_len: int,
    flank_window_max_length: int = 20,
) -> List[Candidate]:
    candidates: List[Candidate] = []
    if overlap_len > 0:
        left_len = len(left_aa)
        right_len = len(right_aa)
        for k in range(0, overlap_len + 1):
            # Assign k of the overlap to the left, and (overlap_len - k) to the right.
            # So trim (overlap_len - k) from left end and trim k from right start.
            left_prefix = left_aa[0: left_len - (overlap_len - k)]
            right_suffix = right_aa[k: right_len]
            stitched_local = left_prefix + right_suffix
            window_seq = _tail(left_prefix, flank_window_max_length) + _head(right_suffix, flank_window_max_length)
            candidates.append(Candidate(split_k=k, window_seq=window_seq, stitched_local=stitched_local))
        return candidates
    # Gap case ⇒ single candidate with X-fill
    fill = "X" * gap_len if gap_len > 0 else ""
    stitched_local = left_aa + fill + right_aa
    window_seq = _tail(left_aa, flank_window_max_length) + _head(right_aa, flank_window_max_length)
    candidates.append(Candidate(split_k=None, window_seq=window_seq, stitched_local=stitched_local))
    return candidates


def score_and_select_best_transition(
    candidates: List[Candidate],
    hmm_file_name: str,
) -> Candidate:
    if len(candidates) == 1:
        return candidates[0]
    try:
        with tempfile.TemporaryDirectory() as tmpdir:
            fasta_path = os.path.join(tmpdir, "cands.faa")
            domtbl_path = os.path.join(tmpdir, "out.domtbl")
            with open(fasta_path, "w") as f:
                for i, cand in enumerate(candidates):
                    f.write(f">cand_{i}\n{cand.window_seq}\n")
            cmd = ["hmmsearch", "--domtblout", domtbl_path, hmm_file_name, fasta_path]
            subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            best_idx = 0
            best_score = float("-inf")
            with open(domtbl_path, "r") as domf:
                for line in domf:
                    if not line or line.startswith("#"):
                        continue
                    parts = line.strip().split()
                    score = None
                    for idx in (13, 7, 8):
                        try:
                            score = float(parts[idx])
                            break
                        except Exception:
                            continue
                    if score is None:
                        score = -1e9
                    name = parts[0]
                    if name.startswith("cand_"):
                        idx = int(name.split("_")[1])
                        if score > best_score:
                            best_score = score
                            best_idx = idx
            return candidates[best_idx]
    except Exception:
        # Deterministic fallback
        best = max(
            enumerate(candidates),
            key=lambda t: (len(t[1].stitched_local), t[1].window_seq),
        )[0]
        return candidates[best]


def _aa_by_match(matches: List[NucMatch]) -> Dict[int, str]:
    mapping: Dict[int, str] = {}
    for m in matches:
        aa_full = m.target_sequence_translated()
        expected = m.query_end - m.query_start + 1
        mapping[id(m)] = aa_full[:expected]
    return mapping


def stitch_cleaned_sequence(
    matches: List[NucMatch],
    best_candidates_by_pair_index: Dict[int, Candidate],
    aa_by_match: Dict[int, str],
) -> str:
    if not matches:
        return ""
    ordered = sorted(matches, key=lambda m: (m.query_start, m.query_end))
    result = aa_by_match[id(ordered[0])]
    pairs = order_matches_for_junctions(ordered)
    for idx, (left, _right, _overlap, _gap) in enumerate(pairs):
        cand = best_candidates_by_pair_index[idx]
        left_len = len(aa_by_match[id(left)])
        result = result[: max(0, len(result) - left_len)] + cand.stitched_local
    return result


# Helper: clone a match
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
        on_reverse_strand=m.on_reverse_strand,
        matched_sequence=m.matched_sequence,
        query_sequence=m.query_sequence,
        target_sequence=m.target_sequence,
        target_sequence_upstream=m.target_sequence_upstream,
        target_sequence_downstream=m.target_sequence_downstream,
    )

# Helper: trim DNA sequence by AA counts
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

# Adjust coordinates and DNA per chosen candidate at each junction
def adjust_target_coordinates(left: NucMatch, right: NucMatch, cand: Candidate) -> Tuple[NucMatch, NucMatch]:
    nl = _clone(left)
    nr = _clone(right)
    if cand.split_k is None:
        # Gap: leave as-is
        return nl, nr
    k = max(0, int(cand.split_k))
    if k == 0:
        return nl, nr
    # Overlap: trim k AA from end of left; trim k AA from start of right (shift start by k).
    # Keep right end unchanged to preserve downstream gaps.
    nl.query_end = max(nl.query_start - 1, nl.query_end - k)
    # Only trim the left block; the right block remains unchanged
    # Do not alter right.query_start or right.target_start; overlap resolution is represented by trimming left only.
    # Adjust target genomic coordinates (5'->3') by 3 bases per trimmed AA for left only
    nl.target_end = max(nl.target_start, nl.target_end - 3 * k)
    # Trim target DNA accordingly for left only
    nl.target_sequence = _trim_dna_back(nl.target_sequence, k)
    # Right unchanged
    return nl, nr


def hmm_clean_protein(
    protein_match: ProteinMatch,
    hmm_file_name: str,
    flank_window_max_length: int = 20,
) -> ProteinMatch:

    # Compute AA per match and junction candidates
    aa_map = _aa_by_match(protein_match.matches)

    pairs = order_matches_for_junctions(protein_match.matches)
    selected: Dict[int, Candidate] = {}
    for idx, (left, right, overlap_len, gap_len) in enumerate(pairs):
        cands = generate_transition_candidates(
            aa_map[id(left)], aa_map[id(right)], overlap_len, gap_len, flank_window_max_length
        )
        best = cands[0] if len(cands) <= 1 else score_and_select_best_transition(cands, hmm_file_name)
        selected[idx] = best
        print("left", left.target_sequence_translated(), left.query_end,
              "right", right.target_sequence_translated(), right.query_start,
              "trim", best)

    # Stitch the final AA from original matches and chosen splits
    cleaned_aa = stitch_cleaned_sequence(protein_match.matches, selected, aa_map)

    new_matches: List[NucMatch] = []
    assert protein_match.matches
    ordered = sorted(protein_match.matches, key=lambda m: (m.query_start, m.query_end))
    current_left = _clone(ordered[0])
    for idx, (_left, right, _1, _2) in enumerate(pairs):
        selected_candidate = selected[idx]
        new_left, new_right = adjust_target_coordinates(current_left, right, selected_candidate)
        new_matches.append(new_left)
        current_left = new_right

        print("adjusting", idx,
              "new left", new_left.target_sequence_translated(), new_left.query_end,
              "new right", new_right.target_sequence_translated(), new_right.query_start,
              "trim", selected_candidate)

    if current_left.query_start <= current_left.query_end:
        new_matches.append(current_left)

    print(f"hmm_clean_protein expected_cleaned_sequence: {cleaned_aa}")
    cleaned_pm = ProteinMatch(
        target_id=protein_match.target_id,
        matches=new_matches,
        query_start=min(m.query_start for m in new_matches),
        query_end=max(m.query_end for m in new_matches),
        target_start=protein_match.target_start,
        target_end=protein_match.target_end,
        covers_start_to_end=protein_match.covers_start_to_end,
        likely_complete=protein_match.likely_complete,
        query_overlap=protein_match.query_overlap,
    )
    # Validate cleaned sequence matches the newly collated sequence from adjusted matches
    assert cleaned_pm.collated_protein_sequence == cleaned_aa, (
        f"Cleaned ProteinMatch collated sequence mismatch: "
        f"{cleaned_pm.collated_protein_sequence} != {cleaned_aa}"
    )
    return cleaned_pm


def hmm_clean(protein_matches: List[ProteinMatch], hmm_dir: str, flank_window_max_length: int = 20) -> List[ProteinMatch]:
    cleaned: List[ProteinMatch] = []
    for pm in protein_matches:
        q_acc = pm.matches[0].query_accession if pm.matches else None
        if not q_acc:
            cleaned.append(pm)
            continue
        hmm_path = os.path.join(hmm_dir, f"{q_acc}.hmm")
        if os.path.exists(hmm_path):
            cleaned.append(hmm_clean_protein(pm, hmm_path, flank_window_max_length))
        else:
            cleaned.append(pm)
    return cleaned


