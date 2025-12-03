import re
import os
import itertools
import subprocess
import tempfile
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple
from .blast import Results, NucMatch, ProteinMatch, order_matches_for_junctions
from Bio.Seq import Seq


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


def parse_hmmsearch_domtbl(domtbl_path):
    expected_header = "# target name  accession  tlen  query name  accession  qlen  E-value  score  bias  #  of  c-Evalue  i-Evalue  score  bias  from  to  from  to"
    idx_target = 0
    idx_eval = 6
    idx_score = 7
    idx_q_from = 15
    idx_q_to = 16
    idx_t_from = 17
    idx_t_to = 18

    # first, sanity check these indices
    expected_header_parts = re.split(r'\s\s+', expected_header)
    assert expected_header_parts[idx_eval] == "E-value"
    assert expected_header_parts[idx_score] == "score"
    assert expected_header_parts[idx_target] == "# target name"
    assert expected_header_parts[idx_q_from] == "from"
    assert expected_header_parts[idx_q_to] == "to"
    assert expected_header_parts[idx_t_from] == "from"
    assert expected_header_parts[idx_t_to] == "to"

    has_headers = False
    expected_header = " ".join(expected_header.split())

    matches = []
    with open(domtbl_path, "r") as domf:
        for line in domf:
            if " ".join(line.split()).startswith(expected_header):
                has_headers = True
            if not line or line.startswith("#") or has_headers is False:
                continue
            parts = line.strip().split()
            match = dict(
                target_name = parts[idx_target],
                evalue = float(parts[idx_eval]),
                score = float(parts[idx_score]),
                hmm_from = int(parts[idx_q_from]),
                hmm_to = int(parts[idx_q_to]),
                target_from = int(parts[idx_t_from]),
                target_to = int(parts[idx_t_to])
            )
            matches.append(match)

    assert has_headers
    return matches


def hmmsearch(hmm_file_name, sequences):
    with tempfile.TemporaryDirectory() as tmpdir:
        fasta_path = os.path.join(tmpdir, "cands.faa")
        domtbl_path = os.path.join(tmpdir, "out.domtbl")
        with open(fasta_path, "w") as f:
            for i, cand in enumerate(sequences):
                f.write(f">cand_{i}\n{cand}\n")
        cmd = ["hmmsearch", "--domtblout", domtbl_path, hmm_file_name, fasta_path]
        run_command(cmd)
        return parse_hmmsearch_domtbl(domtbl_path)


def hmmsearch_find_best_candidate(hmm_file_name, sequences):
    matches = hmmsearch(hmm_file_name, sequences)

    best_idx = None
    best_score = float("-inf")
    best_evalue = None

    for match in matches:
        name = match["target_name"]
        score = match["score"]
        evalue = match["evalue"]

        if name.startswith("cand_"):
            idx = int(name.split("_")[1])
            if score > best_score:
                best_score = score
                best_evalue = evalue
                best_idx = idx

    return best_idx, best_score, best_evalue


def score_and_select_best_transition(
    candidates: List[Candidate],
    hmm_file_name: str,
) -> Candidate:

    if len(candidates) == 1:
        return candidates[0]
    sequences = [c.window_seq for c in candidates]
    best_idx, _1, _2 = hmmsearch_find_best_candidate(hmm_file_name, sequences)
    if best_idx is None:
        print("WARNING: cannot determine best candidate using hmmsearch, default to first transition")
        return candidates[0]
    return candidates[best_idx]


def hmmsearch_score(hmm_file: str, protein_seq: str) -> Tuple[Optional[float], Optional[float]]:
    if not hmm_file or not protein_seq:
        return None, None
    sequences = [protein_seq]
    _, score, eval = hmmsearch_find_best_candidate(hmm_file, sequences)
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


def hmmsearch_to_dna_coords(hmm_file, three_frame_translations):
    assert len(three_frame_translations) == 3

    sequences = [aa for dna_start, dna_end, aa in three_frame_translations]
    hmm_matches = hmmsearch(hmm_file, sequences)

    to_return = []
    for hmm_match in hmm_matches:
        assert hmm_match["target_name"].startswith("cand_")
        frame = int(hmm_match["target_name"][len("cand_"):])
        frame_dna_start = three_frame_translations[frame][0]
        frame_dna_end = three_frame_translations[frame][1]
        assert (abs(frame_dna_end-frame_dna_start)+1)%3 == 0

        aa_start = hmm_match["target_from"]
        aa_end = hmm_match["target_to"]

        # HMM hit must be in fwd direction
        if aa_end < aa_start:
            continue

        if frame_dna_end > frame_dna_start: # fwd strand
            hmm_match["target_from"] = frame_dna_start+(aa_start-1)*3
            hmm_match["target_to"] = frame_dna_start+aa_end*3-1
            to_return.append(hmm_match)

        else: # rev strand 
            hmm_match["target_from"] = frame_dna_start-(aa_start-1)*3
            hmm_match["target_to"] = frame_dna_start-aa_end*3+1
            to_return.append(hmm_match)

    return to_return


def compute_three_frame_translations(full_seq, start, end):
    target_sequence = Results._extract_subsequence(full_seq, start, end)
    if start > end:
        target_sequence = Results._reverse_complement(target_sequence)

    translations = []
    for frame in range(3):
        trim_right = (len(target_sequence)-frame)%3
        if trim_right > 0:
          frame_sequence = target_sequence[frame:-trim_right]
        else:
          frame_sequence = target_sequence[frame:]
        assert len(frame_sequence) % 3 == 0
        aa = Seq(frame_sequence).translate(to_stop=False) # Translate entire sequence, including stops

        if end > start:  # fwd strand
            translations.append((start+frame, end-trim_right, str(aa)))
        else:  # rev strand
            translations.append((start-frame, end+trim_right, str(aa)))

    return translations


def find_matches_at_locus(old_matches, full_seq, start, end, hmm_file, step=2000):

    translations = compute_three_frame_translations(full_seq, start, end)
    hmm_matches = hmmsearch_to_dna_coords(hmm_file, translations)

    new_matches = []
    for hmm_match in hmm_matches:
        match = NucMatch(
            query_accession=old_matches[0].query_accession,
            target_accession=old_matches[0].target_accession,
            query_start=hmm_match["hmm_from"],
            query_end=hmm_match["hmm_to"],
            target_start=hmm_match["target_from"],
            target_end=hmm_match["target_to"],
            e_value=hmm_match["evalue"],
            identity=None,
        )
        new_matches.append(match)

    if not ProteinMatch.can_collate_from_matches(new_matches):
        return None

    old_query_start = min(m.query_start for m in old_matches)
    old_query_end = max(m.query_end for m in old_matches)
    old_nmatches = len(old_matches)

    new_query_start = min(m.query_start for m in new_matches)
    new_query_end = max(m.query_end for m in new_matches)
    new_nmatches = len(new_matches)

    if old_query_start == new_query_start and \
       old_query_end == new_query_end and \
       old_nmatches == new_nmatches:
        return None

    if end > start:
      if start > step or end+step <= len(full_seq):
          more_matches = find_matches_at_locus(new_matches, full_seq, max(0, start-step), min(len(full_seq), end+step), hmm_file)
          return more_matches if more_matches else new_matches
    else:
      if end > step or start+step <= len(full_seq):
          more_matches = find_matches_at_locus(new_matches, full_seq, min(len(full_seq), start+step), max(0, end-step), hmm_file)
          return more_matches if more_matches else new_matches

    return new_matches


def hmm_refine_protein(protein_match, results, hmm_file):
    """
    Further refine protein match using hmmsearch, at the genomic locus
    """

    target_full_sequence = results._target_sequences_by_accession.get(protein_match.target_accession, None)
    assert target_full_sequence is not None

    new_matches = find_matches_at_locus(protein_match.matches, target_full_sequence, protein_match.target_start, protein_match.target_end, hmm_file)
    if new_matches is None:
        return protein_match

    new_pm = ProteinMatch(
        target_id=protein_match.target_id,
        matches=new_matches,
        query_start=min(m.query_start for m in new_matches),
        query_end=max(m.query_end for m in new_matches),
        target_start=min(m.target_start for m in new_matches) if protein_match.target_start < protein_match.target_end else max(m.target_start for m in new_matches),
        target_end=max(m.target_end for m in new_matches) if protein_match.target_start < protein_match.target_end else min(m.target_end for m in new_matches),
        hmm_file=hmm_file
    )
    return new_pm
