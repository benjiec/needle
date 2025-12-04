import os
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple
from Bio.Seq import Seq
from .match import Match, ProteinHit, order_matches_for_junctions, extract_subsequence, extract_subsequence_strand_sensitive
from .hmm import hmmsearch


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
        # print("candidates were", sequences)
        return candidates[0]
    return candidates[best_idx]


def hmmsearch_score(hmm_file: str, protein_seq: str) -> Tuple[Optional[float], Optional[float]]:
    if not hmm_file or not protein_seq:
        return None, None
    sequences = [protein_seq]
    _, score, eval = hmmsearch_find_best_candidate(hmm_file, sequences)
    return score, eval


def aa_by_match(matches: List[Match]) -> Dict[int, str]:
    mapping: Dict[int, str] = {}
    for m in matches:
        aa_full = m.target_sequence_translated()
        mapping[id(m)] = aa_full
    return mapping


def stitch_cleaned_sequence(
    ordered_pairs: List[Tuple[Match, Match, int, int]],
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


def _clone(m: Match) -> Match:
    return Match(
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


def adjust_target_coordinates(left: Match, right: Match, cand: Candidate) -> Tuple[Match, Match]:
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
    protein_hit: ProteinHit,
    hmm_file_name: str,
    overlap_flanking_len: int = 20,
) -> ProteinHit:

    if len(protein_hit.matches) < 2:
        new_protein_hit = ProteinHit(
            matches=protein_hit.matches,
            query_start=protein_hit.query_start,
            query_end=protein_hit.query_end,
            target_start=protein_hit.target_start,
            target_end=protein_hit.target_end,
            hmm_file=hmm_file_name
        )
        return new_protein_hit

    """
    print("")
    print("cleaning", protein_hit.protein_hit_id, protein_hit.query_accession)
    """

    # Compute AA per match and junction candidates
    aa_map = aa_by_match(protein_hit.matches)
    pairs = order_matches_for_junctions(protein_hit.matches)

    selected: Dict[int, Candidate] = {}
    for idx, (left, right, overlap_len, gap_len) in enumerate(pairs):
        """
        print(idx, "left", left.query_start, left.query_end, left.target_start, left.target_end, aa_map[id(left)])
        print(idx, "right", right.query_start, right.query_end, right.target_start, right.target_end, aa_map[id(right)])
        print(idx, "overlap/gap", overlap_len, gap_len)
        """
        cands = generate_transition_candidates(
            aa_map[id(left)], aa_map[id(right)], overlap_len, gap_len, overlap_flanking_len
        )
        # print("choosing candidate for", left.query_start, left.query_end, " and ", right.query_start, right.query_end)
        best = cands[0] if len(cands) <= 1 else score_and_select_best_transition(cands, hmm_file_name)
        selected[idx] = best

    # Stitch the final AA from original matches and chosen splits
    cleaned_aa = stitch_cleaned_sequence(pairs, selected, aa_map)

    # Create new Match objects
    new_matches: List[Match] = []
    assert protein_hit.matches
    current_left = pairs[0][0]
    for idx, (_left, right, _1, _2) in enumerate(pairs):
        selected_candidate = selected[idx]
        new_left, new_right = adjust_target_coordinates(current_left, right, selected_candidate)
        new_matches.append(new_left)
        current_left = new_right
    new_matches.append(current_left)

    cleaned_pm = ProteinHit(
        matches=new_matches,
        query_start=min(m.query_start for m in new_matches),
        query_end=max(m.query_end for m in new_matches),
        target_start=protein_hit.target_start,
        target_end=protein_hit.target_end,
        hmm_file=hmm_file_name
    )
    # Validate cleaned sequence matches the newly collated sequence from adjusted matches
    assert cleaned_pm.collated_protein_sequence == cleaned_aa, (
        f"Cleaned ProteinHit collated sequence mismatch: "
        f"{cleaned_pm.collated_protein_sequence} != {cleaned_aa}"
    )
    return cleaned_pm


def hmm_clean(protein_hits: List[ProteinHit], hmm_dir: str, overlap_flanking_len: int = 20) -> List[ProteinHit]:
    cleaned: List[ProteinHit] = []
    for pm in protein_hits:
        hmm_path = os.path.join(hmm_dir, f"{pm.query_accession}.hmm")
        if os.path.exists(hmm_path):
            cleaned.append(hmm_clean_protein(pm, hmm_path, overlap_flanking_len))
        else:
            cleaned.append(pm)
    return cleaned


def hmmsearch_to_dna_coords(hmm_file, three_frame_translations):
    assert len(three_frame_translations) == 3

    sequences = [aa for dna_start, dna_end, aa in three_frame_translations]
    hmm_matches = hmmsearch(hmm_file, sequences)

    for m in hmm_matches:
        frame = int(m["target_name"][len("cand_"):])

    to_return = []
    for hmm_match in hmm_matches:
        assert hmm_match["target_name"].startswith("cand_")
        frame = int(hmm_match["target_name"][len("cand_"):])
        frame_dna_start = three_frame_translations[frame][0]
        frame_dna_end = three_frame_translations[frame][1]
        assert (abs(frame_dna_end-frame_dna_start)+1)%3 == 0

        aa_start = hmm_match["target_from"]
        aa_end = hmm_match["target_to"]
        aa = extract_subsequence(three_frame_translations[frame][2], aa_start, aa_end)
        hmm_match["matched_sequence"] = aa
        # print(hmm_match)

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

        # print("converted to dna coords", hmm_match)

    return to_return


def compute_three_frame_translations(full_seq, start, end):
    target_sequence = extract_subsequence_strand_sensitive(full_seq, start, end)
    if target_sequence is None:
        print("Cannot extract sequence using", len(full_seq), start, end)

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


def find_matches_at_locus(old_matches, full_seq, start, end, hmm_file, step=5000, force_extend=False):

    # print("HMM search space", start, end)

    translations = compute_three_frame_translations(full_seq, start, end)
    hmm_matches = hmmsearch_to_dna_coords(hmm_file, translations)

    new_matches = []
    for hmm_match in hmm_matches:
        target_sequence = extract_subsequence_strand_sensitive(full_seq, hmm_match["target_from"], hmm_match["target_to"])

        match = Match(
            query_accession=old_matches[0].query_accession,
            target_accession=old_matches[0].target_accession,
            query_start=hmm_match["query_from"],
            query_end=hmm_match["query_to"],
            target_start=hmm_match["target_from"],
            target_end=hmm_match["target_to"],
            e_value=hmm_match["evalue"],
            identity=None,
            target_sequence=target_sequence,
            matched_sequence=hmm_match["matched_sequence"]
        )
        assert match.matched_sequence == match.target_sequence_translated()
        new_matches.append(match)

    if not new_matches:
        return None

    if not ProteinHit.can_collate_from_matches(new_matches):
        return None

    old_query_start = min(m.query_start for m in old_matches)
    old_query_end = max(m.query_end for m in old_matches)
    old_nmatches = len(old_matches)

    new_query_start = min(m.query_start for m in new_matches)
    new_query_end = max(m.query_end for m in new_matches)
    new_nmatches = len(new_matches)

    if force_extend is False and \
       old_query_start == new_query_start and \
       old_query_end == new_query_end and \
       old_nmatches == new_nmatches:
        # print("Cannot further extend protein match")
        return None

    """
    print("Using HMM, found more protein fragments", new_matches[0].target_accession, new_matches[0].query_accession)
    print("Before", old_query_start, old_query_end, old_nmatches)
    print("Now", new_query_start, new_query_end, new_nmatches)
    for nm in new_matches:
        print("    ", nm.target_start, nm.target_end, nm.query_start, nm.query_end)
    """

    if end > start:
      if start > 1 or end < len(full_seq):
          more_matches = find_matches_at_locus(new_matches, full_seq, max(1, start-step), min(len(full_seq), end+step), hmm_file, step=step)
          return more_matches if more_matches else new_matches
    else:
      if end > 1 or start < len(full_seq):
          more_matches = find_matches_at_locus(new_matches, full_seq, min(len(full_seq), start+step), max(1, end-step), hmm_file, step=step)
          return more_matches if more_matches else new_matches

    return new_matches


def hmm_find_protein_around_locus(protein_hit, results, hmm_file):
    """
    Further refine protein match using hmmsearch, at the genomic locus
    """

    target_full_sequence = results._target_sequences_by_accession.get(protein_hit.target_accession, None)
    assert target_full_sequence is not None

    new_matches = find_matches_at_locus(protein_hit.matches, target_full_sequence, protein_hit.target_start, protein_hit.target_end, hmm_file, force_extend=True)
    if new_matches is None:
        return protein_hit

    query_sequence = results._query_sequences_by_accession.get(protein_hit.query_accession, None)
    for match in new_matches:
        match.query_sequence = extract_subsequence(query_sequence, match.query_start, match.query_end)

    new_pm = ProteinHit(
        matches=new_matches,
        query_start=min(m.query_start for m in new_matches),
        query_end=max(m.query_end for m in new_matches),
        target_start=min(m.target_start for m in new_matches) if protein_hit.target_start < protein_hit.target_end else max(m.target_start for m in new_matches),
        target_end=max(m.target_end for m in new_matches) if protein_hit.target_start < protein_hit.target_end else min(m.target_end for m in new_matches),
        hmm_file=hmm_file
    )


    return new_pm


def hmm_find_proteins(protein_hits, results, hmm_dir):
    new_protein_hits = []
    for pm in protein_hits:
        """
        print()
        print(pm.protein_hit_id)
        """
        hmm_path = os.path.join(hmm_dir, f"{pm.query_accession}.hmm")
        new_pm = hmm_find_protein_around_locus(pm, results, hmm_path)
        new_protein_hits.append(new_pm)
    return new_protein_hits
