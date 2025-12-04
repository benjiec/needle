import hashlib
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple
from Bio.Seq import Seq

MAX_AA_OVERLAP = 20


@dataclass
class NucMatch:  # does not support matches across circular boundary
    query_accession: str
    target_accession: str
    query_start: int
    query_end: int
    target_start: int  # 1-based, 5' to 3' of gene, so target_start > target_end on reverse strand
    target_end: int    # 1-based, 5' to 3' of gene, so target_start > target_end on reverse strand
    e_value: float
    identity: float
    matched_sequence: Optional[str] = None

    query_sequence: Optional[str] = None   # 5' to 3'
    target_sequence: Optional[str] = None  # 5' to 3'

    @property
    def on_reverse_strand(self) -> bool:
        return self.target_start > self.target_end

    def target_sequence_translated(self) -> str:
        if not self.target_sequence or self.query_start > self.query_end:
            return ""
        dna = self.target_sequence.upper()
        usable_len = (len(dna) // 3) * 3
        if usable_len == 0:
            return ""
        dna = dna[:usable_len]
        return str(Seq(dna).translate(table="Standard", to_stop=False))


class NonlinearMatchException(Exception):
    pass


def order_matches_for_junctions(matches: List[NucMatch], max_overlap_len: int = MAX_AA_OVERLAP) -> List[Tuple[NucMatch, NucMatch, int, int]]:
    if not matches:
        return []

    ordered = sorted(matches, key=lambda m: (m.query_start, m.query_end))
    pairs: List[Tuple[NucMatch, NucMatch, int, int]] = []
    junctions: List[Tuple[int, int]] = []

    for i in range(len(ordered) - 1):
        left = ordered[i]
        right = ordered[i + 1]

        if right.query_end <= left.query_end or \
           right.query_start <= left.query_start:
            raise NonlinearMatchException("Found contained match")

        if left.on_reverse_strand != right.on_reverse_strand:
            raise NonlinearMatchException("Found fragments on different strands")

        if left.on_reverse_strand is False and \
           left.target_start > right.target_start:
            raise NonlinearMatchException("Consecutive protein fragments are reversed on the source DNA")

        if left.on_reverse_strand is True and \
           left.target_start < right.target_start:
            raise NonlinearMatchException("Consecutive protein fragments are reversed on the source DNA")

        overlap_len = max(0, left.query_end - right.query_start + 1)
        if overlap_len > max_overlap_len:
            raise NonlinearMatchException("Overlap too large, likely a different copy of the protein")

        gap_len = max(0, right.query_start - left.query_end - 1)
        pairs.append((left, right, overlap_len, gap_len))
        if gap_len:
            junctions.append((left.query_end, right.query_start))
        else:
            junctions.append((right.query_start, left.query_end))

    # none of the junctions should overlap, else we can't globally determine
    # the best sequence at each junction

    if len(junctions) > 1:
        cur_right = junctions[0][1]
        for left, right in junctions[1:]:
            if left <= cur_right:
                raise NonlinearMatchException("Junctions overlap")
            cur_right = right

    return pairs


@dataclass
class ProteinMatch:
    target_id: str
    matches: List[NucMatch]
    query_start: int
    query_end: int
    target_start: int  # 1-based, 5' to 3' of gene, so target_start > target_end on reverse strand
    target_end: int    # 1-based, 5' to 3' of gene, so target_start > target_end on reverse strand
    hmm_cleaned_protein_sequence: Optional[str] = None
    hmm_file: Optional[str] = None
    _protein_hit_id: Optional[str] = None

    @property
    def on_reverse_strand(self) -> bool:
        return self.target_start > self.target_end

    @staticmethod
    def can_collate_from_matches(matches) -> bool:
        try:
            pairs = order_matches_for_junctions(matches)
        except NonlinearMatchException:
            return False
        return True

    def can_collate(self) -> bool:
        return self.can_collate_from_matches(self.matches)

    @staticmethod
    def can_produce_single_sequence_from_matches(matches) -> bool:
        try:
            pairs = order_matches_for_junctions(matches)
        except NonlinearMatchException:
            return False
        for left, right, overlap, gaps in pairs:
            if overlap > 0:
                return False
        return True

    def can_produce_single_sequence(self) -> bool:
        return self.can_produce_single_sequence_from_matches(self.matches)

    @property
    def collated_protein_sequence(self) -> str:
        collated = ""
        if len(self.matches) < 2:
            return self.matches[0].target_sequence_translated()
        pairs = order_matches_for_junctions(self.matches)
        cur_left_aa = pairs[0][0].target_sequence_translated()
        for left, right, overlap, gaps in pairs:
            right_aa = right.target_sequence_translated()
            if gaps:
                new_s = cur_left_aa
                new_s += "X" * gaps
                collated += new_s
                cur_left_aa = right_aa
            else:
                new_s = cur_left_aa[0:len(cur_left_aa)-overlap]
                if overlap > 0:
                    new_s += "("+cur_left_aa[len(cur_left_aa)-overlap:]+"/"+right_aa[0:overlap]+")"
                collated += new_s
                cur_left_aa = right_aa[overlap:]
        collated += cur_left_aa
        return collated

    @property
    def protein_hit_id(self) -> str:
        if self._protein_hit_id is not None:
            return self._protein_hit_id
        assert self.matches
        ordered = sorted(self.matches, key=lambda m: (m.query_start, m.query_end))
        first = ordered[0]
        base_q = first.query_accession
        base_t = first.target_accession
        hasher = hashlib.sha1()
        for m in ordered:
            parts = [
                m.query_accession,
                m.target_accession,
                str(m.query_start),
                str(m.query_end),
                str(m.target_start),
                str(m.target_end),
            ]
            hasher.update("|".join(parts).encode("utf-8"))
        digest8 = hasher.hexdigest()[:8]
        self._protein_hit_id = f"{base_q}_{base_t}_{digest8}"
        return self._protein_hit_id

    @property
    def query_accession(self) -> str:
        assert self.matches
        first = sorted(self.matches, key=lambda m: (m.query_start, m.query_end))[0]
        return first.query_accession

    @property
    def target_accession(self) -> str:
        assert self.matches
        first = sorted(self.matches, key=lambda m: (m.query_start, m.query_end))[0]
        return first.target_accession


def group_matches(all_matches, max_intron_length: int = 10_000, max_aa_overlap: int = MAX_AA_OVERLAP) -> List[ProteinMatch]:
    """
    Group Match objects into ProteinMatch objects.
    - Groups by (query_accession, target_id)
    - Within a (query, target) group, further splits into clusters if adjacent
      matches on the target are separated by more than max_intron_length.
    """

    if not all_matches:
        return []

    # Helper to get interval on target chromosome (normalize orientation)
    def target_interval(m: NucMatch) -> (int, int):
        return (min(m.target_start, m.target_end), max(m.target_start, m.target_end))

    grouped: Dict[tuple, List[NucMatch]] = {}
    for m in all_matches:
        key = (m.query_accession, m.target_accession)
        grouped.setdefault(key, []).append(m)

    protein_matches: List[ProteinMatch] = []

    for (query_id, target_id), group in grouped.items():
        # Sort by target interval start to create distance-based clusters
        group_sorted_by_target = sorted(group, key=lambda m: target_interval(m)[0])

        clusters: List[List[NucMatch]] = []
        current_cluster: List[NucMatch] = []
        current_end: Optional[int] = None
        current_query_end: Optional[int] = None

        for m in group_sorted_by_target:
            start_t, end_t = target_interval(m)

            # New cluster
            if not current_cluster:
                current_cluster = [m]
                current_end = end_t
                current_query_end = m.query_end

            else:
                distance = start_t - (current_end or start_t)
                # Too far by target nuc distance, start new cluster
                if distance > max_intron_length:
                    clusters.append(current_cluster)
                    current_cluster = [m]
                    current_end = end_t
                    current_query_end = m.query_end

                # Query rewound, start new cluster
                elif m.query_start < current_query_end and current_query_end-m.query_start > max_aa_overlap:
                    clusters.append(current_cluster)
                    current_cluster = [m]
                    current_end = end_t
                    current_query_end = m.query_end

                # Add to cluster
                else:
                    current_cluster.append(m)
                    current_end = max(current_end or end_t, end_t)
                    current_query_end = m.query_end

        if current_cluster:
            clusters.append(current_cluster)

        # Build ProteinMatch objects for each cluster
        for cluster in clusters:
            # For boolean computations, sort by query_start
            cluster_by_query = sorted(cluster, key=lambda m: (m.query_start, m.query_end))

            # Aggregate min/max coordinates
            q_min = min(m.query_start for m in cluster)
            q_max = max(m.query_end for m in cluster)
            t_min = min(target_interval(m)[0] for m in cluster)
            t_max = max(target_interval(m)[1] for m in cluster)

            # Determine strand by majority orientation among matches:
            # reverse if target_start > target_end on most matches
            reverse_votes = sum(1 for m in cluster if m.on_reverse_strand)
            is_reverse = reverse_votes > (len(cluster) - reverse_votes)
            if is_reverse:
                pm_t_start, pm_t_end = t_max, t_min  # 5'->3' on reverse: higher coord to lower coord
            else:
                pm_t_start, pm_t_end = t_min, t_max

            protein_matches.append(
                ProteinMatch(
                    target_id=target_id,
                    matches=cluster_by_query,
                    query_start=q_min,
                    query_end=q_max,
                    target_start=pm_t_start,
                    target_end=pm_t_end
                )
            )

    return protein_matches


def extract_subsequence(full_sequence: Optional[str], start_1_based: int, end_1_based: int) -> Optional[str]:
    if full_sequence is None:
        return None
    if start_1_based <= 0 or end_1_based <= 0:
        return None
    # coordinates are 1-based inclusive; order can be reversed depending on alignment direction
    left = min(start_1_based, end_1_based)
    right = max(start_1_based, end_1_based)
    if left > len(full_sequence):
        return None
    # slice is exclusive of end; adjust for 1-based inclusive
    return full_sequence[left - 1 : min(right, len(full_sequence))]
