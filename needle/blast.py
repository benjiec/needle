from __future__ import annotations

import csv
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Tuple
import hashlib
from Bio.Seq import Seq


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


def order_matches_for_junctions(matches: List[NucMatch]) -> List[Tuple[NucMatch, NucMatch, int, int]]:
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

        overlap_len = max(0, left.query_end - right.query_start + 1)
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


class Results:
    # NCBI-style canonical headers (preferred)
    H_QSEQID = "qseqid"
    H_SSEQID = "sseqid"
    H_EVALUE = "evalue"
    H_PIDENT = "pident"
    H_QSTART = "qstart"
    H_QEND = "qend"
    H_SSTART = "sstart"
    H_SEND = "send"
    H_SSEQ = "sseq"  # Matched_Sequence in legacy output

    # Producer header order (used by blast-genome.py). Keep minimal set we rely on.
    PRODUCER_HEADER = [
        H_QSEQID,
        H_SSEQID,
        H_EVALUE,
        H_PIDENT,
        H_QSTART,
        H_QEND,
        H_SSTART,
        H_SEND,
        H_SSEQ,
    ]

    # Raw BLAST outfmt fields (as configured in run_blast_search)
    RAW_OUTFMT_FIELDS = [
        "qseqid",
        "sseqid",
        "evalue",
        "pident",
        "qstart",
        "qend",
        "stitle",  # not included in output header, but present in raw
        "sseq",
        "sstart",
        "send",
    ]

    # Mapping from raw BLAST field names to our header names
    RAW_TO_HEADER = {
        "qseqid": H_QSEQID,
        "sseqid": H_SSEQID,
        "evalue": H_EVALUE,
        "pident": H_PIDENT,
        "qstart": H_QSTART,
        "qend": H_QEND,
        "sseq": H_SSEQ,
        "sstart": H_SSTART,
        "send": H_SEND,
    }
    # Inverse mapping: header -> raw key
    HEADER_TO_RAW = {v: k for k, v in RAW_TO_HEADER.items()}

    def __init__(
        self,
        results_tsv_path: str,
        query_fasta_path: Optional[str] = None,
        target_fasta_path: Optional[str] = None
    ) -> None:
        self._results_tsv_path = results_tsv_path
        self._query_fasta_path = query_fasta_path
        self._target_fasta_path = target_fasta_path
        self._matches: Optional[List[NucMatch]] = None

        self._query_sequences_by_accession: Optional[Dict[str, str]] = None
        self._target_sequences_by_accession: Optional[Dict[str, str]] = None

    def matches(self) -> List[NucMatch]:
        if self._matches is None:
            self._parse_once()
        # mypy: _matches is now set
        return self._matches or []

    # ----- Internal helpers -----

    def _parse_once(self) -> None:
        if self._query_fasta_path:
            self._query_sequences_by_accession = self._read_fasta_as_dict(self._query_fasta_path)
        if self._target_fasta_path:
            self._target_sequences_by_accession = self._read_fasta_as_dict(self._target_fasta_path)

        with open(self._results_tsv_path, "r") as tsv_file:
            reader = csv.reader(tsv_file, delimiter="\t")
            header_row = next(reader, None)
            if header_row is None:
                self._matches = []
                return
            header_index = self._build_header_index(header_row)

            matches: List[NucMatch] = []
            for row in reader:
                if not row or all(not cell for cell in row):
                    continue
                match = self._row_to_match(row, header_index)
                # Attach sequences if requested and available
                if match is not None:
                    if self._query_sequences_by_accession is not None:
                        match.query_sequence = self._extract_subsequence(
                            self._query_sequences_by_accession.get(match.query_accession, None),
                            match.query_start,
                            match.query_end,
                        )
                    if self._target_sequences_by_accession is not None:
                        seq = self._extract_subsequence(
                            self._target_sequences_by_accession.get(match.target_accession, None),
                            match.target_start,
                            match.target_end,
                        )
                        if seq is not None and match.on_reverse_strand:
                            seq = self._reverse_complement(seq)
                        match.target_sequence = seq
                    matches.append(match)

            self._matches = matches

    def _build_header_index(self, header_row: List[str]) -> Dict[str, int]:
        # Normalize header names by stripping whitespace
        normalized = [h.strip() for h in header_row]
        index: Dict[str, int] = {name: i for i, name in enumerate(normalized)}

        # Build a mapping for canonical names only; sseq is optional
        mapping: Dict[str, Optional[int]] = {
            self.H_QSEQID: index.get(self.H_QSEQID),
            self.H_SSEQID: index.get(self.H_SSEQID),
            self.H_EVALUE: index.get(self.H_EVALUE),
            self.H_PIDENT: index.get(self.H_PIDENT),
            self.H_QSTART: index.get(self.H_QSTART),
            self.H_QEND: index.get(self.H_QEND),
            self.H_SSTART: index.get(self.H_SSTART),
            self.H_SEND: index.get(self.H_SEND),
            self.H_SSEQ: index.get(self.H_SSEQ),
        }

        # Ensure required fields exist
        required = [
            self.H_QSEQID,
            self.H_SSEQID,
            self.H_EVALUE,
            self.H_PIDENT,
            self.H_QSTART,
            self.H_QEND,
            self.H_SSTART,
            self.H_SEND,
        ]
        missing = [key for key in required if mapping.get(key) is None]
        if missing:
            raise ValueError(f"Missing required columns in results TSV: {', '.join(missing)}")

        # Convert Optional[int] to int where present
        return {k: v for k, v in mapping.items() if v is not None}

    def _row_to_match(self, row: List[str], header_index: Dict[str, int]) -> Optional[NucMatch]:
        try:
            qacc = row[header_index[self.H_QSEQID]].strip()
            sacc = row[header_index[self.H_SSEQID]].strip()
            evalue_str = row[header_index[self.H_EVALUE]].strip()
            pident_str = row[header_index[self.H_PIDENT]].strip()
            qstart_str = row[header_index[self.H_QSTART]].strip()
            qend_str = row[header_index[self.H_QEND]].strip()
            sstart_str = row[header_index[self.H_SSTART]].strip()
            send_str = row[header_index[self.H_SEND]].strip()

            matched_seq = None
            if self.H_SSEQ in header_index and header_index[self.H_SSEQ] < len(row):
                cell = row[header_index[self.H_SSEQ]].strip()
                matched_seq = cell if cell != "" else None

            qstart = int(qstart_str)
            qend = int(qend_str)
            if qstart > qend:
                raise ValueError(f"qstart ({qstart}) must be <= qend ({qend}) in results TSV row: {row}")

            sstart = int(sstart_str)
            send = int(send_str)

            match = NucMatch(
                query_accession=qacc,
                target_accession=sacc,
                query_start=qstart,
                query_end=qend,
                target_start=sstart,
                target_end=send,
                e_value=float(evalue_str.replace(",", "")),
                identity=float(pident_str.replace(",", "")),
                matched_sequence=matched_seq,
            )
            return match
        except (IndexError, ValueError) as exc:
            # Skip malformed rows; callers generally prefer partial results over failure
            # Narrow exception types only
            raise ValueError(f"Malformed row in results TSV: {row}") from exc

    @staticmethod
    def _read_fasta_as_dict(path: str) -> Dict[str, str]:
        sequences_by_accession: Dict[str, str] = {}
        current_acc: Optional[str] = None
        current_seq_parts: List[str] = []

        with open(path, "r") as f:
            for raw_line in f:
                if not raw_line:
                    continue
                line = raw_line.rstrip("\n")
                if not line:
                    continue
                if line.startswith(">"):
                    # Flush previous
                    if current_acc is not None:
                        sequences_by_accession[current_acc] = "".join(current_seq_parts)
                    header_content = line[1:].strip()
                    # Accession is the first whitespace-delimited token
                    accession = header_content.split(None, 1)[0]
                    current_acc = accession
                    current_seq_parts = []
                else:
                    current_seq_parts.append(line.strip())
            # Flush final
            if current_acc is not None:
                sequences_by_accession[current_acc] = "".join(current_seq_parts)

        return sequences_by_accession

    @staticmethod
    def _extract_subsequence(full_sequence: Optional[str], start_1_based: int, end_1_based: int) -> Optional[str]:
        if full_sequence is None:
            return None
        if start_1_based <= 0 or end_1_based <= 0:
            return None
        # BLAST coordinates are 1-based inclusive; order can be reversed depending on alignment direction
        left = min(start_1_based, end_1_based)
        right = max(start_1_based, end_1_based)
        if left > len(full_sequence):
            return None
        # Slice is exclusive of end; adjust for 1-based inclusive
        return full_sequence[left - 1 : min(right, len(full_sequence))]

    @staticmethod
    def _reverse_complement(seq: str) -> str:
        return str(Seq(seq).reverse_complement())


def group_matches(results: Results, max_intron_length: int = 10_000, max_aa_overlap: int = 20) -> List[ProteinMatch]:
    """
    Group Match objects into ProteinMatch objects.
    - Groups by (query_accession, target_id)
    - Within a (query, target) group, further splits into clusters if adjacent
      matches on the target are separated by more than max_intron_length.
    """
    all_matches = results.matches()
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

            # Compute coverage booleans using query length (if available)
            query_len: Optional[int] = None
            if getattr(results, "_query_sequences_by_accession", None) is not None:
                seq_map: Dict[str, str] = results._query_sequences_by_accession or {}
                seq = seq_map.get(query_id)
                if seq:
                    query_len = len(seq)

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


