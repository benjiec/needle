from __future__ import annotations

import csv
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional


@dataclass
class Match:
    query_accession: str
    target_accession: str
    query_start: int
    query_end: int
    target_start: int
    target_end: int
    e_value: float
    identity: float
    matched_sequence: Optional[str] = None
    query_sequence: Optional[str] = None
    target_sequence: Optional[str] = None


@dataclass
class ProteinMatch:
    target_id: str
    matches: List[Match]
    query_start: int
    query_end: int
    target_start: int
    target_end: int
    covers_start_to_end: bool
    likely_complete: bool
    query_overlap: bool


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
        target_fasta_path: Optional[str] = None,
    ) -> None:
        self._results_tsv_path = results_tsv_path
        self._query_fasta_path = query_fasta_path
        self._target_fasta_path = target_fasta_path
        self._matches: Optional[List[Match]] = None

        self._query_sequences_by_accession: Optional[Dict[str, str]] = None
        self._target_sequences_by_accession: Optional[Dict[str, str]] = None

    def matches(self) -> List[Match]:
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

            matches: List[Match] = []
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
                        match.target_sequence = self._extract_subsequence(
                            self._target_sequences_by_accession.get(match.target_accession, None),
                            match.target_start,
                            match.target_end,
                        )
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

    def _row_to_match(self, row: List[str], header_index: Dict[str, int]) -> Optional[Match]:
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

            match = Match(
                query_accession=qacc,
                target_accession=sacc,
                query_start=int(qstart_str),
                query_end=int(qend_str),
                target_start=int(sstart_str),
                target_end=int(send_str),
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


def group_matches(results: Results, max_intron_length: int = 10_000) -> List[ProteinMatch]:
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
    def target_interval(m: Match) -> (int, int):
        return (min(m.target_start, m.target_end), max(m.target_start, m.target_end))

    grouped: Dict[tuple, List[Match]] = {}
    for m in all_matches:
        key = (m.query_accession, m.target_accession)
        grouped.setdefault(key, []).append(m)

    protein_matches: List[ProteinMatch] = []

    for (query_id, target_id), group in grouped.items():
        # Sort by target interval start to create distance-based clusters
        group_sorted_by_target = sorted(group, key=lambda m: target_interval(m)[0])

        clusters: List[List[Match]] = []
        current_cluster: List[Match] = []
        current_end: Optional[int] = None

        for m in group_sorted_by_target:
            start_t, end_t = target_interval(m)
            if not current_cluster:
                current_cluster = [m]
                current_end = end_t
            else:
                # Distance from previous end to next start
                distance = start_t - (current_end or start_t)
                if distance > max_intron_length:
                    clusters.append(current_cluster)
                    current_cluster = [m]
                    current_end = end_t
                else:
                    current_cluster.append(m)
                    current_end = max(current_end or end_t, end_t)
        if current_cluster:
            clusters.append(current_cluster)

        # Build ProteinMatch objects for each cluster
        for cluster in clusters:
            # For boolean computations, sort by query_start
            cluster_by_query = sorted(cluster, key=lambda m: (min(m.query_start, m.query_end), max(m.query_start, m.query_end)))

            # Compute query overlap
            query_overlap = False
            prev_q_start = None
            prev_q_end = None
            for m in cluster_by_query:
                qs = min(m.query_start, m.query_end)
                qe = max(m.query_start, m.query_end)
                if prev_q_start is not None:
                    if qs <= (prev_q_end or qs):
                        query_overlap = True
                        break
                prev_q_start, prev_q_end = qs, qe

            # Compute coverage booleans using query length (if available)
            query_len: Optional[int] = None
            if getattr(results, "_query_sequences_by_accession", None) is not None:
                seq_map: Dict[str, str] = results._query_sequences_by_accession or {}
                seq = seq_map.get(query_id)
                if seq:
                    query_len = len(seq)

            # Aggregate min/max coordinates
            q_min = min(min(m.query_start, m.query_end) for m in cluster)
            q_max = max(max(m.query_start, m.query_end) for m in cluster)
            t_min = min(target_interval(m)[0] for m in cluster)
            t_max = max(target_interval(m)[1] for m in cluster)

            # covers_start_to_end: union of matches includes position 1 and query_len
            if query_len is None:
                covers_start_to_end = False
            else:
                # Create coverage of query ends only (start and end must be covered)
                covers_start = any(min(m.query_start, m.query_end) <= 1 <= max(m.query_start, m.query_end) for m in cluster)
                covers_end = any(min(m.query_start, m.query_end) <= query_len <= max(m.query_start, m.query_end) for m in cluster)
                covers_start_to_end = covers_start and covers_end

            # likely_complete: if covers_start_to_end and each gap between adjacent
            # query segments is <= 10 amino acids.
            likely_complete = False
            if covers_start_to_end:
                # Check gaps on query axis
                ok = True
                for i in range(len(cluster_by_query) - 1):
                    a = cluster_by_query[i]
                    b = cluster_by_query[i + 1]
                    a_end = max(a.query_start, a.query_end)
                    b_start = min(b.query_start, b.query_end)
                    gap = max(0, b_start - a_end - 1)
                    if gap > 10:
                        ok = False
                        break
                likely_complete = ok

            protein_matches.append(
                ProteinMatch(
                    target_id=target_id,
                    matches=cluster_by_query,
                    query_start=q_min,
                    query_end=q_max,
                    target_start=t_min,
                    target_end=t_max,
                    covers_start_to_end=covers_start_to_end,
                    likely_complete=likely_complete,
                    query_overlap=query_overlap,
                )
            )

    return protein_matches


