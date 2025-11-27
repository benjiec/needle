from __future__ import annotations

import csv
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional
from Bio.Seq import Seq


@dataclass
class NucMatch:
    query_accession: str
    target_accession: str
    query_start: int
    query_end: int
    target_start: int
    target_end: int
    e_value: float
    identity: float
    on_reverse_strand: bool = False
    matched_sequence: Optional[str] = None

    # note: target_sequence and query_sequence match in direction
    query_sequence: Optional[str] = None
    target_sequence: Optional[str] = None
    target_sequence_upstream: Optional[str] = None
    target_sequence_downstream: Optional[str] = None

    def target_sequence_translated(self) -> str:
        if not self.target_sequence:
            return ""
        dna = self.target_sequence.upper()
        usable_len = (len(dna) // 3) * 3
        if usable_len == 0:
            return ""
        dna = dna[:usable_len]
        return str(Seq(dna).translate(table="Standard", to_stop=False))


@dataclass
class ProteinMatch:
    target_id: str
    matches: List[NucMatch]
    query_start: int
    query_end: int
    target_start: int
    target_end: int
    covers_start_to_end: bool
    likely_complete: bool
    query_overlap: bool

    def target_protein_sequence(self, results: "Results") -> str:
        # Determine query length if possible
        query_acc = self.matches[0].query_accession if self.matches else None
        query_len = None
        if query_acc and getattr(results, "_query_sequences_by_accession", None) is not None:
            seq_map: Dict[str, str] = results._query_sequences_by_accession or {}
            seq = seq_map.get(query_acc)
            if seq:
                query_len = len(seq)
        if query_len is None:
            # Fallback to range spanned by matches
            query_len = max((m.query_end for m in self.matches), default=0)
        # Collect candidates per query position (1-based indexing)
        choices: List[List[str]] = [[] for _ in range(query_len)]
        for m in self.matches:
            aa = m.target_sequence_translated()
            expected = m.query_end - m.query_start + 1
            limit = min(len(aa), expected)
            for i in range(limit):
                pos = m.query_start + i
                if 1 <= pos <= query_len:
                    choices[pos - 1].append(aa[i])
        # Render with -, single letter, or {a/b/c}
        out: List[str] = []
        for cands in choices:
            if not cands:
                out.append("-")
            else:
                uniq = sorted(set(cands))
                out.append(uniq[0] if len(uniq) == 1 else "{" + "/".join(uniq) + "}")
        return "".join(out)

    def pprint_target_protein_sequence(self) -> str:
        lines: List[str] = []
        # Sort along inferred 5'->3' orientation:
        # Forward strand (target_start < target_end): ascending by target_start
        # Reverse strand (target_start > target_end): descending by target_start
        reverse_sort = self.target_start > self.target_end
        for m in sorted(self.matches, key=lambda x: x.target_start, reverse=reverse_sort):
            indent = " " * (m.query_start - 1)
            lines.append(f"{indent}{m.target_sequence_translated()}")
        return "\n".join(lines)


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
        target_sequence_flanking_window: int = 10,
    ) -> None:
        self._results_tsv_path = results_tsv_path
        self._query_fasta_path = query_fasta_path
        self._target_fasta_path = target_fasta_path
        self._matches: Optional[List[NucMatch]] = None
        self._flank_window: int = max(0, int(target_sequence_flanking_window))

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
                        # Flanking sequences in same orientation as target_sequence
                        upstream, downstream = self._compute_flanks(
                            self._target_sequences_by_accession.get(match.target_accession, None),
                            match.target_start,
                            match.target_end,
                            match.on_reverse_strand,
                            self._flank_window,
                        )
                        match.target_sequence_upstream = upstream
                        match.target_sequence_downstream = downstream
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
            reverse = sstart > send
            t_start = min(sstart, send)
            t_end = max(sstart, send)

            match = NucMatch(
                query_accession=qacc,
                target_accession=sacc,
                query_start=qstart,
                query_end=qend,
                target_start=t_start,
                target_end=t_end,
                e_value=float(evalue_str.replace(",", "")),
                identity=float(pident_str.replace(",", "")),
                on_reverse_strand=reverse,
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

    @staticmethod
    def _compute_flanks(
        full_sequence: Optional[str],
        start_1_based: int,
        end_1_based: int,
        is_reverse: bool,
        flank_window: int,
    ) -> (Optional[str], Optional[str]):
        if full_sequence is None or flank_window <= 0:
            return None, None
        n = len(full_sequence)
        left = min(start_1_based, end_1_based)
        right = max(start_1_based, end_1_based)
        # Forward orientation in genome coordinates
        up_start = max(1, left - flank_window)
        up_end = max(0, left - 1)
        down_start = min(n + 1, right + 1)
        down_end = min(n, right + flank_window)

        upstream_genomic = full_sequence[up_start - 1 : up_end] if up_end >= up_start else ""
        downstream_genomic = full_sequence[down_start - 1 : down_end] if down_end >= down_start else ""

        if is_reverse:
            # Orientation matches target_sequence (reverse-complemented)
            upstream_oriented = str(Seq(downstream_genomic).reverse_complement()) if downstream_genomic else ""
            downstream_oriented = str(Seq(upstream_genomic).reverse_complement()) if upstream_genomic else ""
            return upstream_oriented, downstream_oriented

        return upstream_genomic, downstream_genomic


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
            cluster_by_query = sorted(cluster, key=lambda m: (m.query_start, m.query_end))

            # Compute query overlap
            query_overlap = False
            prev_q_start = None
            prev_q_end = None
            for m in cluster_by_query:
                qs = m.query_start
                qe = m.query_end
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
                    target_start=pm_t_start,
                    target_end=pm_t_end,
                    covers_start_to_end=covers_start_to_end,
                    likely_complete=likely_complete,
                    query_overlap=query_overlap,
                )
            )

    return protein_matches


