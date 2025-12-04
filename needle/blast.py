import csv
from typing import Dict, List, Optional
from Bio.Seq import Seq

from .match import NucMatch, extract_subsequence, extract_subsequence_strand_sensitive


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
                        match.query_sequence = extract_subsequence(
                            self._query_sequences_by_accession.get(match.query_accession, None),
                            match.query_start,
                            match.query_end,
                        )
                    if self._target_sequences_by_accession is not None:
                        seq = extract_subsequence_strand_sensitive(
                            self._target_sequences_by_accession.get(match.target_accession, None),
                            match.target_start,
                            match.target_end,
                        )
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
