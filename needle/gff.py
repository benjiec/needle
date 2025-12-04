from typing import Any, Dict, Iterable, List, Tuple
from BCBio import GFF
from .match import Match, ProteinHit


def _cds_transcript_order_key(strand: str):
    # 5'->3' transcript order: '+' ascending coordinates, '-' descending
    if strand == "-":
        return lambda feat: (-_get_start(feat), -_get_end(feat))
    return lambda feat: (_get_start(feat), _get_end(feat))


def _compute_aa_span(length_nt: int, phase: int, cumulative_aa: int) -> Tuple[int, int, int]:
    usable_nt = max(0, length_nt - max(0, phase))
    codons = usable_nt // 3
    if codons <= 0:
        return cumulative_aa + 1, cumulative_aa, cumulative_aa  # empty contribution
    aa_start = cumulative_aa + 1
    aa_end = aa_start + codons - 1
    return aa_start, aa_end, cumulative_aa + codons


def _safe_int(v: Any, default: int = 0) -> int:
    try:
        return int(v)
    except Exception:
        return default


def _get_attr_dict(feat: Any) -> Dict[str, Any]:
    # BCBio.GFF provides Bio.SeqFeature with .qualifiers: Dict[str, List[str]]
    q = getattr(feat, "qualifiers", None)
    if isinstance(q, dict):
        return q
    return {}


def _get_type(feat: Any) -> str:
    # BCBio SeqFeature has .type
    v = getattr(feat, "type", None)
    return v if isinstance(v, str) else ""


def _get_start(feat: Any) -> int:
    # BCBio SeqFeature location is 0-based half-open; convert to 1-based inclusive
    loc = getattr(feat, "location", None)
    if loc is not None:
        try:
            return int(loc.start) + 1
        except Exception:
            pass
    return _safe_int(getattr(feat, "start", None))


def _get_end(feat: Any) -> int:
    loc = getattr(feat, "location", None)
    if loc is not None:
        try:
            return int(loc.end)
        except Exception:
            pass
    return _safe_int(getattr(feat, "end", None))


def _get_strand(feat: Any) -> str:
    # BCBio uses strand in location: +1 / -1 / None
    loc = getattr(feat, "location", None)
    s = None
    if loc is not None:
        s = getattr(loc, "strand", None)
    if s == -1:
        return "-"
    return "+"


def _get_frame(feat: Any) -> int:
    # BCBio stores phase as qualifier 'phase' or 'codon_start' (1..3)
    q = _get_attr_dict(feat)
    phase_vals = q.get("phase") or q.get("Phase") or []
    if phase_vals:
        try:
            return int(phase_vals[0])
        except Exception:
            return 0
    codon_start = q.get("codon_start") or q.get("Codon_start") or []
    if codon_start:
        try:
            cs = int(codon_start[0])
            return {1: 0, 2: 1, 3: 2}.get(cs, 0)
        except Exception:
            return 0
    return 0


def _get_seqid_from_record_or_feat(record: Any, feat: Any) -> str:
    # SeqRecord.id is the seqname/contig
    rec_id = getattr(record, "id", None)
    if isinstance(rec_id, str) and rec_id:
        return rec_id
    return ""


def _flatten_features_from_bcbb_records(records: Iterable[Any]) -> Iterable[Tuple[Any, Any]]:
    """
    Given records from bcbb.gff parser, yield (record, feature) pairs.
    Traverse SeqRecord.features recursively to collect CDS features.
    """
    def _iter_feats(rec, feats):
        for f in feats or []:
            yield rec, f
            sub = getattr(f, "sub_features", None) or getattr(f, "sub_features", None)
            if sub:
                yield from _iter_feats(rec, sub)
    for rec in records:
        feats = getattr(rec, "features", None) or []
        yield from _iter_feats(rec, feats)


def _parse_with_bcbb(gff_path: str) -> Iterable[Tuple[Any, Any]]:
    # Use BCBio.GFF (bcbb) parser exclusively
    parser = GFF.GFFParser()
    with open(gff_path, "r") as in_handle:
        records = parser.parse(in_handle)
        return list(_flatten_features_from_bcbb_records(records))

def parse_gff_to_hits(gff_path: str, protein_id_attr: str = "protein_id") -> List[ProteinHit]:
    """
    Parse a GFF (GFF3 preferred) file and convert CDS features into Match objects,
    grouped by the same attribute 'protein_id' into ProteinHit objects.
    Uses chapmanb/bcbb (BCBio.GFF) parser exclusively.
    """
    records_and_feats: Iterable[Tuple[Any, Any]] = _parse_with_bcbb(gff_path)

    def _group(records_feats: Iterable[Tuple[Any, Any]]) -> Tuple[Dict[str, List[Tuple[Any, Any]]], Dict[str, str], Dict[str, str]]:
        cds_by_protein: Dict[str, List[Tuple[Any, Any]]] = {}
        seqid_by_protein: Dict[str, str] = {}
        strand_by_protein: Dict[str, str] = {}
        for rec, feat in records_feats:
            if _get_type(feat) != "CDS":
                continue
            attrs = _get_attr_dict(feat)
            values = attrs.get(protein_id_attr, [])
            protein_id = None
            if isinstance(values, list):
                protein_id = values[0] if values else None
            elif isinstance(values, str):
                protein_id = values
            if not protein_id:
                continue
            cds_by_protein.setdefault(protein_id, []).append((rec, feat))

            seqid = _get_seqid_from_record_or_feat(rec, feat)
            strand = _get_strand(feat)
            prev_seqid = seqid_by_protein.get(protein_id)
            prev_strand = strand_by_protein.get(protein_id)
            if prev_seqid is None:
                seqid_by_protein[protein_id] = seqid
            elif prev_seqid != seqid:
                raise ValueError(f"Inconsistent seqid for protein_id '{protein_id}': {prev_seqid} vs {seqid}")
            if prev_strand is None:
                strand_by_protein[protein_id] = strand
            elif prev_strand != strand:
                raise ValueError(f"Inconsistent strand for protein_id '{protein_id}': {prev_strand} vs {strand}")
        return cds_by_protein, seqid_by_protein, strand_by_protein

    cds_by_protein, seqid_by_protein, strand_by_protein = _group(records_and_feats)

    protein_hits: List[ProteinHit] = []

    for protein_id, rf_list in cds_by_protein.items():
        strand = strand_by_protein[protein_id]
        seqid = seqid_by_protein[protein_id]

        feats_only = [f for _r, f in rf_list]
        feats_sorted = sorted(feats_only, key=_cds_transcript_order_key(strand))

        matches: List[Match] = []
        cumulative_aa = 0
        min_coord = None
        max_coord = None

        for f in feats_sorted:
            start = _get_start(f)
            end = _get_end(f)
            if not start or not end:
                continue
            length_nt = abs(end - start) + 1
            phase = _get_frame(f)

            aa_start, aa_end, cumulative_aa = _compute_aa_span(length_nt, phase, cumulative_aa)
            if aa_end < aa_start:
                continue

            if strand == "+":
                t_start, t_end = start, end
            else:
                t_start, t_end = end, start  # 5'->3' of gene: start > end on reverse

            matches.append(
                Match(
                    query_accession=protein_id,
                    target_accession=seqid,
                    query_start=aa_start,
                    query_end=aa_end,
                    target_start=t_start,
                    target_end=t_end,
                    e_value=0.0,
                    identity=0.0,
                )
            )

            if min_coord is None or min(start, end) < min_coord:
                min_coord = min(start, end)
            if max_coord is None or max(start, end) > max_coord:
                max_coord = max(start, end)

        if not matches:
            continue

        if strand == "+":
            prot_t_start, prot_t_end = min_coord, max_coord
        else:
            prot_t_start, prot_t_end = max_coord, min_coord

        protein_hits.append(
            ProteinHit(
                matches=matches,
                query_start=1,
                query_end=cumulative_aa,
                target_start=prot_t_start,
                target_end=prot_t_end,
                _protein_hit_id=protein_id,
            )
        )

    return protein_hits


