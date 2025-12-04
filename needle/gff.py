from typing import Dict, List, Tuple

from .match import Match, ProteinHit


def _cds_transcript_order_key(strand: str):
    # 5'->3' transcript order: '+' ascending coordinates, '-' descending
    if strand == "-":
        return lambda feat: (-feat.start, -feat.end)
    return lambda feat: (feat.start, feat.end)


def _compute_aa_span(length_nt: int, phase: int, cumulative_aa: int) -> Tuple[int, int, int]:
    usable_nt = max(0, length_nt - max(0, phase))
    codons = usable_nt // 3
    if codons <= 0:
        return cumulative_aa + 1, cumulative_aa, cumulative_aa  # empty contribution
    aa_start = cumulative_aa + 1
    aa_end = aa_start + codons - 1
    return aa_start, aa_end, cumulative_aa + codons


def parse_gff_to_hits(gff_path: str, protein_id_attr: str = "protein_id") -> List[ProteinHit]:
    """
    Parse a GFF (GFF3 preferred) file and convert CDS features into Match objects,
    grouped by the same attribute 'protein_id' into ProteinHit objects.
    - query_accession = protein_id
    - target_accession = seqid (contig/chromosome)
    - Match.target_start/target_end honor strandness:
      '+' strand => start <= end; '-' strand => start > end (5'->3' of gene)
    - ProteinHit.protein_hit_id is set explicitly from protein_id
    """
    try:
        import gffutils  # type: ignore
    except Exception as exc:
        raise ImportError(
            "gffutils is required to parse GFF. Install with 'pip install gffutils'."
        ) from exc

    db = gffutils.create_db(
        gff_path,
        dbfn=":memory:",
        force=True,
        keep_order=True,
        merge_strategy="merge",
        sort_attribute_values=True,
    )

    # Collect CDS grouped by protein_id
    cds_by_protein: Dict[str, List] = {}
    seqid_by_protein: Dict[str, str] = {}
    strand_by_protein: Dict[str, str] = {}

    for feat in db.features_of_type("CDS", order_by=("seqid", "start")):
        attrs = getattr(feat, "attributes", {}) or {}
        values = attrs.get(protein_id_attr, [])
        if not values:
            # Skip CDS without the desired protein_id attribute
            continue
        protein_id = values[0]
        cds_by_protein.setdefault(protein_id, []).append(feat)

        # Ensure consistent seqid/strand within a protein_id
        prev_seqid = seqid_by_protein.get(protein_id)
        prev_strand = strand_by_protein.get(protein_id)
        if prev_seqid is None:
            seqid_by_protein[protein_id] = feat.seqid
        elif prev_seqid != feat.seqid:
            raise ValueError(f"Inconsistent seqid for protein_id '{protein_id}': {prev_seqid} vs {feat.seqid}")
        if prev_strand is None:
            strand_by_protein[protein_id] = feat.strand
        elif prev_strand != feat.strand:
            raise ValueError(f"Inconsistent strand for protein_id '{protein_id}': {prev_strand} vs {feat.strand}")

    protein_hits: List[ProteinHit] = []

    for protein_id, feats in cds_by_protein.items():
        if not feats:
            continue
        strand = strand_by_protein[protein_id]
        seqid = seqid_by_protein[protein_id]
        feats_sorted = sorted(feats, key=_cds_transcript_order_key(strand))

        matches: List[Match] = []
        cumulative_aa = 0

        min_coord = None
        max_coord = None

        for f in feats_sorted:
            start = int(f.start)
            end = int(f.end)
            length_nt = abs(end - start) + 1
            # gffutils uses .frame for GFF phase; often provided as a string "0"/"1"/"2" or ".".
            # Convert to int and default to 0 if missing/invalid.
            try:
                phase = 0 if f.frame in (None, ".", "") else int(f.frame)
            except Exception:
                phase = 0

            aa_start, aa_end, cumulative_aa = _compute_aa_span(length_nt, phase, cumulative_aa)
            if aa_end < aa_start:
                # No usable codons from this fragment; skip
                continue

            if strand == "+":
                t_start, t_end = start, end
            else:
                t_start, t_end = end, start  # 5'->3' of gene: start > end on reverse

            m = Match(
                query_accession=protein_id,
                target_accession=seqid,
                query_start=aa_start,
                query_end=aa_end,
                target_start=t_start,
                target_end=t_end,
                e_value=0.0,
                identity=0.0,
            )
            matches.append(m)

            if min_coord is None or min(start, end) < min_coord:
                min_coord = min(start, end)
            if max_coord is None or max(start, end) > max_coord:
                max_coord = max(start, end)

        if not matches:
            continue

        # Protein-level coordinates
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


