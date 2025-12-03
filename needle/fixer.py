import csv
import itertools
from needle.blast import Results, ProteinMatch
from needle.hits import hmmsearch
from Bio.Seq import Seq


def hmmsearch_matches(old_pm, full_target_len, hmm_file, sequences):
    idx_qstart = 999
    idx_qend = 999
    idx_tstart = 999
    idx_tend = 999

    expected_headers = "# target name        accession   tlen query name           accession   qlen   E-value  score"
    expected_eval_idx = 6
    expected_score_idx = 7
    assert expected_headers.replace("# target name", "target_name").replace("query name", "query_name").split()[expected_score_idx] == "score"
    assert expected_headers.replace("# target name", "target_name").replace("query name", "query_name").split()[expected_eval_idx] == "E-value"

    out = hmmsearch(hmm_file, sequences, score=False)
    matches = []
    has_headers = False

    for line in out.split("\n"):
        if line.startswith(expected_headers):
            has_headers = True
        if not line or line.startswith("#") or has_headers is False:
            continue

        parts = line.strip().split()
        frame = int(parts[0].replace("cand_",""))
        reverse = False if frame in (0, 1, 2) else True
        t_start = ((parts[idx_tstart]-1)*3+frame)+1
        t_end = parts[idx_tend]*3+frame
        if reverse is True:
            tt_end = full_target_len-(t_start-1)
            tt_start = full_target_len-t_end+1
            t_end = tt_end
            t_start = tt_start

        match = NucMatch(
            query_accession=old_pm.query_accession,
            target_accession=old_pm.target_id,
            query_start=parts[idx_qstart],
            query_end=parts[idx_qend],
            target_start=t_start,
            target_end=t_end,
            e_value=parts[idx_evalue],
            identity=None,
            on_reverse_strand=reverse
        )
        matches.append(match)

    return matches


def find_matches_at_locus(old_pm, old_matches, full_seq, start, end, hmm_file, step=5000):

    translations = all_frame_aa_sequences(full_seq, start, end)
    new_matches = hmmsearch_matches(old_pm, len(full_seq), hmm_file, translations)
    if not ProteinMatch.can_produce_single_sequence_from_matches(new_matches):
        return None

    old_query_start = min(m.query_start for m in old_matches)
    old_query_end = max(m.query_end for m in old_matches)
    old_nmatches = len(old_matches)

    new_query_start = min(m.query_start for m in nuc_matches)
    new_query_end = max(m.query_end for m in nuc_matches)
    new_nmatches = len(nuc_matches)

    if old_query_start == new_query_start and \
       old_query_end == new_query_end and \
       old_nmatches == new_nmatches:
        return None

    if start > step or end+step <= len(full_seq):
        more_matches = find_matches_at_locus(old_pm, new_matches, full_seq, max(0, start-step), min(len(full_seq), end+step))
        return more_matches if more_matches else new_matches

    return new_matches


def hmmsearch_at_locus(protein_match, fasta_file_path, hmm_file):
    sequence_by_acc = Results._read_fasta_as_dict(fasta_file_path)
    assert protein_match.target_id in sequence_by_acc
    target_full_sequence = sequence_by_acc[protein_match.target_id]

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
