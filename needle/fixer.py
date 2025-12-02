import csv
import itertools
from needle.blast import Results
from needle.hits import hmmsearch
from Bio.Seq import Seq


def align_against_locus(match_tsv, fasta_file_path, hmm_file):

    data = []
    with open(match_tsv, 'r', newline='') as tsv_fh:
        reader = csv.DictReader(tsv_fh, delimiter='\t')
        for row in reader:
            data.append(row)
    
    sequence_by_acc = Results._read_fasta_as_dict(fasta_file_path)

    hits = []
    for hit in data:
        if hit["target_accession"] not in sequence_by_acc:
            continue
        if hit["protein_hit_id"].split("_")[0] not in hmm_file:
            continue
        hits.append(hit)

    kf = lambda h: h["protein_hit_id"]
    hits = sorted(hits, key=kf)

    for k, grouped_hits in itertools.groupby(hits, key=kf):
        protein_hit_id = k
        grouped_hits = list(grouped_hits)
        target_id = [h["target_accession"] for h in grouped_hits][0] 
        target_start = min([int(h["target_start"]) for h in grouped_hits])
        target_end = max([int(h["target_end"]) for h in grouped_hits])
        target_sequence = Results._extract_subsequence(sequence_by_acc[target_id], target_start, target_end)

        dna_seq = Seq(target_sequence)
        translations = []
        for frame in range(3):
            aa = dna_seq[frame:].translate(to_stop=False) # Translate entire sequence, including stops
            translations.append(str(aa))
        rev_comp_seq = dna_seq.reverse_complement()
        for frame in range(3):
            aa = rev_comp_seq[frame:].translate(to_stop=False)
            translations.append(str(aa))

        for h in grouped_hits:
            print(h["protein_hit_id"], h["target_accession"], h["target_start"], h["target_end"], h["query_start"], h["query_end"],
                  h["target_sequence"],
                  Seq(Results._extract_subsequence(sequence_by_acc[target_id], int(h["target_start"]), int(h["target_end"]))).translate(to_stop=False))

        out = hmmsearch(hmm_file, translations, score=False)
        print(out)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("match_tsv")
    parser.add_argument("fasta_fn")
    parser.add_argument("hmm_fn")
    args = parser.parse_args()

    align_against_locus(args.match_tsv, args.fasta_fn, args.hmm_fn) 
