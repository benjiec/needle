import re
import os
import subprocess
import tempfile


def run_command(cmd: str):
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def parse_hmmsearch_domtbl(domtbl_path):
    expected_header = "# target name  accession  tlen  query name  accession  qlen  E-value  score  bias  #  of  c-Evalue  i-Evalue  score  bias  from  to  from  to"
    idx_target = 0
    idx_eval = 6
    idx_score = 7
    idx_q_from = 15
    idx_q_to = 16
    idx_t_from = 17
    idx_t_to = 18

    # first, sanity check these indices
    expected_header_parts = re.split(r'\s\s+', expected_header)
    assert expected_header_parts[idx_eval] == "E-value"
    assert expected_header_parts[idx_score] == "score"
    assert expected_header_parts[idx_target] == "# target name"
    assert expected_header_parts[idx_q_from] == "from"
    assert expected_header_parts[idx_q_to] == "to"
    assert expected_header_parts[idx_t_from] == "from"
    assert expected_header_parts[idx_t_to] == "to"

    has_headers = False
    expected_header = " ".join(expected_header.split())

    matches = []
    with open(domtbl_path, "r") as domf:
        for line in domf:
            if " ".join(line.split()).startswith(expected_header):
                has_headers = True
            if not line or line.startswith("#") or has_headers is False:
                continue
            parts = line.strip().split()
            match = dict(
                target_name = parts[idx_target],
                evalue = float(parts[idx_eval]),
                score = float(parts[idx_score]),
                hmm_from = int(parts[idx_q_from]),
                hmm_to = int(parts[idx_q_to]),
                target_from = int(parts[idx_t_from]),
                target_to = int(parts[idx_t_to])
            )
            matches.append(match)

    assert has_headers
    return matches


def hmmsearch(hmm_file_name, sequences):
    with tempfile.TemporaryDirectory() as tmpdir:
        fasta_path = os.path.join(tmpdir, "cands.faa")
        domtbl_path = os.path.join(tmpdir, "out.domtbl")
        with open(fasta_path, "w") as f:
            for i, cand in enumerate(sequences):
                f.write(f">cand_{i}\n{cand}\n")
        cmd = ["hmmsearch", "--domtblout", domtbl_path, hmm_file_name, fasta_path]
        run_command(cmd)
        return parse_hmmsearch_domtbl(domtbl_path)
