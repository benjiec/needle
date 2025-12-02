import time
import argparse
from needle.ncbi import fetch_and_append_taxonomy

parser = argparse.ArgumentParser()
parser.add_argument("accessions_fn")
args = parser.parse_args()

with open(args.accessions_fn, "r") as f:
    lines = f.readlines()
    for line in lines:
        accession = line.strip()
        assert len(accession.split()) == 1
        print(accession)
        fetch_and_append_taxonomy(accession, "data/genomes.tsv")
        time.sleep(1)
