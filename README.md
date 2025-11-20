# Needle

Tools to search for protein sequences in a genome, brute force. Designed to
work with novel genomes with cryptic gene structures, unconventional intron
models, etc.

## Setup

Install NCBI Docker image

`docker pull ncbi/blast`


## Scripts

Find protein sequence in a genome

`python3 scripts/blast-genome-search.py query.faa GCF_002042975.1 results.tsv --min-word-size 2`


