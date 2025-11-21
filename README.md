# Needle

Tools to search for protein sequences in a genome, brute force. Designed to
work with novel genomes with cryptic gene structures, unconventional intron
models, etc.

## Setup

Install NCBI Docker image

```
docker pull ncbi/blast
```

Create Python virtualenv

```
python3 -v venv .venv
pip3 install -r requirements.txt
```


## Scripts

Always activate the virtualenv first

```
source .venv/bin/activate
```

Find protein sequence in a genome

```
python3 needle/blast-genome.py example/query_on_GCF_002042975.1.faa GCF_002042975.1 results.tsv --min-word-size 2
python3 needle/show-protein-results.py examples/query_on_GCF_002042975.1.faa ncbi-downloads/ncbi_dataset/data/GCF_002042975.1/GCF_002042975.1_ofav_dov_v1_genomic.fna results.tsv 
```
