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


## Workflow and Scripts

The general workflow looks like the following

  * Start with "Ortholog Reference" composed of 3 TSVs and a FASTA
    * A TSV listing Module ID and Module Name
    * A TSV listing Ortholog ID and Ortholog Name
    * A TSV listing Module ID and Module Definition
      * Module definition is a comma delimited list of either Ortholog ID or Module ID
    * A FASTA of Ortholog consensus amino acid sequences

  * Generate "Blast Results" TSV (see header definition in `needle/blast.py`)

  * Generate "Ortholog Hits" database composed of the following files
    * A FASTA file of protein hits
    * A TSV of protein hit ID, ortholog ID, target genome accession, full protein hit stats
    * An accompanying TSV of protein hits broken down by fragments (likely exons)
      * Each row lists protein hit ID, fragment coordinate, contig coordinates and strand

  * Generate "Ortholog MSA" database composed of the following files
    * For each Ortholog, MSAs in Stockholm Format

  * Generate HTML+JS assets for visualizing the final data


Always activate the virtualenv first

```
source .venv/bin/activate
```

### Generating Blast Search Results

Find protein sequence in a genome

```
python3 needle/blast-genome.py example/query_on_GCF_002042975.1.faa GCF_002042975.1 results.tsv --min-word-size 2
python3 needle/show-protein-results.py examples/query_on_GCF_002042975.1.faa ncbi-downloads/ncbi_dataset/data/GCF_002042975.1/GCF_002042975.1_ofav_dov_v1_genomic.fna results.tsv 
```
