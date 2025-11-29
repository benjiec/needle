# Needle

Tools to search and curate protein sequences in a genome, brute force. Designed
to work with novel genomes with cryptic gene structures, unconventional intron
models, etc.


## Setup

Install NCBI Docker image

```
docker pull ncbi/blast
```

Install MMSeqs2 Docker image

```
docker pull ghcr.io/soedinglab/mmseqs2
```

Setup SwissProt DB for MMSeqs2

```
scripts/mmseqs-swissprot-setup
```

Install `HMMer` package. E.g. on MacOS run `brew install hmmer`.

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

### Generating Amino Acid FASTA Files for KEGG Modules

See below. The first argument is the module ID. The second argument is the FASTA file path. E.g.

```
PYTHONPATH=. python3 scripts/generate-module-fasta.py M00009 m00009.faa
```

### Generating Blast Search Results

Find protein sequence in a genome, using the above example of TCP cycle module

```
python3 needle/blast-genome.py m00009.faa GCF_002042975.1 results.tsv --min-word-size 2
PYTHONPATH=. python3 scripts/show-protein-results.py m00009.faa ncbi-downloads/ncbi_dataset/data/GCF_002042975.1/GCF_002042975.1_ofav_dov_v1_genomic.fna results.tsv kegg-downloads/profiles
```

### Generating Ortholog Hits Database

The following command creates output files with prefix "m00009-hits"

```
PYTHONPATH=. python3 scripts/export-protein-results.py m00009.faa ncbi-downloads/ncbi_dataset/data/GCF_002042975.1/GCF_002042975.1_ofav_dov_v1_genomic.fna results.tsv kegg-downloads/profiles GCF_002042975.1 m00009-hits
```

Search in SwissProt for related proteins

```
scripts/mmseqs-swissprot-search m00009-hits_proteins.faa m00009-hits_swissprot_search.tsv
```
