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

Install Muscle Docker image

```
docker pull pegi3s/muscle
```

Install `HMMer` package. E.g. on MacOS run `brew install hmmer`.

Create Python virtualenv

```
python3 -v venv .venv
pip3 install -r requirements.txt
```

## Starting Point - Modules and Orthologs from KEGG

### Download a list of KO (KEGG Ortholog) numbers and names

```
curl https://rest.kegg.jp/list/ko -o data/ko.txt
echo "Ortholog ID\tOrtholog Name" | cat - data/ko.txt > data/ko.tsv
rm data/ko.txt
```

### Download a list of KEGG modules

```
curl https://rest.kegg.jp/list/module -o data/modules.txt
echo "Module ID\tModule Name" | cat - data/modules.txt > data/modules.tsv
rm data/modules.txt
```

### Fetch KO numbers for all the modules

The following creates `data/module_ko.tsv`

```
python3 scripts/fetch-kegg-module-ko.py
```

### Create list of consensus protein sequences for all the KO numbers

Download the HMM profiles from `https://www.genome.jp/ftp/db/kofam/`. The
`profiles.tar.gz` file is large, so this may take awhile.

Concatenate all the .hmm files together, e.g.

```
cat profiles/*.hmm > ko_full.hmm
```

Generate consensus protein sequence as a FASTA file, with

```
/opt/homebrew/Cellar/hmmer/3.4/bin/hmmemit -c ko_full.hmm > ko.fasta
```


## Workflow and Scripts

The general workflow looks like the following

  * Select a module to use
  * Generate "Blast Results" TSV (see header definitions in `needle/blast.py`) against interested genome accessions
  * Collate the blast results to "Ortholog Hits" database, consists of
    * A TSV of protein hit ID, ortholog ID, target genome accession, full protein hit stats
    * An accompanying TSV of protein hits broken down by fragments (likely exons)
      * Each row lists protein hit ID, fragment coordinate, contig coordinates and strand
    * FASTA files of detected protein sequences, split by Ortholog ID
  * Align detected protein sequences to create an MSA for each ortholog
    * Optionally add in protein sequences from SwissProt
  * Generate HTML+JS assets for visualizing the final data


Always activate the virtualenv first

```
source .venv/bin/activate
```

### Generating Query .faa for a KEGG Module

See below. The first argument is the module ID. The second argument is the FASTA file path. E.g.

```
PYTHONPATH=. python3 scripts/generate-module-fasta.py m00009
```

### Generate Ortholog Hits Database against a Genome Accession

```
./scripts/search-genome m00009 GCF_002042975.1
```

### Search in SwissProt for related proteins

```
scripts/mmseqs-swissprot-search proteins.faa results.tsv
```
