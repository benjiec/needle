# Needle

Tools to search and curate protein sequences in a genome, brute force, starting
with protein sequences and raw genomic DNA data. Designed to work with novel
genomes with cryptic gene structures, unconventional intron models, etc. For
example, for some coral and many symbiodinium species, training gene prediction
models suffers from insufficient gold-set data, and conducting RNAseq
experiments is diffuclt.


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

## Prepare Initial Data

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

### Define custom modules

Add custom modules, with KO numbers, using the following approach. Define the
modules first in the custom files, then add them to master lists.

```
cat custom_module_ko.tsv >> module_ko.tsv 
cat custom_modules.tsv >> modules.tsv  
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

Put all the HMM profiles under `kegg-downloads/profiles`. E.g.
`kegg-downloads/profiles/K00030.hmm` should exist. A couple of the scripts
below depends on this directory.

### Prepare List of Genomes

There is a list of Coral genomes in `data/genomes_coral.txt`, and a list of
Symbiodinium (algae) genomes in `data/genomes_algae.txt`. Run the following
command to generate `data/genomes.tsv`, which includes genome name and taxonomy
information.

```
PYTHONPATH=. python3 scripts/fetch-genomes.py data/genomes_coral.txt
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

The following script puts outputs in `data/m0009_results` directory

```
./scripts/search-genome m00009 GCF_002042975.1
```

Or if you have a list of genome accessions in a file, e.g. `genomes.txt`, then do

```
./scripts/search-genomes m00009 genomes.txt
```

### Generating MSAs

To generate MSAs and their visuals, run the following command for each module.
This command runs the four commands below this, for each KO number.

```
./scripts/generate-msas m00009
```

To generate MSAs, for a module and a KO, use the following commands. Each
script puts a MSA in FASTA format in `data/m00009_results/m00009-alignments`

```
./scripts/muscle-ko m00009 K00030
./scripts/hmmalign-ko m00009 K00030
```

Note: *s (STOP codons) from BLAST search are preserved through MUSCLE by first
replacing them with Zs, then run MUSCLE, then replace Zs from the output of
MUSCLE back to *s.

SVG files (which can be opened via Chrome and other browsers) visualizing the
MSAs can generated with the following commands. In this example, the SVG files
are in `data/m00009_results/m00009-msas`

```
./scripts/mk-msa-vis m00009 K00030 muscle
./scripts/mk-msa-vis m00009 K00030 hmmalign
```

### Web Summaries of Modules and Orthologs

Put the module IDs to be shown in `pages/modules.tsv`, then either start a
local server at the root of the repository working directory, e.g.

```
python3 -m http.server 8080
```

and visit `localhost:8080/pages/`. 

Or, commit any changes, merge into the "pages" branch, and push to github, and
then visit `https://benjiec.github.io/needle/pages/`.


### Compare / Sanity check Annotated Proteins against Needle Hits

Use the following script to compare, for a given HMM model, how NCBI annotated
proteins (i.e. in `protein.faa` and `genomic.gff`) compare against protein
found by Needle.

```
PYTHONPATH=. python3 scripts/compare-gff-with-match.py K00024 GCF_002042975.1 data/m00009_results/matches.tsv
```


### Search in SwissProt for related proteins

```
scripts/mmseqs-swissprot-search proteins.faa results.tsv
```
