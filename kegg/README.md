# Preparing Data from KEGG

This directory contains scripts to create the `data/kegg_module_proteins.tsv`
file.

Roughly, assuming you are in the top level directory for the repository

1. Download a list of KO (KEGG Ortholog) numbers and names

```
curl https://rest.kegg.jp/list/ko -o kegg/ko_list.txt
echo "KO ID\tKO Name" | cat - kegg/ko_list.txt > kegg/ko_list.tsv
rm kegg/ko_list.txt
```

2. Download a list of KEGG modules

```
curl https://rest.kegg.jp/list/module -o kegg/modules.txt
echo "Module ID\tModule Name" | cat - kegg/modules.txt > kegg/modules.tsv
rm kegg/modules.txt
```

3. Fetch KO numbers for all the modules

```
python3 kegg/fetch_modules.py
```

This last command creates `kegg/module_ko_list.txt`

4. Create list of consensus protein sequences for all the KO numbers

Download the HMM profiles from `https://www.genome.jp/ftp/db/kofam/`. The
`profiles.tar.gz` file is large, so this may take awhile.

Install `HMMer` package. E.g. on MacOS run `brew install hmmer`.

Generate consensus protein sequence as a FASTA file, with

```
/opt/homebrew/Cellar/hmmer/3.4/bin/hmmemit -c Kofam.hmm > kofam_consensus.fasta
```
