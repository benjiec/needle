# Preparing Data from KEGG

This directory contains the following files

  * `ko_list.tsv` - Ortholog ID and name
  * `modules.tsv` - Module ID and for each module, a comma-delimited list of KO IDs
  * `ko_full_consensus.fasta` - FASTA file with consensus sequence for each KO

To re-generate these files, do the following.

### Download a list of KO (KEGG Ortholog) numbers and names

```
curl https://rest.kegg.jp/list/ko -o kegg/ko_list.txt
echo "Ortholog ID\tOrtholog Name" | cat - kegg/ko_list.txt > kegg/ko_list.tsv
rm kegg/ko_list.txt
```

### Download a list of KEGG modules

```
curl https://rest.kegg.jp/list/module -o kegg/modules.txt
echo "Module ID\tModule Name" | cat - kegg/modules.txt > kegg/modules.tsv
rm kegg/modules.txt
```

### Fetch KO numbers for all the modules

```
python3 kegg/fetch_modules.py
```

This last command creates `kegg/module_ko_list.txt`

### Create list of consensus protein sequences for all the KO numbers

Install `HMMer` package. E.g. on MacOS run `brew install hmmer`.

Download the HMM profiles from `https://www.genome.jp/ftp/db/kofam/`. The
`profiles.tar.gz` file is large, so this may take awhile.

Concatenate all the .hmm files together, e.g.

```
cat profiles/*.hmm > ko_full.hmm
```

Generate consensus protein sequence as a FASTA file, with

```
/opt/homebrew/Cellar/hmmer/3.4/bin/hmmemit -c ko_full.hmm > ko_full_consensus.fasta
```
