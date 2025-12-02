from needle.ortholog import load_modules
import argparse

parser = argparse.ArgumentParser(
    description="Generate FASTA file for KEGG module"
)
parser.add_argument("module_id", help="Module ID")
args = parser.parse_args()

modules = load_modules(
  'data/modules.tsv',
  'data/ko.tsv',
  'data/module_ko.tsv',
  'data/ko.fasta.gz',
  module_id=args.module_id.upper()
)

assert len(modules.keys()) == 1
module = modules[args.module_id.upper()]
module.create_consensus_aa_fasta(f"data/{args.module_id.lower()}_query.faa")
