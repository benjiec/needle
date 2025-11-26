from needle.ortholog import load_modules
import argparse

parser = argparse.ArgumentParser(
    description="Generate FASTA file for KEGG module"
)
parser.add_argument("module_id", help="Module ID")
parser.add_argument("output_fasta_filename", help="output FASTA filename")
args = parser.parse_args()

modules = load_modules(
  'kegg/modules.tsv',
  'kegg/ko_list.tsv',
  'kegg/module_ko_list.tsv',
  'kegg/ko_all_consensus.fasta.gz',
  module_id=args.module_id
)

assert len(modules.keys()) == 1
module = modules[args.module_id]
module.create_consensus_aa_fasta(args.output_fasta_filename)
