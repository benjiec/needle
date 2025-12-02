import argparse

parser = argparse.ArgumentParser()
parser.add_argument("fn")
parser.add_argument("query")
parser.add_argument("--repl", default="")

args = parser.parse_args()

with open(args.fn, "r") as f:
  lines = f.readlines()
  for line in lines:
    line = line.strip()
    if line.startswith(">"):
      print(line)
    else:
      print(line.replace(args.query,args.repl))
