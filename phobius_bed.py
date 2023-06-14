import argparse
import re
from pathlib import Path

parser = argparse.ArgumentParser(description='Make BED file for cleaved portion of proteins')
parser.add_argument('-i', '--phobius_file', action='store', type=str, required=True, help='Phobius short output file')
args = parser.parse_args()

file = Path(args.phobius_file)

with open(file, 'r') as phobiusout, open(f"{file.stem}.bed", "w") as bedfile:
    for ln in phobiusout:
        ln = ln.split()
        if ln[2] == "Y":
            start = 0
            end = re.search(r"-\d+c\d+/\d", ln[3])
            end = end.group()
            end = end.split('c')[1]
            end = end.split('/')[0]
            signal = f"{ln[0]}\t{start}\t{end}"
            bedfile.write(f"{signal}\n")
    phobiusout.close()
