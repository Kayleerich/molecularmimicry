import argparse
from pathlib import Path


filetype = ["blast-query", "blast-subject", "phobius", "segmasker"]
parser = argparse.ArgumentParser(usage='%(prog)s [-h] -i INFILE -t {'+','.join(filetype)+'}', description='Make BED file from given file')

parser.add_argument('-i', '--input', metavar='', action='store', type=str, required=True, help='Input file')
parser.add_argument('-t', '--filetype', metavar='', choices = filetype, type=str, required=True, help='Input file type, choices are: '+', '.join(filetype))

args = parser.parse_args()

file = Path(args.input)
intype = args.filetype

def make_bed(file, intype):
	"""Make bed file from BLAST tabular output, SEGmasker interval output or Phobius short output"""
	import re

	with open(f"{file.stem}.bed", 'w') as bedfile:
		if "blast" in intype:
			with open(file, "r") as blastfile:
				for ln in blastfile: 
					ln = ln.split()
					if intype == "blast-subject":
						#bed file for subject hits, columns 8/9 are subject start/end
						bed = f"{ln[1]}\t{int(ln[8]) - 1}\t{int(ln[9])}"
					else:
						#bed file for query hits, columns 6/7 are query start/end
						bed = f"{ln[0]}\t{int(ln[6]) - 1}\t{int(ln[7])}"
					#bed_list.append(bed)
					bedfile.write(f"{bed}\n")
				blastfile.close()

		elif intype == "phobius":
			with open(file, 'r') as phobiusout:
				for ln in phobiusout:
					ln = ln.split()
					if ln[2] == "Y":
						#find the signal peptide start and finish string
						start = re.search(r"n\d+-", ln[3]) 
						end = re.search(r"-\d+c", ln[3])
						#make signal readable 
						start = start.group()
						end = end.group()
						start = int(start[1:-1]) - 1
						end = int(end[1:-1]) - 1
						signal = f"{ln[0]}\t{start}\t{end}" #add gene name to first column 
						bedfile.write(f"{signal}\n")
				phobiusout.close()
		
		elif intype == "segmasker":
			with open(file, 'r') as segintervals:
				for ln in segintervals:
					if ln.startswith(">"): 
						name = ln.split("|")[1]
						name = name.split(" ")[0]
					else:
						start = ln.split(" ")[0]
						end = int(ln.split(" ")[2]) + 1
						bedfile.write(f"{name}\t{start}\t{end}\n")

		bedfile.close()


make_bed(file, intype)
