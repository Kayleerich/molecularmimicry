import argparse
from pathlib import Path
import pandas as pd
import numpy as np

parser = argparse.ArgumentParser(usage='%(prog)s [-h] -i INFILE (-f FILE1 FILE2 ... | -l LIST) [-o OUTID] [-s]', \
    description='Retrieve and filter highest-scoring HSPs from BLAST tabular output')

parser.add_argument('-i', '--hostfile', metavar='', action='store', type=str, required=True, help='Input BLAST tabular file for host species')
parser.add_argument('-o', '--outid', metavar='', action='store', type=str, required=False, help='Output files identifier')
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument('-f', '--files', metavar='', action='store', type=str, nargs='+', help='List of BLAST tabular files for control species')
group.add_argument('-l', '--list', metavar='', action='store', type=str, help='File containing list of BLAST tabular files for control species')
parser.add_argument('-s', '--silent', action='store_true', help='Silence terminal output')

args = parser.parse_args()

hostblastfile = Path(args.hostfile)

if args.outid == None:
    outid = ""
    compoutid = "-comparison"
    hostoutid = "-host_only"
else:
    outid = f"-{args.outid}"
    compoutid = f"{outid}-comparison"
    hostoutid = f"{outid}-host_only"

if args.list == None: 
    controlblastfiles = args.files
else:
    controlblastfiles = open(args.list).read().split('\n')
    controlblastfiles = [x for x in controlblastfiles if x]

if args.silent == False:
    print(f"Host species file: {hostblastfile}")
    print(f"Control species file(s): {', '.join(controlblastfiles[0:])}\n")


#get HSPs from control blastp files
controlHSP_dictionary = {}

for name in controlblastfiles:
    file = Path(name)
    with open(file, 'r') as controlblast:
        for cln in controlblast:
            cln = cln.split()
            cquery = cln[0]
            csubject0 = [cln[1], float(cln[2]), float(cln[11])]
            cpercid0 = float(cln[2])
            cscore0 = float(cln[11])
            if cquery not in controlHSP_dictionary.keys():
                csubject, cpercid, cscoremax = csubject0, cpercid0, cscore0
                controlHSP_dictionary[cquery] = list(csubject)
            else:
                centry0 = controlHSP_dictionary[cquery]
                chit, cpercid, cscoremax = centry0
                if cscore0 > cscoremax:
                    csubject, cpercid, cscoremax = csubject0, cpercid0, cscore0
                    controlHSP_dictionary[cquery] = list(csubject)
                elif cscore0 == cscoremax:
                    if cpercid0 > cpercid:
                        csubject, cpercid, cscoremax = csubject0, cpercid0, cscore0
                        controlHSP_dictionary[cquery] = list(csubject)
    controlblast.close()

#get HSPs from host blastp file
hostHSP_dictionary = {}
fullblastln = []

with open (hostblastfile, 'r') as hostblast:
    for ln in hostblast:
        hln = ln.split()
        hquery = hln[0]
        hsubject0 = [hln[1], float(hln[2]), float(hln[11])]
        hpercid0 = float(hln[2])
        hscore0 = float(hln[11])
        if hquery not in hostHSP_dictionary.keys():
            hsubject, hpercid, hscoremax = hsubject0, hpercid0, hscore0
            hostHSP_dictionary[hquery] = list(hsubject)
        else:
            hentry0 = hostHSP_dictionary[hquery]
            hhit, hpercid, hscoremax = hentry0
            if hscore0 > hscoremax:
                hsubject, hpercid, hscoremax = hsubject0, hpercid0, hscore0
                hostHSP_dictionary[hquery] = list(hsubject)
            elif hscore0 == hscoremax:
                if hpercid0 > hpercid:
                    hsubject, hpercid, hscoremax = hsubject0, hpercid0, hscore0
                    hostHSP_dictionary[hquery] = list(hsubject)
hostblast.close()

#compare HSPs from controls and host blastp results
hspcomparisons = []
hsphostonly = []

with open(f"topHSPs{compoutid}.tsv", 'w') as comparefile, open(f"topHSPs{hostoutid}.tsv", 'w') as hostonlyfile:
    for query in hostHSP_dictionary:
        hostHSP = hostHSP_dictionary.get(query)
        if query in controlHSP_dictionary.keys():
            controlHSP = controlHSP_dictionary.get(query)
            hspcomparisons.append([query, hostHSP[0], hostHSP[1], hostHSP[2], controlHSP[0], controlHSP[1], controlHSP[2]])
            comparefile.write(f"{query}\t{hostHSP[0]}\t{hostHSP[1]}\t{hostHSP[2]}\t{controlHSP[0]}\t{controlHSP[1]}\t{controlHSP[2]}\n")
        else:
            hsphostonly.append([query, hostHSP[0], hostHSP[1], hostHSP[2]])
            hostonlyfile.write(f"{query}\t{hostHSP[0]}\t{hostHSP[1]}\t{hostHSP[2]}\n")
comparefile.close()
hostonlyfile.close()


#filter k-mer HSPs
comparisonfile = pd.DataFrame(hspcomparisons)
comparisonfile.columns =  ["queryseqid", "hostseqid", "hostpident", "hostbitscore", "controlseqid", "controlpident", "controlbitscore"]

#calculate differences between host and control bitscore and %ID
def diff(a, b):
	return a - b
comparisonfile["bitscorediff"] = np.vectorize(diff)(comparisonfile["hostbitscore"], comparisonfile["controlbitscore"])
comparisonfile["iddiff"] = np.vectorize(diff)(comparisonfile["hostpident"], comparisonfile["controlpident"])

#calculate mean and standard deviation of differences
bitdiffmean = comparisonfile["bitscorediff"].mean()
iddiffmean = comparisonfile["iddiff"].mean()
bitstd = comparisonfile["bitscorediff"].std()
idstd = comparisonfile["iddiff"].std()

#get minimum of host bitscore and %ID
min_bit_diff = comparisonfile["bitscorediff"].mean() + (2 * comparisonfile["bitscorediff"].std())
min_id_diff = comparisonfile["iddiff"].mean() + comparisonfile["iddiff"].std()

#filter k-mer hits by minimum difference in bitscore and %ID
comparison_selected = comparisonfile[(comparisonfile["bitscorediff"] >= min_bit_diff) & (comparisonfile["iddiff"] >= min_id_diff)]

#get minimum bitscore and %ID
min_bit = comparison_selected["hostbitscore"].min()
min_id = comparison_selected["hostpident"].min()

#read in host-only HSP file
hostfile = pd.DataFrame(hsphostonly)
hostfile.columns = ["queryseqid", "hostseqid", "hostpident", "hostbitscore"]

#filter k-mer hits by minimum host bitscore and %ID
all_selected = pd.concat([hostfile[(hostfile["hostbitscore"] >= min_bit) & (hostfile["hostpident"] >= min_id)], comparison_selected], join='inner', ignore_index=True)

#write lists of selected k-mers and host proteins to file
allselectedlist = all_selected.values.tolist()
selected_query = set(all_selected['queryseqid'].unique())
selected_host = set(all_selected['hostseqid'].unique())

with open(f"selected_kmers_list{outid}.txt", 'w') as kmersfile:
    for ln in selected_query:
        kmersfile.write(f"{ln}\n")
kmersfile.close()

with open(f"host_proteins_list{outid}.txt", 'w') as hostlist:
    for ln in selected_host:
        hostlist.write(f"{ln}\n")
hostlist.close()

#get full blast lines for each query
with open(hostblastfile, 'r') as blastfile, open(f"blast_results-selected_HSPs{outid}.tsv", 'w') as kmerblast:
    for ln in blastfile:
        ln = ln.split()
        for selectedln in allselectedlist:
            if ln[0] == selectedln[0]:
                if ln[1] == selectedln[1]:
                    if float(ln[2]) == selectedln[2]:
                        if float(ln[11]) == selectedln[3]:
                            full = '\t'.join(ln)
                            kmerblast.write(f"{full}\n")
kmerblast.close()
blastfile.close()

#get rest of relevant stats
comparison_count = comparisonfile["hostseqid"].count()
comparisonselect_count = comparison_selected["hostseqid"].count()
hostselect_count = len(hostfile[(hostfile["hostbitscore"] >= min_bit) & (hostfile["hostpident"] >= min_id)])
host_count = hostfile["hostseqid"].count()
kmer_count = all_selected["hostseqid"].count()

#print relevant stats in separate file
with open(f"stats{outid}.txt", 'w') as statsfile:
	statsfile.write(f"Number of starting comparison HSPs\t{comparison_count}\n"\
		f"Number of starting host-only HSPs\t{host_count}\n\n"\

		f"--Number of k-mers selected--\n"\
		f"K-mers selected from comparison HSPs\t{comparisonselect_count}\n"\
        f"K-mers selected from host-only HSPs\t{hostselect_count}\n"\
        f"Total number of selected k-mers\t{kmer_count}\n\n"\

		f"--Bitscore--\n"\
		f"Mean\t{round(bitdiffmean, 4)}\n"\
		f"Std\t{round(bitstd, 4)}\n"\
		f"Minimum bitscore difference\t{round(min_bit_diff, 4)}\n"\
		f"Minimum bitscore\t{min_bit}\n\n"\

		f"--%ID--\n"\
		f"Mean\t{round(iddiffmean, 4)}\n"\
		f"Std\t{round(idstd, 4)}\n"\
		f"Minimum %ID difference\t{round(min_id_diff, 4)}\n"\
		f"Minimum %ID\t{min_id}")
statsfile.close()

if args.silent == False:
    print(f"Done. Files saved as:\n\
        topHSPs{compoutid}.tsv\n\
        topHSPs{hostoutid}.tsv\n\
        blast_results-selected_HSPs{outid}.tsv\n\
        selected_kmers_list{outid}.txt\n\
        host_proteins_list{outid}.txt\n\
        stats{outid}.txt")