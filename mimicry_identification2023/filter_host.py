import argparse
from pathlib import Path

parser = argparse.ArgumentParser(usage='%(prog)s [-h] -i INFILE', description='Retrieve each query that surpasses threshhold from BLAST tabular output v. 2')
parser.add_argument('-i', '--infile', metavar='', action='store', type=str, required=True, help='Input BLAST tabular file for parasite v. host species')
args = parser.parse_args()

hostblastfile = args.infile

host_id_dict = {14: 12, 13: 12, 12: 11, 11: 10, 10: 10}

def host_blast_filter(file):
    hostblast = Path(file)
    hostHSPs_list = []
    with open(file, 'r') as hostblast:
        for ln in hostblast:
            ln = ln.split()
            length = int(ln[3])
            if length in host_id_dict:
                minaa = host_id_dict.get(length)
                mism = int(ln[4])
                if  (length - mism) >= minaa:
                    hostHSPs_list.append('\t'.join(ln))
    hostHSPs = sorted(set(hostHSPs_list))

    return hostHSPs

if Path(hostblastfile).is_file():
    hostHSPs = host_blast_filter(hostblastfile)
    for HSP in hostHSPs:
        print(HSP)
