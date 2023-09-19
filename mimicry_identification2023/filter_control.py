import argparse
from pathlib import Path

parser = argparse.ArgumentParser(usage='%(prog)s [-h] -i INFILE', description='Retrieve each query that surpasses threshhold from BLAST tabular output')
parser.add_argument('-i', '--infile', metavar='', action='store', type=str, required=True, help='Input BLAST tabular file for parasite v. control species')
args = parser.parse_args()

controlblast = args.infile

control_id_dict = {14: 12, 13: 10, 12: 9, 11: 9, 10: 8, 9: 8, 8: 7, 7: 6, 6: 6}

def control_blast_filter(file):
    controlblast = Path(file)
    controlHSPs_list = []
    with open(file, 'r') as controlblast:
        for ln in controlblast:
            ln = ln.split()
            query = ln[0]
            length = int(ln[3])
            if length in control_id_dict:
                minaa = control_id_dict.get(length)
                mism = int(ln[4])
                if  (length - mism) >= minaa:
                    controlHSPs_list.append(query)
    controlHSPs = sorted(set(controlHSPs_list))

    return controlHSPs

controlHSPs = control_blast_filter(controlblast)

for HSP in controlHSPs:
    print(HSP)
