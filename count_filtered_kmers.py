import argparse
from pathlib import Path


parser = argparse.ArgumentParser(usage='%(prog)s [-h] -i INFILE', description='Count number of HSPs for each identity')
parser.add_argument('-i', '--infile', metavar='', action='store', type=str, required=True, help='Input BLAST tabular file for parasite v. host species')
args = parser.parse_args()

blastfile = args.infile

def count_filtered_kmers(blastfile):
    """Reports number identities in alignments and the number of HSPs for each"""
    query_id_dict = {}
    with open(blastfile, 'r') as filteredblast:
        for ln in filteredblast:
            ln = ln.split()
            length = int(ln[3])
            mism = int(ln[4])
            id = (length - mism)
            query_id_dict.setdefault(ln[0], []).append(id)

    id_dict = {}
    for query, num in query_id_dict.items():
        maxid = max(query_id_dict.get(query))
        id_dict.setdefault(maxid, []).append(query)

    return id_dict

id_dict = count_filtered_kmers(blastfile)
for id, num in id_dict.items():
    print(f"{id}\t{len(num)}")
##    Uncomment for full lists of kmers:
    # print(f"{id}\t{len(num)}\t{','.join(num)}")