import pandas as pd
import numpy as np
import os
from os.path import exists as file_exists
import glob
import argparse
from pathlib import Path
import sys

parser = argparse.ArgumentParser(usage='%(prog)s [-h] -i BLAST_FILE -p PATHOGEN_SPECIES -s HOST_SPECIES [-f FILE_IDENTIFIER -d POPS_PATH -q QSASA_THRESHHOLD]', description='Get paired pathogen-host regions from BLAST file')

parser.add_argument('-i', '--input', metavar='', action='store', type=str, required=True, help='Input BLAST file')
parser.add_argument('-p', '--pathogen_species', metavar='', action='store', type=str, required=True, help='Pathogen species name/subdirectory identifier')
parser.add_argument('-s', '--host_species', metavar='', action='store', type=str, required=True, help='Host species name/subdirectory identifier')
parser.add_argument('-k', '--kmer_size', metavar='', action='store', type=int, required=True, help='K-mer size used')
parser.add_argument('-d', '--pops_path', metavar='', action='store', type=str, required=False, help='Path to POPS parent directory')
parser.add_argument('-f', '--fileid', metavar='', action='store', type=str, required=False, help='File identifier')
parser.add_argument('-q', '--min_qsasa', metavar='', action='store', type=int, required=False, help='Minimum average QSASA threshold, default 0.75')

args = parser.parse_args()

blastfilename = Path(args.input)
p_name =  args.pathogen_species
h_name = args.host_species

if args.min_qsasa:
    min_qsasa = args.min_qsasa
else:
    min_qsasa = 0.75

k=args.kmer_size
if args.pops_path:
    pops_path = Path(args.pops_path)
else:
    pops_path = Path(".")
hpops_path = Path(f"{pops_path}/{h_name}/pops")
if not os.path.isdir(hpops_path):
    hpops_path.mkdir(parents=True, exist_ok=True)
    
ppops_path = Path(f"{pops_path}/{p_name}/pops")
if not os.path.isdir(ppops_path):
    ppops_path.mkdir(parents=True, exist_ok=True)

if args.fileid:
    fileid = args.fileid
else:
    fileid = f"{p_name}.{h_name}"

ppops_files_all = glob.glob(f"{ppops_path}/pops_*.out")
if not ppops_files_all:
    sys.exit(f"Warning: no pops files found for pathogen in {ppops_path}")
ppops_files = [file for file in ppops_files_all if "-model_v2" not in file]

hpops_files_all = glob.glob(f"{hpops_path}/pops_*.out")
if not hpops_files_all:
    sys.exit(f"Warning: no pops files found for host in {hpops_path}")
hpops_files = [file for file in hpops_files_all if "-model_v2" not in file]

coords_list_dicts = {}

## create dictionary of all ranges
with open(blastfilename, 'r') as blast_file:
    for ln in blast_file:
        ln = ln.split()
        pair = ln[0].split('/')[0] + "." + ln[1] ## might need to be changed depending on sequence names: ln[0].split('|')[1] + "." + ln[1]

        kcoords = ln[0].split('/')[1]
        kstart = int(kcoords.split(':')[0]) - 1
        kend = int(kcoords.split(':')[1])

        qstart = (kstart + int(ln[6]))
        qend = kend - (k - int(ln[7])) + 1
        query_range = [qstart, qend]

        try:
            coords_list_dicts[pair]['query'].append(query_range)
        except KeyError:
            coords_list_dicts[pair] = {'query':[], 'host':[]}
            coords_list_dicts[pair]['query'] = [query_range]

        hstart = (int(ln[8]) - 1)
        hend = int(ln[9]) + 1
        host_range = [hstart, hend]

        try:
            coords_list_dicts[pair]['host'].append(host_range)
        except KeyError:
            coords_list_dicts[pair]['host'] = [host_range]

for pair in coords_list_dicts.keys():
    length = len(coords_list_dicts[pair]['query'])
    coords_list_dicts[pair]['length'] = len(coords_list_dicts[pair]['query'])

def is_overlaping(a, b):
    if (b[0] >= a[0] and b[0] <= a[1]):
        return True
    elif (a[0] >= b[0] and a[0] <= b[1]):
        return True
    elif (a[1] >= b[0] and a[1] <= b[1]):
        return True
    elif (b[1] >= a[0] and b[1] <= a[1]):
        return True
    else:
        return False

def merge_range(range1, range2):
    range_start = min(min(range1), min(range2))
    range_end = max(max(range1), max(range2))
    new_range = (range_start, range_end)
    return new_range

def merge_ranges_dict(coords_dict):
    merged_coords_dict = {}
    for pair in coords_dict.keys():

        query_lists = coords_dict[pair]['query']
        host_lists = coords_dict[pair]['host']
        range_num = coords_dict[pair]['length']

        ## assign variables and add first range to range_lists
        q_range = tuple(query_lists[0])
        h_range = tuple(host_lists[0])

        query_range_list = [q_range]
        host_range_list = [h_range]
        len_query_range_list = len(query_range_list)

        i = 0 ## position of range to compare from coords_dict
        while i < range_num:
            j = 0 ## position of range to compare from range_lists

            while j < len_query_range_list: ## check ranges added to range_lists
                if (q_range == query_range_list[j]) and (h_range == host_range_list[j]):
                    ## if exact same range already in range_lists -- go to next in query_lists
                    if i < range_num:
                        q_range = tuple(query_lists[i])
                        h_range = tuple(host_lists[i])
                        i += 1
                    break

                else:
                    any_overlap = []
                    for prev in range(len(query_range_list)):
                        tf = is_overlaping(q_range, query_range_list[prev]) and is_overlaping(h_range, host_range_list[prev])
                        any_overlap.append(tf)

                    if any(any_overlap):
                        ## if any overlapping ranges are found, find index from query_range_list for overlaps
                        idx = sorted([n for n, o in enumerate(any_overlap) if o == True], reverse=True)

                        ## get overlapping ranges and remove them from range_lists
                        q_ovrlp = [query_range_list.pop(n) for n in idx]
                        h_ovrlp = [host_range_list.pop(n) for n in idx]

                        ## get start and stop values from overlapping ranges
                        qmin = min(q_ovrlp)[0]
                        qmax = max(q_ovrlp, key = lambda t: t[1])[1]
                        q_ovrlp = (qmin, qmax)

                        hmin = min(h_ovrlp)[0]
                        hmax = max(h_ovrlp, key = lambda t: t[1])[1]
                        h_ovrlp = (hmin, hmax)

                        ## overwrite range variables with new start/stop coordinates and add to range_lists
                        q_range = merge_range(q_range, q_ovrlp)
                        h_range = merge_range(h_range, h_ovrlp)
                        query_range_list.append(q_range)
                        host_range_list.append(h_range)
                        ## update length variable of range_lists (j)
                        len_query_range_list = len(query_range_list)

                    else:
                        ## no overlaps in current range_lists, append current ranges
                        query_range_list.append(q_range)
                        host_range_list.append(h_range)

                    if i < range_num:
                        q_range = tuple(query_lists[i])
                        h_range = tuple(host_lists[i])
                        i += 1
                    j += 1
            len_query_range_list = len(query_range_list)

        paired_coords_set = set()

        for item in zip(query_range_list, host_range_list):
            paired_coords_set.add(item)

        for item in paired_coords_set:
            query_range = [item[0][0], item[0][1]]
            host_range = [item[1][0], item[1][1]]
            try:
                merged_coords_dict[pair]['query'].append(query_range)
                merged_coords_dict[pair]['host'].append(host_range)

            except KeyError:
                merged_coords_dict[pair] = {'query':[], 'host':[]}
                merged_coords_dict[pair]['query'] = [query_range]
                merged_coords_dict[pair]['host'] = [host_range]

    for pair in merged_coords_dict.keys():
        merged_coords_dict[pair]['length'] = len(merged_coords_dict[pair]['query'])

    return merged_coords_dict

def region_avg_qsasa(ppops_files, kmer_coords):
    colnames = ['ResidNe', 'Chain', 'ResidNr', 'iCode', 'Phob/A^2', 'Phil/A^2', 'SASA/A^2', 'Q(SASA)', 'N(overl)', 'Surf/A^2']
    region_vals_lists = []
    no_hpops = set()
    p = 0
    h = 0

    for pfile in ppops_files:
        ## check each pathogen pops file, get protein name
        prot_name = pfile.split('pops_')[1]
        prot_name = prot_name.split('.')[0]
        prot_name = prot_name.split('-')[0]
        ## get all pairs that include pathogen protein name
        all_pair_keys = [pair_key for pair_key, val in kmer_coords.items() if prot_name in pair_key]
        if len(all_pair_keys) > 0:
            ## if there is at least one pathogen-host pair for that protein, check number of pops files
            pfiles = [file for file in ppops_files_all if prot_name in file]
            if len(pfiles) > 1:
                print(f"Multiple pops files for pathogen protein {prot_name}")

            ## open pops file, make dataframe
            prot_sasa = pd.read_csv(f'{pfile}', sep="\s+", header=None, engine='python', skiprows=3, skipfooter=3)
            prot_sasa.columns = colnames
            for pair in all_pair_keys:
                q_region_vals_lists = []
                h_region_vals_lists = []
                ## get host protein file, need to add check if exists
                h_prot_name = pair.split('.')[1]
                ## find all pops files
                hfiles = [file for file in hpops_files_all if h_prot_name in file]
                if len(hfiles) == 0:
                    no_hpops.add(h_prot_name)
                    #print(f"No pops files for host protein {h_prot_name}")
                elif len(hfiles) > 1:
                    print(f"Multiple pops files for host protein {h_prot_name}")

                ## assign variables
                query_ranges = kmer_coords[pair]['query']
                host_ranges = kmer_coords[pair]['host']
                for i in range(len(query_ranges)):
                    q_range = query_ranges[i]
                    h_range = host_ranges[i]

                    q_region_vals = prot_sasa.iloc[q_range[0]:q_range[1]]['Q(SASA)'].to_list()
                    q_region_len = (q_range[1] - q_range[0])
                    if not q_region_vals:
                        print(f"Pathogen protein {prot_name} region {q_range[0]}:{q_range[1]} out of pops file range")
                        p += 1
                        q_region_vals = 0
                        q_region_avg = "NaN"
                    else:
                        q_region_len = len(q_region_vals) #(q_range[1] - q_range[0])
                        q_region_avg = (sum(q_region_vals) / q_region_len)
                    q_region_vals_lists.append(q_region_vals)
                    

                    h_region_len = (h_range[1] - h_range[0])
                    if len(hfiles) != 0:
                        hfile = f"{hpops_path}/pops_{h_prot_name}.out"
                        h_prot_sasa = pd.read_csv(f'{hfile}', sep="\s+", header=None, engine='python', skiprows=3, skipfooter=3)
                        h_prot_sasa.columns = colnames
                        h_region_vals = h_prot_sasa.iloc[h_range[0]:h_range[1]]['Q(SASA)'].to_list()
                        if not h_region_vals:
                            print(f"Host protein {h_prot_name} region {h_range[0]}:{h_range[1]} out of pops file range")
                            h += 1
                            h_region_vals = 0
                            h_region_avg = "NaN"
                        else:
                            h_region_len = len(h_region_vals)
                            h_region_avg = (sum(h_region_vals) / h_region_len)                    
                        h_region_vals_lists.append(h_region_vals)


                    else: 
                        h_region_vals = 0
                        h_region_vals_lists.append(h_region_vals)
                        h_region_avg = "NaN"

                    region_vals_lists.append([prot_name, q_region_avg, q_region_len, q_range[0], (q_range[1] - 1), 
                           h_prot_name, h_region_avg, h_region_len, h_range[0], (h_range[1] - 1)])

    return region_vals_lists, p, h

merged_coords_dict = merge_ranges_dict(coords_list_dicts)
region_vals_lists, p, h = region_avg_qsasa(ppops_files, merged_coords_dict)

all_qsasa_df = pd.DataFrame(region_vals_lists, columns = ['q_prot', 'q_avg', 'q_length', 'q_start', 'q_end', 'h_prot', 'h_avg', 'h_length', 'h_start', 'h_end'])
all_qsasa_df.to_csv(f"{fileid}_paired_regions.tsv", sep="\t", header=False, index=False)

qsasa_df_filtered = all_qsasa_df[(all_qsasa_df['h_avg'] != 'NaN') & (all_qsasa_df['q_avg'] != 'NaN')]
qsasa_df_filtered = qsasa_df_filtered[(qsasa_df_filtered['q_avg'] > min_qsasa) & (qsasa_df_filtered['h_avg'] > min_qsasa)] 
min_qsasa_str=str(min_qsasa).replace('.', '')
qsasa_df_filtered.to_csv(f"{fileid}_paired_qsasa{min_qsasa_str}.tsv", sep="\t", header=False, index=False)

print(f"\n\n{all_qsasa_df.shape[0]} paired regions saved to {fileid}_paired_regoins.tsv")
print(f"{qsasa_df_filtered.shape[0]} paired regions with minimum QSASA of {min_qsasa} saved to {fileid}_paired_qsasa{min_qsasa_str}.tsv")

# pathogen_blast_prots = set()
# host_blast_prots = set()
# for pair in merged_coords_dict.keys():
#     pathogen_blast_prots.add(pair.split('.')[0])
#     host_blast_prots.add(pair.split('.')[1])

# pathogen_qsasa_prots = set(all_qsasa_df['q_prot'].to_list())
# all_pathogen_prots = [file.split('_')[1] for file in ppops_files]
# all_pathogen_prots = set([file.split('.')[0] for file in all_pathogen_prots])
# no_pqsasa = set(pathogen_blast_prots.difference(pathogen_qsasa_prots))
# no_ppops = set(pathogen_blast_prots.difference(all_pathogen_prots))


# with open(f"{fileid}_missing_pathogen_files.tsv", "w") as missing_ppops:
#     for name in no_ppops:
#         missing_ppops.write(f"{name}\n")

# host_qsasa_prots = set(all_qsasa_df['h_prot'].to_list())
# all_host_prots = [file.split('_')[1] for file in hpops_files]
# all_host_prots = set([file.split('.')[0] for file in all_host_prots])
# no_hqsasa = set(host_blast_prots.difference(host_qsasa_prots))
# no_hpops = set(host_blast_prots.difference(all_host_prots))

# with open(f"{fileid}_missing_host_files.tsv", "w") as missing_hpops:
#     for name in no_hpops:
#         missing_hpops.write(f"{name}\n")

# print(f"{len(no_ppops)} {p_name} proteins missing pops file, {p} regions out of pops file range")
# print(f"{len(no_hpops)} {h_name} proteins missing pops file, {h} regions out of pops file range")
# print(f"Names of proteins missing pops files saved to {fileid}_missing_pathogen_files.tsv and {fileid}_missing_host_files.tsv\n")
