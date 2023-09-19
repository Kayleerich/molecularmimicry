import argparse
from pathlib import Path
from Bio import SeqIO

parser = argparse.ArgumentParser(description='Create sequence fragments from FASTA')

parser.add_argument('-i', '--input', metavar='INFILE', action='store', type=str, required=True, help='Input fasta file')
parser.add_argument('-k', '--kmer_size', metavar='NUM', action='store', type=int, required=True, help='Sequence fragment/k-mer size')
parser.add_argument('-n', '--num_masked', metavar='NUM', action='store', type=int, help='Max number of masked residues per k-mer, default = k/2')
parser.add_argument('-m', '--mask', metavar='CHAR', action='store', type=str, default='X', help='Mask character, default = X')
parser.add_argument('-s', '--silent', action='store_true', help='Silence terminal output') 

args = parser.parse_args()

file = Path(args.input)
k = args.kmer_size
if args.num_masked == None:
    n = int(k/2)
else: 
    n = args.num_masked
mask = args.mask

if args.silent == False:
    print(f"Input file: {file.name}\nK-mer size: {k}\nMaximum number of masked residues per k-mer: {n}\nMasking character: {mask}")


def fragment_fasta(file, k, n = int(k)/2, mask = "X"):
    """
    Read fasta file and create k-mers for each sequence.
    Can take hard-masked sequences and remove k-mers that exceed the threshold given for number of masked redisues. 
    By default k-mers containing more than half masked residues will be removed, assumes masking character to be "X".
    """

    k = int(k)
    n = int(n)

    with open(f"{file.stem}-{str(k)}mers.fasta", 'w') as fragfile:
        records = SeqIO.parse(file, 'fasta')
        for record in records:
            sequence = record.seq
            for i in range(0, len(sequence)):
                seq_fragment = sequence[i:i+k]
                header = f"{record.description}/{i}:{i+k}"
                record.id = header
				        #remove fragments that have too many masked characters
                if seq_fragment.count(mask) > n:
                    continue
                #remove fragments that are less than 14aa long
                elif len(seq_fragment) != k:
                    continue
                else:
                    fragfile.write(f">{header}\n{seq_fragment}\n")

fragment_fasta(file, k, n, mask)

if args.silent == False:
    print(f"\nDone, k-mer file saved as: {file.stem}-{str(k)}mers.fasta\n")
