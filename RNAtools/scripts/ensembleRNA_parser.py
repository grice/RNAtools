from RNAtools import CT
from Bio import SeqIO
import sys

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("USAGE: python ensembleRNA_parser.py ref.fa ensemble.csv outDir/")  # noqa
        sys.exit()

    # read the fasta file
    with open(sys.argv[1], "r") as f:
        x = list(SeqIO.parse(f, "fasta"))
        name = x[0].id
        seq = x[0].seq
        print(name, seq)

    # read the ensemble csv file
    # output is fasta_id + _ + fraction
    with open(sys.argv[2], "r") as g:
        for n, line in enumerate(g):
            if n == 0:
                continue
            spline = line.rstrip().split(",")
            frac = spline[2]

            # filter for ensembles that are at least 5%
            if int(frac) < 50:
                continue

            dotbracket = spline[3]

            # make a ct object, read dotbracket, write
            x = CT()
            x.dot2ct(seq, dotbracket, name=name)
            x.writeCT("{2}/{0}_{1}.ct".format(name, frac, sys.argv[3]))

            print(frac, x.pairList())
