from RNAtools import CT
import argparse

def parseArgs():
    """
    function to parse command line parameters
    """

    # init argument parser object
    prs = argparse.ArgumentParser(description='Convert dotbracket format to ct')

    # add some options
    prs.add_argument('infile', action='store', type=str,
                     help='The name of the dot-bracket file to convert.')
    prs.add_argument('outfile', action='store', type=str,
                     help='The destination CT filename.')

    # parse and return
    out=prs.parse_args()

    return out.infile, out.outfile

def main():

    infile, outfile = parseArgs()

    with open(infile, "rU") as f:
        name = ""
        sequence = ""
        dotbracket = ""

        # read the file
        for i, line in enumerate(f):
            line = line.rstrip()
            if i == 0:
                name = line.lstrip(">")
            if i == 1:
                sequence = line
            if i == 2:
                dotbracket = line

        print(name)
        print(sequence)
        print(dotbracket)

        # construct the ct object
        rna = CT()

        # give it the dotbracket line
        rna.dot2ct(sequence, dotbracket, name=name)

        # print the pairlist
        print(rna.pairList())

        # write it to disk
        rna.writeCT(outfile)


if __name__ == "__main__":
    main()
