from RNAtools import CT
import argparse

def parseArgs():
    """
    function to parse command line parameters
    """

    # init argument parser object
    prs = argparse.ArgumentParser(description='Clean up the input CT file by\
                removing non-canonical base pairs and 1 bp helices')

    # add some options
    prs.add_argument('infile', action='store', type=str,
                     help='The name of the CT file to clean.')
    prs.add_argument('outfile', action='store', type=str,
                     help='The destination CT filename.')

    # parse and return
    out=prs.parse_args()

    return out.infile, out.outfile

def main():

    infile, outfile = parseArgs()

    infile_ct = CT(infile)
    outfile_ct = infile_ct.cleanCT()

    outfile_ct.writeCT(outfile)

if __name__ == "__main__":
    main()
