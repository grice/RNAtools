from RNAtools import dotPlot, readSHAPE
import matplotlib.cm as cm
import matplotlib.colors as colors

import argparse

def parseArgs():
    """
    function to parse command line parameters
    """

    # init argument parser object
    prs = argparse.ArgumentParser(description='Creates a colorfile for colorRNA from dp and two column formats')

    # add some options
    prs.add_argument('infile', action='store', type=str,
                     help='.dp partition file or .txt')
    prs.add_argument('outfile', action='store', type=str,
                     help='.colorfile for colorRNA')
    prs.add_argument('--palette', default="viridis", type=str,
                     help="Matplotlib color palette")
    prs.add_argument("--vmin", default=0.0, type=float,
                     help='Min value for color scale')
    prs.add_argument("--vmax", default=1.0, type=float,
                     help='Max value for color scale')
    prs.add_argument('--twoCol', action="store_true",
                     help="input is two column .shape format")

    # parse and return
    out=prs.parse_args()

    return out.infile, out.outfile, out.palette, out.vmin, out.vmax, out.twoCol

def main():

    infile, outfile, palette, vmin, vmax, twoCol = parseArgs()

    if twoCol:
        data = readSHAPE(infile)
    else:
        data = dotPlot(infile).calcShannon()
    
    norm = colors.Normalize(vmin=vmin, vmax=vmax)
    cmap = cm.get_cmap(palette)

    line = ""
    for i, val in enumerate(data):
        s = list(map(lambda x: str(int(x*255)), cmap(norm(val))))
        line += "{0} {1} {1}\n".format(i+1, ",".join(s[:3]))
    
    w = open(outfile, "w")
    w.write(line)
    w.close()

if __name__ == "__main__":
    main()
