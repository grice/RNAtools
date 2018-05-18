# GPL 2.0
#
# written by Gregg Rice, gmr@unc.edu, rice.gregg@gmail.com
# copywrite 2013-2015
# all rights reserved
#

##################################################################################
# render one or two CT secondary structures with helices shown as semicircular arcs
# This simple version of the script does not plot any additional data above the structure(s).
# functionsplotArcRibbons by Steve Busan, modified in 2014 with permission for python hooks Gregg Rice
# om function arcplot
# ----

import RNAtools as RNA
import numpy as np
import sys, argparse, os, traceback, copy, re, math

from matplotlib import rc
#rc('text',usetex=True)

def rgb_int2pct(rgbArr):
    """
    converts RGB 255 values to a 0 to 1 scale

    (255,255,255) --> (1,1,1)
    """
    out = []
    for rgb in rgbArr:
        out.append((rgb[0]/255.0, rgb[1]/255.0, rgb[2]/255.0))
    return out

def parseArgs():
    prs = argparse.ArgumentParser()
    prs.add_argument("dotplot", type=str, help='input dotplot file to get pairing probabilities')
    prs.add_argument("outputPDF",type=str, help="Name of the output PDF graphic")
    prs.add_argument("--referenceCT",type=str, help="reference ct file used for getting the RNA sequence")
    prs.add_argument("--pkDS", type=str, help="Double stranded pk file to optionally add in black arcs for pseudoknots. Two column space seperated text file with pairs in columns")
    prs.add_argument("--secondaryStructure", action="store_true", default=False, help="Plot the secondary structure on top")

    out = prs.parse_args()

    return out

def plotArcRibbons(pairedNuc, ax, color, alpha=1.0, flip=False):
    """
    pairedNuc is a vector length N (seqeunce length) with non-zero values setting the connections
    this is a non-redundant .ct column. Example below
    0 3 0 0 0
      |___|   <-- arc between pos 2 and 4

    ax is the plot axis object

    color is the color of the arc, and alpha is the transparency of the arc

    flip will plot the arcs as rainbows rather than smiles
    """

    from matplotlib.path import Path
    import matplotlib.patches as patches

    handleLengthFactor = 4*(math.sqrt(2)-1)/3
    #vert = 100.0

    i = 0
    while i < len(pairedNuc):
        if pairedNuc[i] > i+1:
            outerPair = [i+0.5,pairedNuc[i]+0.5]
            # find the right side of helix
            lastPairedNuc = pairedNuc[i]
            offset = 1
            while i+offset<len(pairedNuc) and abs(pairedNuc[i+offset]-lastPairedNuc) == 1:
                lastPairedNuc = pairedNuc[i+offset]
                offset += 1
            innerPair = [i+offset+0.5, pairedNuc[i+offset-1]-0.5]
            i += offset-1
            outerRadius = (outerPair[1]-outerPair[0])/2.0
            innerRadius = (innerPair[1]-innerPair[0])/2.0
            #print "innerPair %s, outerPair %s"%(str(innerPair),str(outerPair))
            verts = [
            (outerPair[0], 0), # outer left

            (outerPair[0], -handleLengthFactor*outerRadius), # outer left control 1
            (outerPair[0]+outerRadius-handleLengthFactor*outerRadius, -outerRadius), # outer left control 2
            (outerPair[0]+outerRadius, -outerRadius), # outer center
            (outerPair[0]+outerRadius+handleLengthFactor*outerRadius, -outerRadius), # outer right control 1
            (outerPair[1], -handleLengthFactor*outerRadius), # outer right control 2

            (outerPair[1], 0), # outer right
            (innerPair[1], 0), # inner right

            (innerPair[1], -handleLengthFactor*innerRadius), # inner right control 1
            (innerPair[0]+innerRadius+handleLengthFactor*innerRadius, -innerRadius), # inner right control 2
            (innerPair[0]+innerRadius, -innerRadius), # inner center
            (innerPair[0]+innerRadius-handleLengthFactor*innerRadius, -innerRadius), # inner right control 1
            (innerPair[0], -handleLengthFactor*innerRadius), # inner right control 2

            (innerPair[0], 0), # inner left
            (outerPair[0], 0) # outer left duplicate point
            ]

            if flip:
                flip_offset=2

                verts = [
                (outerPair[0], 0 + flip_offset), # outer left

                (outerPair[0], handleLengthFactor*outerRadius + flip_offset), # outer left control 1
                (outerPair[0]+outerRadius-handleLengthFactor*outerRadius, outerRadius + flip_offset), # outer i left control 2
                (outerPair[0]+outerRadius, outerRadius + flip_offset), # outer center
                (outerPair[0]+outerRadius+handleLengthFactor*outerRadius, outerRadius + flip_offset), # outer right control 1
                (outerPair[1], handleLengthFactor*outerRadius + flip_offset), # outer right control 2

                (outerPair[1], 0 + flip_offset), # outer right
                (innerPair[1], 0 + flip_offset), # inner right

                (innerPair[1], handleLengthFactor*innerRadius + flip_offset), # inner right control 1
                (innerPair[0]+innerRadius+handleLengthFactor*innerRadius, innerRadius + flip_offset), # inner right control 2
                (innerPair[0]+innerRadius, innerRadius + flip_offset), # inner center
                (innerPair[0]+innerRadius-handleLengthFactor*innerRadius, innerRadius + flip_offset), # inner right control 1
                (innerPair[0], handleLengthFactor*innerRadius + flip_offset), # inner right control 2

                (innerPair[0], 0 + flip_offset), # inner left
                (outerPair[0], 0 + flip_offset) # outer left duplicate point
                ]

            for n in range(len(verts)):
                verts[n] = [verts[n][0], verts[n][1]-1.2]#/vert-0.013]
            codes = [
            Path.MOVETO,

            Path.CURVE4,
            Path.CURVE4,
            Path.CURVE4,
            Path.CURVE4,

            Path.LINETO,
            Path.LINETO,
            Path.LINETO,

            Path.CURVE4,
            Path.CURVE4,
            Path.CURVE4,
            Path.CURVE4,

            Path.LINETO,
            Path.LINETO,
            Path.CLOSEPOLY,
            ]
            path = Path(verts, codes)
            patch = patches.PathPatch(path, facecolor=color, linewidth=0, edgecolor='none', alpha=alpha)
            ax.add_patch(patch)
            patch.set_clip_on(False)
        i += 1

def arcplot(outPath="arcs.pdf",title="",seq=["A"],pairedNucArr=[], arcColors = [],alpha=[1.0], maxDistance=None, ct=None, secondary_structure=False):
    """
    draws arcs for many sets of nucleotides, pairedNucArr is a 2D array containing all sets of plotting
    elements. arcColor is the same length as pairedNucArr but contains color strings

    pairedNucArr is M x N where M is the number of plotted elements and N is the length of the RNA to plot
    arcColors is a M x 3 array of RGB tuples for each element
    alpha is a M x 1 array of alpha values to plot

    secondary structure is a booleon that sets whether to plot the structure in the ct file above the arcs representation
    seq is length N

    function duplicated by gmr for plotting several structures on the same plot
    """

    def findMaxDistance(pairingArrays):
        """
        finds the maximum pairing distance from a 2D array of all arc elements
        """
        maxDistance = 0

        for pairedNucA in pairingArrays:
            for i in range(len(pairedNucA)):
                fromNuc = i+1
                toNuc = pairedNucA[i]
                if toNuc == 0:
                    toNuc = fromNuc
                dist = toNuc-fromNuc
                if dist > maxDistance:
                    maxDistance = dist

        return maxDistance

    import matplotlib as mp
    mp.use('Agg')
    mp.rcParams['xtick.major.size'] = 8
    mp.rcParams['xtick.major.width'] = 2.5
    mp.rcParams['xtick.direction'] = 'out'
    mp.rcParams['xtick.minor.size'] = 4
    mp.rcParams['xtick.minor.width'] = 1


    import matplotlib.pyplot as plot
    import matplotlib.patches as patches
    import matplotlib.gridspec as gridspec

    num = range(1,len(seq)+1)

    # adjust the scale factor to fit on a page --gmr
    if len(num)>10000:
        scaleFactor = 0.0005
    if len(num)>500:
        scaleFactor = 0.005
    else:
        scaleFactor = 0.05

    # find longest base-pair and scale height of plot to fit this arc
    if not maxDistance:
        maxDistance = findMaxDistance(pairedNucArr)

    # set up figure dimensions
    figWidth = len(seq)*scaleFactor
    figHeight = maxDistance/2.0*scaleFactor

    # double the figure height if plotting the secondary structure above the plot
    if secondary_structure:
        figHeight = maxDistance * scaleFactor

    fig = plot.figure(figsize=(figWidth, figHeight)) # 500*scaleFactor


    if not secondary_structure:
        ax2 = plot.subplot(111)
    else:
        ax2 = plot.subplot(212)

    plot.xlim(0,len(seq))
    #ax2 = plot.gca()

    # ticks on top
    ax2.get_xaxis().tick_top()

    # determine tick locations based on sequece length
    from matplotlib.ticker import MultipleLocator, FormatStrFormatter
    if len(seq) <= 500:
        majorLocator = MultipleLocator(100)
        minorLocator = MultipleLocator(10)
        interval = 10
    elif len(seq) <= 10000:
        majorLocator = MultipleLocator(500)
        minorLocator = MultipleLocator(100)
        interval = 5
    else:
        majorLocator = MultipleLocator(2500)
        minorLocator = MultipleLocator(500)
        interval = 5

    majorFormatter = FormatStrFormatter('%i')
    minorFormatter = FormatStrFormatter('%i')

    ax2.xaxis.set_major_locator(majorLocator)
    ax2.xaxis.set_major_formatter(majorFormatter)
    ax2.xaxis.set_minor_locator(minorLocator)
    ax2.xaxis.set_minor_formatter(minorFormatter)


    plot.subplots_adjust(hspace=0)
    plot.ylim((-float(maxDistance)/2.0,0))
    #plot.ylim((-1.0,0))

    ax2.set_frame_on(False)
    ax2.axes.get_yaxis().set_visible(False)
    #ax2.axes.get_xaxis().tick_bottom()
    xlabels = ax2.axes.get_xaxis().get_majorticklabels()
    xlabels[0].set_visible(False)
    for label in xlabels:
        label.set_weight('bold')
        label.set_size(14)
        #label.set_rotation(30)
    xlabels = ax2.axes.get_xaxis().get_minorticklabels()
    xlabels[0].set_visible(False)
    labelCount = 0
    for label in xlabels:
        label.set_size(7)
        #label.set_rotation(30)
        # also need to hide minor tick labels that overlap major tick labels
        if labelCount%interval==1:
            label.set_visible(False)
        labelCount += 1

    xticks = ax2.axes.get_xaxis().get_major_ticks()
    xticks[0].set_visible(False)
    xticks = ax2.axes.get_xaxis().get_minor_ticks()
    xticks[0].set_visible(False)


    #bothColor = "green"
    #aColor = "red"
    #bColor = "purple"

    # plot the arcs
    for arc in range(len(pairedNucArr)):
        #if arc % 500 == 0:
        #    print arc, len(pairedNucArr)
        plotArcRibbons(pairedNucArr[arc], ax2, arcColors[arc], alpha=alpha[arc])

    #plotArcRibbons(aOnlyPaired, ax2, aColor, alpha=alpha)
    #plotArcRibbons(bOnlyPaired, ax2, bColor, alpha=alpha)
    #plotArcRibbons(bothPaired, ax2, bothColor, alpha=alpha)

    xmax, xmin, ymin, ymax = plot.axis()

    # put nuc sequence on axis
    if len(seq) <= 500:
        fontProp = mp.font_manager.FontProperties(family = "monospace",
                                              style="normal",
                                              weight="extra bold",
                                              size="4")
        for i in range(len(seq)):
            nuc = seq[i]
            if nuc == "T":
                nuc = "U"
            col = "black"
            plot.annotate(nuc, xy=(i+0.5,ymax),fontproperties=fontProp,color=col,annotation_clip=False,verticalalignment="top")

    if secondary_structure:
        ax1 = plot.subplot(211)
        ax1.set_frame_on(False)
        ax1.axes.get_yaxis().set_visible(False)
        ax2.axes.get_xaxis().set_visible(False)
        plot.xlim(0,len(seq))
        plot.ylim(0,maxDistance/2.0)

        ax1.spines["bottom"].set_position(("axes", -0.025))

        # get the helices from the ct file
        helix = ct.extractHelices(fillPairs=False)

        plot_set = []
        plot_color = []
        ss_color = (0.8,0.8,0.8)

        for k in helix:
            temp = np.zeros_like(ct.ct)
            for nt in helix[k]:
                temp[nt[0]-1] = nt[1]
                plot_set.append(temp)
                plot_color.append(ss_color)
            #print(helix[k])
        #print(helix)
        for hel, col in zip(plot_set, plot_color):
            plotArcRibbons(hel, ax1, col, alpha=1, flip=True)
        #plotArcRibbons(helix, ax1, (0.8,0.8,0.8), alpha=1)
        # xticks = ax1.axes.get_xaxis().get_major_ticks()
        # for k in xticks:
        #     k.set_visible(False)
        # xticks = ax1.axes.get_xaxis().get_minor_ticks()
        # for k in xticks:
        #     k.set_visible(False)

    plot.savefig(outPath,dpi=100,bbox_inches="tight")



def splitPlot(dpObj, ctObj, pk=None, outFile="arcs.pdf", secondary_structure=False):
    """
    splitplot takes a RNAtools dotplot object and a ct object to generate a the figure

    pk is a list pseudoknotted pairs to plot in the bottom panel

    outfile is the destination for the matplotlib output, usually a pdf file or a png

    secondary_structure is a booleon that sets whether to plot the secondary structure above the arc plot
    """
    x = dpObj
    y = ctObj

    # binning is in log10 scale
    #binning = [0.0,0.09691,0.5228,1.0,2.0]
    binning = [1.5228, 1.0, 0.5228, 0.09691, 0.0]
    # 3% 10% 30% 80% 100%

    alphaList = [0.7, 0.7, 0.7, 0.3]
    alphaList = [0.3, 0.7, 0.7, 0.7]
    alphaList = [1.0,1.0,1.0,1.0,]
    #alphaList = [1.0, 1.0, 1.0, 1.0]
    #colorList  = ["red", "orange", "yellow", "green","blue", "violet"]
    #colorList  = [(215, 25, 28), (253, 174, 97), (171, 221, 164), (43, 131, 186)]
    colorList  = [ (43, 131, 186), (171, 221, 164), (253, 174, 97), (215, 25, 28)]
    colorList  = [ (150,150,150), (255,204,0),  (72,143,205) ,(81, 184, 72)  ]
    colorList = rgb_int2pct(colorList)

    nucArr = []
    colors = []
    alpha  = []

    # bin the pairs by cutoff
    for i in range(0, len(binning)-1):

        probPairs = x.requireProb(binning[i],binning[i+1]).pairList()

        for pair in probPairs:

            temp = np.zeros_like(y.ct)
            temp[pair[0]-1] = pair[1]

            #tempCT = RNA.CT()
            #tempCT.pair2CT([pair],y.seq)

            #nucArr.append(tempCT.stripCT())
            nucArr.append(temp)

            #add a color from the choice list
            colors.append(colorList[i])
            alpha.append(alphaList[i])

    if pk:
        for pair in pk:
            temp = np.zeros_like(y.ct)
            temp[pair[0]-1] = pair[1]
            nucArr.append(temp)
            colors.append((0,0,0))
            alpha.append(0.8)
    #nucArr.append(y.stripCT())
    #colors.append("gray")

    arcplot(outPath=outFile, pairedNucArr=nucArr, arcColors=colors,seq=y.seq, alpha=alpha, maxDistance=None, secondary_structure=secondary_structure, ct=y)

def readPKfile(fIN):
    out = []
    for i in open(fIN).readlines():
        if i.lstrip()[0] == '#':
            continue
        pk = map(int, i.rstrip().split())
        out.append((pk[0],pk[1]))
    return out


def main():
    # parse the command line options
    arg = parseArgs()

    # read in the reference files
    x = RNA.dotPlot( arg.dotplot )
    if arg.referenceCT:
        y = RNA.CT( arg.referenceCT )
    else:
        # if we dont have a sequence create an empty CT file
        y = RNA.CT()
        # make the sequence all N as a placeholder
        seq = "N"*x.length
        y.pair2CT([], seq)

    #print(y)


    # add correction for slipped base pairs
    #x.averageSlippedBPs(y,predictedOnly=True)

    # read in the pseudoknots file if it exists
    if arg.pkDS:
        pkntsFile = arg.pkDS
        pknts = readPKfile(pkntsFile)
    else:
        pknts = None

    # main plotting loop
    splitPlot(x, y, pk=pknts, outFile=arg.outputPDF, secondary_structure=arg.secondaryStructure)

if __name__ == '__main__':
    main()
