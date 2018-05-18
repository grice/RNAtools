##########################################
#    Varna SVG color script
#
#      Author: Gregg Rice
#              gmr@unc.edu
#
# Affiliation: Weeks Lab, UNC Chemistry
#
#        Date: Feb 24 2015
#     Version: 0.92
#
# released under GPL 2.0
##########################################

simple plotting script to make SVG files compatable
with illustrator from common RNA structure formats

see also run_examples.sh


1a) The simplest case. VARNA file with SHAPE information embedded

python colorVSVG.py 5S.VARNA foo.svg

1b) XRNA file with SHAPE information embedded. Add the -x flag

python colorVSVG.py 16SrRNA_1X1M7.xrna -x foo.svg

2) Fancy formatting options...


change the nucleotide number start
--offset N (default = 1)


compare the structure with a predicted CT file. RNAstructure compatable CT format
--ct FILE


Add SHAPE coloring information. Will override values stored in the input file
-s SHAPE


Add circles behind the nucleotides from a defined color file
-c circlefile.txt

format of the circlefile:
#nuc outline(RGB) fill(RGB
1 0,0,0 150,150,150
2 0,0,0 150,150,150
3 0,0,0 150,150,150
4 0,0,0 150,150,150
5 0,0,0 150,150,150
...etc


Add enzymatic probing information
-e enzymefile.txt

format of the enzymefile, first line has a header
# 0=no arrow, 1=lo, 2=mid, 3=high
nucleotide reactivity
11 3
12 2
45 1
46 2
63 1
64 2

Add extra lines,
-l linesfile.txt
format of the linesfile, no header, 
has three columns: fromNuc toNuc correlationValue
1 45 0.5
1 60 0.5
45 60 -0.5
