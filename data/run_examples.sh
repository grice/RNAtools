
# simple run, no color options
colorRNA 5S.VARNA foo_simple.svg

# input file is an xRNA file, -x flag
colorRNA 16SrRNA_1X1M7.xrna foo_16S.svg

# compare two structures
colorRNA 5S.VARNA foo_compare.svg --ct 5S.correct.ct

# add colored circles behind each nucleotide, color is defined in the -c file
colorRNA 5S.VARNA foo_circles.svg -c 5S.gray.txt 

# add enzyme probing data
colorRNA 5S.VARNA foo_enzyme.svg -e 5S.enzyme.txt

# add correlation lines, eg from a DMS correlation experiment
colorRNA 5S.VARNA foo_corr.svg -l 5S.lines.txt

# run with lots of options
colorRNA 5S.VARNA foo_complex_plotting.svg -c 5S.gray.txt -e 5S.enzyme.txt -l 5S.lines.txt
