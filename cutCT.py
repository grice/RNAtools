
from RNAtools import CT
import sys

if len(sys.argv)<5:
    print 'Usage: inp.ct out.ct start end'
    sys.exit()

x = CT(sys.argv[1])
x.cutCT( int(sys.argv[3]), int(sys.argv[4]) ).writeCT(sys.argv[2])

