
from RNAtools import CT
import sys


x = CT(sys.argv[1])
x.cutCT(1,int(sys.argv[2])).writeCT(sys.argv[3])
