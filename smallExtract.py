from RNAtools import CT

import sys

rna = CT(sys.argv[1])
pk1,pk2 = rna.extractPK()

print pk1, pk2
