from RNAtools import CT

def extractPK(rna):
    """
    returns the pk1 and pk2 pairs from a CT file. Ignores single base
    pairs. PK1 is defined as the helix crossing the most other helices.
    if there is a tie, the most 5' helix called pk1
    
    returns pk1 and pk2 as a list of base pairs e.g [(1,10),(2,9)...
    """
    
    def checkOverlap(h1,h2):
        # only need to check one set of pairs from each
        # of the helicies. Test is to see if they form
        # a cross hatching pattern
        if max(h1[0]) > min(h2[0]) > min(h1[0]):
            if max(h2[0]) > max(h1[0]):return True
        if max(h2[0]) > min(h1[0]) > min(h2[0]):
            if max(h1[0]) > max(h2[0]):return True
        else: return False
    
    # first step is to find all the helicies
    helicies = {}
    nt = 0
    heNum = 0
    
    while nt < len(rna.ct):
        if rna.ct[nt] == 0:
            nt += 1
            continue
        else:
            # skip dups
            if nt > rna.ct[nt]:
                nt += 1
                continue
            
            tempPairs = []
            
            while rna.ct[nt] != 0:
                tempPairs.append((nt+1,rna.ct[nt]))
                nt +=1
            
            # remove single bp helicies
            if len(tempPairs) <= 1:
                continue
            helicies[heNum] = tempPairs
            heNum += 1
    
    # do the helicies overlap? Check for a crosshatching pattern
    # append them to a list if they have it.
    overlaps = [] # stores the helix number
    
    for i in xrange(heNum):
        for j in xrange(i+1,heNum):
            if checkOverlap(helicies[i],helicies[j]):
                overlaps.append((i,j))
    
    # determine which is pk1
    allHelix = []
    for i,j in overlaps:
        allHelix.append(i), allHelix.append(j)
    pk1Helix = max(set(allHelix), key=allHelix.count)
    pk2Helix = filter(lambda x: x != pk1Helix, allHelix)
    
    # construct list of base pairs
    pk1 = helicies[pk1Helix]
    pk2 = []
    for i in pk2Helix:
        for j in helicies[i]:
            pk2.append(j)
    
    return pk1, pk2

if __name__ == '__main__':
    import sys
    
    rna = CT(sys.argv[1])
    
    pk1, pk2 = extractPK(rna)
    print pk1
    print pk2
    
    
