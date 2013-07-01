
from RNAtools import CT

if __name__ == '__main__':
    import sys
    from itertools import combinations
    a = sys.argv[1]
    rna = CT(a)
    
    pairs = combinations(range(1,len(rna.ct)+1),2)
    #pairs = range(1,len(rna.ct)+1)
    dist = []
    for i,j in pairs:
        count = rna.contactDistance(i,j)
        #print i,j,count
        dist.append(count)
        #print '*'*20
    
    import matplotlib.pyplot as plt
    
    plt.hist(dist,bins=20)
    plt.show()
