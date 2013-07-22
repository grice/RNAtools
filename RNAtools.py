########################################################
#                                                      #
#      RNAtools.py CT and dotplot data structures      #
#               v 0.7.0     22 July 2013               #
#                                                      #
#          author gregg rice - gmr@unc.edu             #
#                                                      #
########################################################

import sys

class CT:
    def __init__(self,fIN=None):
        #if givin an input file construct the ct object automatically
        if fIN:
            self.name = fIN
            self.num,self.seq,self.ct = self.readCT(fIN)
    
    def __str__(self):
        a = '{ Name= %s, len(CT)= %s }' % (self.name, str(len(self.ct)))
        return a
    
    def readCT(self,z):
        #reads a ct file, !reqires header!
        num,seq,bp = [],[],[]
        linesToRead = int(open(z).readlines()[0].rstrip().split()[0])
        #print linesToRead
        for i in open(z).readlines()[1:linesToRead+1]:
            a = i.rstrip().split()
            num.append(int(a[0])),seq.append(str(a[1])),bp.append(int(a[4]))
        return num,seq,bp
    
    def writeCT(self,fOUT):
        #writes a ct file
        w = open(fOUT,'w')
        line = '{0:6d} {1}\n'.format(len(self.num),self.name)
        for i in range(len(self.num)-1):
        #line+='{0:5d} {1} {2:5d} {3:5d} {4:5d} {5:5d} {0:5d}\n'\
        #.format(self.num[i],self.seq[i],self.num[i]-1,self.num[i]+1,self.ct[i])
            line += '{0:5d} {1} {2:5d} {3:5d} {4:5d} {0:5d}\n'.format(self.num[i],self.seq[i],self.num[i]-1,self.num[i]+1,self.ct[i])
        
        #last line is different
        i = len(self.num)-1
        line += '{0:5d} {1} {2:5d} {3:5d} {4:5d} {0:5d}\n'.format(self.num[i],self.seq[i],self.num[i]-1,0,self.ct[i])
        w.write(line)
        w.close()
    
    def copy(self):
        """
        returns a copy of the ct object
        """
        out = CT()
        out.name = self.name[:]
        out.num = self.num[:]
        out.seq = self.seq[:]
        out.ct = self.ct[:]
        return out
    
    def pairList(self):
        #returns a list of base pairs i<j
        out = []
        for nt in range(len(self.ct)):
            if self.ct[nt] != 0 and self.ct[nt] > self.num[nt]:
                out.append((self.num[nt],self.ct[nt]))
        return out
    
    def pair2CT(self, pairs, seq, name=None, skipConflicting=False):
        """
        pairs are an array of bp tuples ( i < j )
           i.e. [(4,26),(5,25),...]
        length is implied from the given sequence
        """
        length = len(seq)
        
        self.num = range(1,length+1)
        self.seq=seq
        
        #give it a name if it has one
        if name:
            self.name=name
        else:
            self.name='RNA_'+str(length)
        
        self.ct = []
        for i in range(length):self.ct.append(0)
        
        for i,j in pairs:
            if self.ct[i-1]!=0:
                print 'Warning: conflicting pairs, (%s - %s) : (%s - %s)' % (str(i),str(j),str(self.ct[i-1]),str(i))
                if skipConflicting: continue
            if self.ct[j-1]!=0:
                print 'Warning: conflicting pairs, (%s - %s) : (%s - %s)' % (str(i),str(j),str(j),str(self.ct[j-1]))
                if skipConflicting: continue
            self.ct[i-1]=j
            self.ct[j-1]=i
    
    def cutCT(self,start,end):
        """
        inclusively cuts a ct file and removes pairs outside the cutsite
        """
        start = int(start)
        end = int(end)
        
        out = CT()
        out.seq = self.seq[start-1:end]
        out.num = range(1,end-start+2)
        out.name = self.name + '_cut_'+str(start)+'_'+str(end)
        
        out.ct = []
        temp  = self.ct[start-1:end]
        # renumber from 1
        for nt in temp:
            nt_out = nt-(start-1)
            # cut out pairings that lie outside the window
            if nt_out <= 0 or nt_out > end-(start-1):
                nt_out = 0
            out.ct.append(nt_out)
        
        return out
    
    def contactDistance(self,i,j):
        """
        calculates the contact distance between pairs i,j in
        in the RNA using the RNAtools CT object. Method for
        calculating the contact is described in Hajdin et al
        (2013).
        """    
        #correct for indexing offset
        i, j = i-1, j-1
        
        #error out if nucleotide out of range
        if max(i,j) > len(self.ct):
            print 'Error!, nucleotide {0} out of range!'.format(max(i,j)+1)
            return
        
        #i must always be less than j, correct for this
        if j < i:
            i, j = j, i
        
        
        count = 0.0
        tempBack = 10000.0
        k = int(i)
        
        def backTrace(rna,j,k):
            bcount = 2
            k -= 1
            if k-1 == j:
                return bcount
            bcount += 1
            while k > j:
                if rna.ct[k] == 0:
                    k -= 1
                    bcount += 1
                # if not single stranded, exit the backtrace
                else:
                    return None
            return bcount
        # search forward through sequence
        while k < j:
            #debuging stuff, prevent infinite loop
            #if count > 200:
            #    print i,j
            #    break
            
            #nonpaired nucleotides are single beads
            if self.ct[k] == 0:
                k += 1
                #count += 6.5
                count += 1
            
            #branches through helices can't be skipped
            #treated as single beads
            elif self.ct[k] > j + 1:
                # try backtracing a few if it is close (within 10)
                if self.ct[k] - 10 < j:
                    back = backTrace(self, j,self.ct[k])
                    # if the backtracing is able to reach
                    # your nt, add its length to the count
                    # and break
                    if back:
                        #print 'backitup'
                        # store the backtracing count for later
                        # if it ends up being lower than the final
                        # we will use it instead
                        if count + back < tempBack:
                            tempBack = count + back
                k += 1
                #count += 6.5
                count += 1
            
            # simple stepping
            elif self.ct[k] < i + 1:
                k += 1
                #count += 6.5
                count += 1
            
            #handle branching, jumping the base of a helix is only 1
            else:
                # one count for the step to the next
                count += 1
                k = self.ct[k]
                
                # if by jumping you land on your nt
                # stop counting
                if k - 1 == j:
                    break
                
                # one count for jumping the helix
                #count += 15.0
                count += 1
            #print i,k,j
        finalCount = min(count, tempBack)
        return finalCount
    
    def extractHelicies(self,fillPairs=True):
        """
        returns a list of helicies in a given CT file and the nucleotides
        making up the list as a dict object. e.g. {0:[(1,50),(2,49),(3,48)}
        
        defaults to filling 1,1 and 2,2 mismatches
        """
        # first step is to find all the helicies
        rna = self.copy()
        if fillPairs:
            rna = rna.fillPairs()
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
                stillHelix = True
                
                previous = rna.ct[nt]
                
                while stillHelix:
                    # see if the helix juts up against another one
                    if abs(rna.ct[nt] - previous ) > 1:
                        break
                    #add the pairs
                    elif rna.ct[nt] != 0:
                        tempPairs.append((nt+1,rna.ct[nt]))
                        previous = rna.ct[nt]
                        nt+=1
                    #handle slip case
                    else:
                        if rna.ct[nt+1] == rna.ct[nt-1]-1:
                            nt+=1
                            #print 'slip'
                        else:
                            break
                
                # remove single bp helicies
                if len(tempPairs) <= 1:
                    continue
                helicies[heNum] = tempPairs
                heNum += 1
        return helicies
    
    def fillPairs(self):
        """
        fills 1,1 and 2,2 mismatches in an RNA structure
        """
        rna = self.copy()
        # fill in 1,1 mismatch, 2,2 mismatch
        for i in xrange(len(rna.ct)-3):
            if rna.ct[i+1] == 0:
                if rna.ct[i] - rna.ct[i+2] == 2:
                    rna.ct[i+1] = rna.ct[i] - 1
            if rna.ct[i+1] + rna.ct[i+2] == 0:
                if rna.ct[i] - rna.ct[i+3] == 3:
                    rna.ct[i+1] = rna.ct[i] - 1
                    rna.ct[i+2] = rna.ct[i] - 2
        return rna
    
    
    def extractPK(self,fillPairs=True):
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
        
        
        # make a copy so we don't destroy the original object
        rna = self.copy()
        
        
        #self.writeCT('foo.ct')
        #rna.writeCT('bar.ct')
        
        # get the helicies by calling the extract helix function
        helicies = rna.extractHelicies(fillPairs=fillPairs)
        heNum = len(helicies)

        
        # do the helicies overlap? Check for a crosshatching pattern
        # append them to a list if they have it.
        overlaps = [] # stores the helix number
        
        for i in xrange(heNum):
            for j in xrange(i+1,heNum):
                if checkOverlap(helicies[i],helicies[j]):
                    overlaps.append((i,j))
        
        # if there are no overlapping bps, return none
        if len(overlaps) == 0:
            return None, None
        
        
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


def padCT(targetCT, referenceCT,giveAlignment=False):
    """Aligns the target CT to the reference CT and pads the referece
    CT file with 000s in order to input into CircleCompare"""
    out = CT()
    out.seq = referenceCT.seq
    out.num = referenceCT.num
    
    #align target to reference
    seed = 200
    if len(targetCT.seq) <= seed:
        seed = len(targetCT.seq) - 1
    pos = 0
    maxScore = 0
    #print len(referenceCT.seq)
    for i in range(len(referenceCT.seq)-seed):
        a, b = referenceCT.seq[i:i+seed], targetCT.seq[:seed]
        s = 0
        # s = # of identical nts across the alignment
        for k, l in zip(a, b):
            if k == l:
                s += 1
        if s == seed:
            pos = i
            maxScore += 1
    # handle the exception when target and reference do not match
    if maxScore != 1:
        print 'reference and target do not match <EXIT>'
        sys.exit()
    
    #create the renumbered ct to fit within the reference
    ct = []
    for i in range(len(referenceCT.seq)):
        # if the target falls within the range of the reference ct
        #     then change the numbers
        # else pad the files with 000's
        if i >= pos and i < pos+len(targetCT.seq):
            val = targetCT.ct[i-pos]
            if val > 0:
                val += pos
            ct.append(val)
        else:
            ct.append(0)
    
    # set to the new ct file and return it
    out.ct = ct
    out.name = targetCT.name+'_renum_'+str(pos)
    if giveAlignment:
        return out, pos
    else:
        return out
    
 


class dotPlot:
    def __init__(self,fIN=None):
        #if givin an input file construct the dotplot object automatically
        self.name = None
        self.length = None
        self.dp = {}
        for elem in ['i','j','logBP']:
            self.dp[elem] = []        
        if fIN:
            self.name = fIN
            self.dp,self.length = self.readDP(fIN)
    
    def __str__(self):
        a = '{ Name= %s, len(RNA)= %s, entries(dotPlot)= %s }' % (self.name, str(self.length),str(len(self.dp['i'])))
        return a
    
    def readDP(self,fIN):
        out = dotPlot()

        out.length = int(open(fIN).readlines()[0].lstrip().split()[0])
        for n in open(fIN).readlines()[2:]:
            line = n.rstrip().split()
            out.dp['i'].append(int(line[0]))
            out.dp['j'].append(int(line[1]))
            out.dp['logBP'].append(float(line[2]))
        return out.dp,out.length
    
    def pairList(self):
        # returns a list of base pairs i< j from the dotplot
        out = []
        for n in range(len(self.dp['i'])):
            out.append((self.dp['i'][n],self.dp['j'][n]))
        return out
    
    def requireProb(self,logBP):
        #require probability at least as large as cutoff
        
        logBP = float(logBP)
        
        out = dotPlot()
        out.length = self.length
        out.name = self.name
        
        dp = self.dp
        
        for n in range(len(dp['i'])):
            if dp['logBP'][n] <= logBP:
                out.dp['i'].append(dp['i'][n])
                out.dp['j'].append(dp['j'][n])
                out.dp['logBP'].append(dp['logBP'][n])
        return out
    
    def trimEnds(self, trimSize, which='Both'):
        out = dotPlot()
        out.length = self.length
        out.name = self.name
        
        dp = self.dp
        
        if which == '5prime':
            for n in range(len(dp['i'])):
                if dp['i'][n] >= trimSize and \
                            dp['j'][n] >= trimSize:
                    out.dp['i'].append(dp['i'][n])
                    out.dp['j'].append(dp['j'][n])
                    out.dp['logBP'].append(dp['logBP'][n])
            return out
        
        if which =='3prime':
            for n in range(len(dp['i'])):
                if (self.length-trimSize) >= dp['i'][n] and \
                            (self.length-trimSize) >= dp['j'][n] :
                    out.dp['i'].append(dp['i'][n])
                    out.dp['j'].append(dp['j'][n])
                    out.dp['logBP'].append(dp['logBP'][n])
            return out
        
        for n in range(len(dp['i'])):
            if (self.length-trimSize) >= dp['i'][n] >= trimSize and \
                        (self.length-trimSize) >= dp['j'][n] >= trimSize:
                out.dp['i'].append(dp['i'][n])
                out.dp['j'].append(dp['j'][n])
                out.dp['logBP'].append(dp['logBP'][n])
        
        return out
