import sys
import numpy as np
import re


class CT:

    def __init__(self, fIN=None, suboptimal=False):
        """
        if givin an input file .ct construct the ct object automatically

        subptimal flag allows access to multiple structures in the same ct flag
        """
        if fIN:
            self.readCT(fIN, suboptimal=suboptimal)

    def __str__(self):
        """
        overide the default print statement for the object
        """
        a = '{ Name= %s, len(CT)= %s }' % (self.name, str(len(self.ct)))
        return a

    def readCT(self, fIN, suboptimal=False):
        """
        #reads a ct file, !reqires header!
        """
        num, seq, bp = [], [], []
        try:
            linesToRead = int(open(fIN).readlines()[0].rstrip().split()[0])
        except:
            print("file is not in .ct format. Requires Header")
            return 0

        # Gets the free energy
        energyMatch = re.compile(r'.*ENERGY = ([0-9\.\-]+).*')
        header = open(fIN).readlines()[0].strip()

        if energyMatch.match(header):
            self.energy = float(energyMatch.match(header).group(1))
        else:
            self.energy = np.nan

        for i in open(fIN).readlines()[1:linesToRead + 1]:
            a = i.rstrip().split()
            num.append(int(a[0])), seq.append(str(a[1])), bp.append(int(a[4]))

        # print linesToRead
        self.name = fIN
        self.num = num
        self.seq = seq
        self.ct = bp
        self.substructure = []
        # print(self.energy)

        substructures = []
        if suboptimal:

            # define the loop reset
            st_reset = linesToRead + 1
            with open(fIN) as f:
                for i, line in enumerate(f):
                    # header line of a new structure
                    if i % st_reset == 0:
                        conformer = []
                        curr_name = line
                    else:
                        a = line.rstrip().split()
                        conformer.append(int(a[4]))

                    # last line of a structure, append the structure
                    if i % st_reset == linesToRead:
                        s = CT()
                        s.name = str(curr_name)
                        s.num = self.num
                        s.seq = self.seq
                        s.ct = conformer
                        substructures.append(s)
            self.substructure = np.array(substructures)

    def writeCT(self, fOUT):
        """
        writes a ct file from the ct object
        """

        # handle empty ct object case
        if not self.ct:
            print("empty ct object. Nothing to write")
            return

        w = open(fOUT, 'w')
        line = '{0:6d} {1}\n'.format(len(self.num), self.name)
        for i in range(len(self.num) - 1):
            line += '{0:5d} {1} {2:5d} {3:5d} {4:5d} {0:5d}\n'.format(
                self.num[i], self.seq[i], self.num[i] - 1, self.num[i] + 1, self.ct[i]   # noqa: E501
            )

        # last line is different
        i = len(self.num) - 1
        line += '{0:5d} {1} {2:5d} {3:5d} {4:5d} {0:5d}\n'.format(
            self.num[i], self.seq[i], self.num[i] - 1, 0, self.ct[i]
        )  # noqa
        w.write(line)
        w.close()

    def copy(self):
        """
        returns a deep copy of the ct object
        """
        out = CT()
        out.name = self.name[:]
        out.num = self.num[:]
        out.seq = self.seq[:]
        out.ct = self.ct[:]
        return out

    def pairList(self):
        """
        #returns a list of base pairs i<j as a array of tuples:
        [(19,50),(20,49)....]
        """
        out = []
        for nt in range(len(self.ct)):
            if self.ct[nt] != 0 and self.ct[nt] > self.num[nt]:
                out.append((self.num[nt], self.ct[nt]))
        return out

    def dot2ct(self, sequence, dotstring, name=None):
        """
        function to convert dotbracket to a CT object

        pseudknots are indicated with distinct brackets
           e.g. (( [[ <<
        """

        pairs = []

        for i, char_i in enumerate(dotstring):
            # go through sequence until we have a
            # left base pair
            lmatch = "([{<"
            rmatch = ")]}>"

            # print i, char_i

            if char_i in lmatch:
                count = 1

                # set the chr type to match
                pos = lmatch.find(char_i)

                # set the chr that increases the counter
                adder = lmatch[pos]
                # set the chr the closes the basepair
                search = rmatch[pos]

                # go through the rest of the sequence
                for j, char_j in enumerate(dotstring[i + 1:]):
                    if char_j == adder:
                        count += 1
                    elif char_j == search:
                        count -= 1

                    # we have a match append the pair
                    if count == 0:
                        # pairs are 1-indexed
                        # right pair is (i+1) + (j+1)
                        pairs.append((i + 1, i + j + 2))
                        break

        # out = CT()
        self.pair2CT(pairs, sequence, name=name)

    def pair2CT(self, pairs, seq, name=None, skipConflicting=True):
        """
        constructs a ct object from a list of base pairs and a sequence

        pairs are an array of bp tuples ( i < j )
           i.e. [(4,26),(5,25),...]
        length is implied from the given sequence
        """
        length = len(seq)

        self.num = list(range(1, length + 1))
        self.seq = seq

        # give it a name if it has one
        if name:
            self.name = name
        else:
            self.name = 'RNA_' + str(length)

        self.ct = []
        for i in range(length):
            self.ct.append(0)

        for i, j in pairs:
            if self.ct[i - 1] != 0:
                print(
                    'Warning: conflicting pairs, ({0} - {1}) : ({2} - {3})'.format(
                        str(i), str(j), str(self.ct[i - 1]), str(i)
                    )
                )  # noqa
                if skipConflicting:
                    continue

            if self.ct[j - 1] != 0:
                print(
                    'Warning: conflicting pairs, ({0} - {1}) : ({2} - {3})'.format(
                        str(i), str(j), str(j), str(self.ct[j - 1])
                    )
                )  # noqa
                if skipConflicting:
                    continue

            self.ct[i - 1] = j
            self.ct[j - 1] = i

    def cleanCT(self):
        """
        removes non-canonical base pairs and removes 1-bp helicies
        from the RNA

        any pair that is not AU GC or GU will be removed.
        """
        # copy the CT object
        out = self.copy()

        # pull the number pairs, and the sequence pairs
        # example of data structure:
        #     pairs = [(1,50)]
        #     seq_pairs = ["A", "U"]

        pairs = out.pairList()
        seq_pairs = []

        for i, j in pairs:
            k, l = out.seq[i - 1], out.seq[j - 1]
            seq_pairs.append((k, l))

        allowedPairs = ["AU", "UA", "GC", "CG", "GU", "UG", "AT", "TA", "GT", "TG"]

        # go through the pairlist and check that pair is in the
        # allowed list, append to a new array
        out_pairs = []
        for k, l in zip(pairs, seq_pairs):
            kind = "".join(l)
            if kind in allowedPairs:
                out_pairs.append(k)

        out_ct = CT()
        out_ct.pair2CT(out_pairs, out.seq, name=out.name)

        # go through cleaned list and remove 1 bp long helicies
        # use the .bp representation of the pairs
        out_pairs_clean = []
        for i, j in enumerate(out_ct.ct):
            # if there is no pair ignore this position
            if j == 0:
                continue

            # skip redunantly annotated pairs
            if (i + 1) > j:
                continue

            # catch the 5p end of the RNA
            if i == 0:
                # if the end is a 1 bp stack the
                # next pair will be 0, don't add it to the
                # cleaned pairs
                if out_ct.ct[i + 1] == 0:
                    continue

            # catch the 3p end of the rna
            if i + 1 == len(out_ct.ct):
                prev_right = out_ct.ct[i - 1]
                next_right = out_ct.ct[i + 1]
                if prev_right == 0 and next_right == 0:
                    continue

            # catch the general case for an internal 1-bp helix

            # fix 1-index
            across = out_ct.ct[i] - 1

            prev_left = out_ct.ct[across - 1]

            # edge case, if nt is paired to the last nt
            # there is no next_left, catch and set to non-zero
            if across + 1 >= len(out_ct.ct):
                next_left = 1
            else:
                next_left = out_ct.ct[across + 1]

            prev_right = out_ct.ct[i - 1]
            next_right = out_ct.ct[i + 1]

            # print("#"*30)
            # print(i, prev_left, prev_right)
            # print(i, next_left, next_right)
            if prev_left == 0 and next_left == 0:
                continue

            elif prev_right == 0 and next_right == 0:
                continue

            # passing all the checks add it to the new ct file
            else:
                out_pairs_clean.append((i + 1, out_ct.ct[i]))

        out_pairs_clean_ct = CT()
        out_pairs_clean_ct.pair2CT(out_pairs_clean, out.seq, name=out.name)
        return out_pairs_clean_ct

    def cutCT(self, start, end):
        """
        inclusively cuts a ct file and removes pairs outside the cutsite
        """
        start = int(start)
        end = int(end)

        out = CT()
        out.seq = self.seq[start - 1:end]
        out.num = list(range(1, end - start + 2))
        out.name = self.name + '_cut_' + str(start) + '_' + str(end)

        out.ct = []
        temp = self.ct[start - 1:end]
        # renumber from 1
        for nt in temp:
            nt_out = nt - (start - 1)
            # cut out pairings that lie outside the window
            if nt_out <= 0 or nt_out > end - (start - 1):
                nt_out = 0
            out.ct.append(nt_out)

        return out

    def stripCT(self):
        """
        returns an array the length of the ct object
        in a non-redundant form - e.g. only gives pairs i<j
        """
        pairs = self.pairList()
        halfPlexCT = np.zeros_like(self.ct)
        for i, j in pairs:
            halfPlexCT[i - 1] = j
        return halfPlexCT

    def contactDistance(self, i, j):
        """
        calculates the contact distance between pairs i,j in
        in the RNA using the RNAtools CT object. Method for
        calculating the contact is described in Hajdin et al
        (2013).
        """
        # correct for indexing offset
        i, j = i - 1, j - 1

        # error out if nucleotide out of range
        if max(i, j) > len(self.ct):
            print('Error!, nucleotide {0} out of range!'.format(max(i, j) + 1))
            return

        # i must always be less than j, correct for this
        if j < i:
            i, j = j, i

        count = 0.0
        tempBack = 10000.0
        k = int(i)

        def backTrace(rna, j, k):
            bcount = 2
            k -= 1
            if k - 1 == j:
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
            # debuging stuff, prevent infinite loop
            # if count > 200:
            #    print i,j
            #    break

            # nonpaired nucleotides are single beads
            if self.ct[k] == 0:
                k += 1
                # count += 6.5
                count += 1

            # branches through helices can't be skipped
            # treated as single beads
            elif self.ct[k] > j + 1:
                # try backtracing a few if it is close (within 10)
                if self.ct[k] - 10 < j:
                    back = backTrace(self, j, self.ct[k])
                    # if the backtracing is able to reach
                    # your nt, add its length to the count
                    # and break
                    if back:
                        # print 'backitup'
                        # store the backtracing count for later
                        # if it ends up being lower than the final
                        # we will use it instead
                        if count + back < tempBack:
                            tempBack = count + back
                k += 1
                # count += 6.5
                count += 1

            # simple stepping
            elif self.ct[k] < i + 1:
                k += 1
                # count += 6.5
                count += 1

            elif self.ct[k] < k + 1:
                k += 1
                count += 1
            # print "Backtracking prevented, going forward to "+str(k+1)

            # handle branching, jumping the base of a helix is only 1
            else:
                # one count for the step to the next
                count += 1
                k = self.ct[k] - 1

                # if by jumping you land on your nt
                # stop counting
                if k - 1 == j:
                    break

        # one count for jumping the helix
        # count += 15.0
        # count += 1
        # print i,k,j
        finalCount = min(count, tempBack)
        return finalCount

    def extractHelices(self, fillPairs=True):
        """
        returns a list of helices in a given CT file and the nucleotides
        making up the list as a dict object. e.g. {0:[(1,50),(2,49),(3,48)}

        defaults to filling 1,1 and 2,2 mismatches
        """
        # first step is to find all the helices
        rna = self.copy()
        if fillPairs:
            rna = rna.fillPairs()
        helices = {}
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
                    if abs(rna.ct[nt] - previous) > 1:
                        break

                    # add the pairs
                    elif rna.ct[nt] != 0:
                        tempPairs.append((nt + 1, rna.ct[nt]))
                        previous = rna.ct[nt]
                        nt += 1
                    # handle slip case
                    else:
                        if rna.ct[nt + 1] == rna.ct[nt - 1] - 1:
                            nt += 1
                        # print 'slip'
                        else:
                            break

                # remove single bp helices
                if len(tempPairs) <= 1:
                    continue

                helices[heNum] = tempPairs
                heNum += 1
        return helices

    def fillPairs(self):
        """
        fills 1,1 and 2,2 mismatches in an RNA structure
        """
        rna = self.copy()
        # fill in 1,1 mismatch, 2,2 mismatch
        for i in range(len(rna.ct) - 3):
            if rna.ct[i + 1] == 0:
                if rna.ct[i] - rna.ct[i + 2] == 2:
                    rna.ct[i + 1] = rna.ct[i] - 1
            if rna.ct[i + 1] + rna.ct[i + 2] == 0:
                if rna.ct[i] - rna.ct[i + 3] == 3:
                    rna.ct[i + 1] = rna.ct[i] - 1
                    rna.ct[i + 2] = rna.ct[i] - 2
        return rna

    def analyzeSubstructures(self):
        """
        analyzes and calculates the fraction of each
        suboptimal structure in the population
        """
        import md5

        if len(self.substructure) == 0:
            print("ERROR: No substructures present")
            return 0

        xcount = {}
        xstruct = {}
        for structure in self.substructure:
            hashed = md5.md5(structure).hexdigest()
            # print hashed
            # print structure
            if hashed not in xcount:
                xcount[hashed] = 1
                xstruct[hashed] = structure
            else:
                xcount[hashed] += 1
        return xcount, xstruct

    def extractPK(self, fillPairs=True):
        """
        returns the pk1 and pk2 pairs from a CT file. Ignores single base
        pairs. PK1 is defined as the helix crossing the most other helices.
        if there is a tie, the most 5' helix called pk1

        returns pk1 and pk2 as a list of base pairs e.g [(1,10),(2,9)...
        """

        def checkOverlap(h1, h2):
            # only need to check one set of pairs from each
            # of the helices. Test is to see if they form
            # a cross hatching pattern
            if max(h1[0]) > min(h2[0]) > min(h1[0]):
                if max(h2[0]) > max(h1[0]):
                    return True

            if max(h2[0]) > min(h1[0]) > min(h2[0]):
                if max(h1[0]) > max(h2[0]):
                    return True

            else:
                return False

        # make a copy so we don't destroy the original object
        rna = self.copy()

        # self.writeCT('foo.ct')
        # rna.writeCT('bar.ct')

        # get the helices by calling the extract helix function
        helices = rna.extractHelices(fillPairs=fillPairs)
        heNum = len(helices)

        # do the helices overlap? Check for a crosshatching pattern
        # append them to a list if they have it.
        overlaps = []  # stores the helix number

        for i in range(heNum):
            for j in range(i + 1, heNum):
                if checkOverlap(helices[i], helices[j]):
                    overlaps.append((i, j))

        # print overlaps
        # print '#'*30
        # if there are no overlapping bps, return none
        if len(overlaps) == 0:
            return None, None

        # determine which is pk1
        allHelix = []
        for i, j in overlaps:
            allHelix.append(i), allHelix.append(j)
        pk1Helix = max(set(allHelix), key=allHelix.count)
        pk2Helix = [x for x in allHelix if x != pk1Helix]

        # construct list of base pairs
        pk1 = helices[pk1Helix]
        pk2 = []
        for i in pk2Helix:
            for j in helices[i]:
                pk2.append(j)

        return pk1, pk2

    def padCT(self, referenceCT, giveAlignment=False):
        """
        utilizes the global padCT method on this object
        """
        return padCT(self, referenceCT, giveAlignment)

    def readSHAPE(self, fIN):
        """
        utilizes the global readSHAPE method and appends a the data to the
        object as .shape
        """
        self.shape = readSHAPE(fIN)
        if len(self.shape) < len(self.ct):
            print("warning! shape array is smaller than the CT range")

    def writeSHAPE(self, fOUT):
        """
        utilizes the global writeSHAPE method, and writes the .shape array
        attached to the object to a file
        """
        try:
            writeSHAPE(self.shape, fOUT)
        except:
            print("No SHAPE data present")
            return


def padCT(targetCT, referenceCT, giveAlignment=False):
    """Aligns the target CT to the reference CT and pads the referece
    CT file with 000s in order to input into CircleCompare"""
    out = CT()
    out.seq = referenceCT.seq
    out.num = referenceCT.num

    # align target to reference
    seed = 200
    if len(targetCT.seq) <= seed:
        seed = len(targetCT.seq) - 1
    pos = 0
    maxScore = 0
    # print len(referenceCT.seq)
    for i in range(len(referenceCT.seq) - seed):
        a, b = referenceCT.seq[i:i + seed], targetCT.seq[:seed]
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
        print('reference and target do not match <EXIT>')
        sys.exit()

    # create the renumbered ct to fit within the reference
    ct = []
    for i in range(len(referenceCT.seq)):
        # if the target falls within the range of the reference ct
        #     then change the numbers
        # else pad the files with 000's
        if i >= pos and i < pos + len(targetCT.seq):
            val = targetCT.ct[i - pos]
            if val > 0:
                val += pos
            ct.append(val)
        else:
            ct.append(0)

    # set to the new ct file and return it
    out.ct = ct
    out.name = targetCT.name + '_renum_' + str(pos)
    if giveAlignment:
        return out, pos

    else:
        return out


def readSHAPE(fIN):
    """
    reads an RNA structure .shape or .map file. Returns an array of the SHAPE
    data
    """
    shape = []
    for i in open(fIN, "rU").readlines():
        line = i.rstrip().split()[1]
        shape.append(float(line))
    return shape


def writeSHAPE(shape, fOUT):
    """
    writes the data from a shape array into the file fOUT
    """
    w = open(fOUT, "w")

    for i in range(len(shape)):
        line = "{0}\t{1}\n".format(i + 1, shape[i])
        w.write(line)
    w.close()


def readSeq(fIN, type='RNAstructure'):
    """
    reads an RNAstructure sequence file format and converts it to an
    array of nucleotides. e.g.: ['A','G','C','C'...]

    also returns the name of the sequence from the file
    """

    # strip the input file of comments
    seqRaw = []
    for i in open(fIN, "rU").read().split():
        if len(i) == 0:
            continue

        if i[0] == ";":
            continue

        seqRaw.append(i)

    name = seqRaw[0]
    seqJoin = ''.join(seqRaw[1:])
    seq = []
    for i in seqJoin:
        if i == '1':
            break

        seq.append(i)
    return name


class dotPlot:

    def __init__(self, fIN=None):
        # if givin an input file construct the dotplot object automatically
        self.name = None
        self.length = None
        self.dp = {}
        for elem in ['i', 'j', 'logBP']:
            self.dp[elem] = np.array([])
        if fIN:
            self.name = fIN
            self.dp, self.length = self.readDP(fIN)

    def __str__(self):
        a = '{ Name= %s, len(RNA)= %s, entries(dotPlot)= %s }' % (
            self.name, str(self.length), str(len(self.dp['i']))
        )
        return a

    def readDP(self, fIN):
        out = dotPlot()

        i = []
        j = []
        logBP = []

        out.length = int(open(fIN).readlines()[0].lstrip().split()[0])
        for n in open(fIN).readlines()[2:]:
            line = n.rstrip().split()
            # out.dp['i'] = np.append(out.dp['i'],int(line[0]))
            # out.dp['j'] = np.append(out.dp['j'],int(line[1]))
            # out.dp['logBP'] = np.append(out.dp['logBP'],float(line[2]))
            i.append(int(line[0]))
            j.append(int(line[1]))
            logBP.append(float(line[2]))
        # ln+=1
        # if ln%500 == 0:
        #    print ln
        # load in as python object first, MUCH faster
        out.dp['i'] = np.append(out.dp['i'], i)
        out.dp['j'] = np.append(out.dp['j'], j)
        out.dp['logBP'] = np.append(out.dp['logBP'], logBP)
        return out.dp, out.length

    def writeDP(self, fOUT):
        """
        writes a DP file back to disk
        """
        w = open(fOUT, "w")

        # file header
        w.write("{0}\ni\tj\t-log10(Probability)\n".format(self.length))

        # resort the array first by j
        jsort = np.argsort(self.dp['j'], kind='mergesort')
        i, j, logbp = self.dp['i'][jsort], self.dp['j'][jsort], self.dp['logBP'][
            jsort
        ]  # noqa

        # then by i
        isort = np.argsort(self.dp['i'], kind='mergesort')
        i, j, logbp = self.dp['i'][isort], self.dp['j'][isort], self.dp['logBP'][
            isort
        ]  # noqa

        self.dp['i'] = i
        self.dp['j'] = j
        self.dp['logBP'] = logbp

        # main file
        for n in range(len(self.dp['i'])):
            line = "{0}\t{1}\t{2}\n".format(
                int(self.dp['i'][n]), int(self.dp['j'][n]), self.dp['logBP'][n]
            )  # noqa
            w.write(line)
        # close file
        w.close()

    def pairList(self):
        # returns a list of base pairs i< j from the dotplot
        out = []
        for n in range(len(self.dp['i'])):
            out.append((int(self.dp['i'][n]), int(self.dp['j'][n])))
        return out

    def requireProb(self, minlogBP, maxlogBP=0.0):
        """
        require probability at least as large as cutoff

        setting a minlogBP of 3 will give a dp object with probibilites of at
        least 0.001
        """

        minlogBP = float(minlogBP)
        maxlogBP = float(maxlogBP)

        out = dotPlot()
        out.length = self.length
        out.name = self.name

        dp = self.dp

        # select which nts are between a certain cutoff
        probFilter = (dp['logBP'] < minlogBP) * (dp['logBP'] > maxlogBP)

        out.dp['logBP'] = dp['logBP'][probFilter]
        out.dp['i'] = dp['i'][probFilter]
        out.dp['j'] = dp['j'][probFilter]

        # for n in range(len(dp['i'])):
        #    if dp['logBP'][n] <= logBP:
        #        out.dp['i'].append(dp['i'][n])
        #        out.dp['j'].append(dp['j'][n])
        #        out.dp['logBP'].append(dp['logBP'][n])
        return out

    def trimEnds(self, trimSize, which='Both'):
        out = dotPlot()
        out.length = self.length
        out.name = self.name

        dp = self.dp

        if which == '5prime':
            dpfilter = (dp['i'] >= trimSize) * (dp['j'] >= trimSize)
            out.dp['i'] = np.array(dp['i'][dpfilter])
            out.dp['j'] = np.array(dp['j'][dpfilter])
            out.dp['logBP'] = np.array(dp['logBP'][dpfilter])
            # for n in range(len(dp['i'])):
            #    if dp['i'][n] >= trimSize and \
            #                dp['j'][n] >= trimSize:
            #        out.dp['i'].append(dp['i'][n])
            #        out.dp['j'].append(dp['j'][n])
            #        out.dp['logBP'].append(dp['logBP'][n])
            return out

        if which == '3prime':
            dpfilter = (self.length - trimSize >= dp['i']) * (
                self.length - trimSize >= dp['j']
            )  # noqa
            out.dp['i'] = np.array(dp['i'][dpfilter])
            out.dp['j'] = np.array(dp['j'][dpfilter])
            out.dp['logBP'] = np.array(dp['logBP'][dpfilter])
            # for n in range(len(dp['i'])):
            #    if (self.length-trimSize) >= dp['i'][n] and \
            #                (self.length-trimSize) >= dp['j'][n] :
            #        out.dp['i'].append(dp['i'][n])
            #        out.dp['j'].append(dp['j'][n])
            #        out.dp['logBP'].append(dp['logBP'][n])
            return out

        dpfilter1 = (dp['i'] >= trimSize) * (dp['j'] >= trimSize)
        dpfilter2 = (self.length - trimSize >= dp['i']) * (
            self.length - trimSize >= dp['j']
        )  # noqa
        dpfilter = dpfilter1 * dpfilter2
        out.dp['i'] = np.array(dp['i'][dpfilter])
        out.dp['j'] = np.array(dp['j'][dpfilter])
        out.dp['logBP'] = np.array(dp['logBP'][dpfilter])

        return out

    def calcShannon(self, printOut=False, toFile=None):
        dp = self.dp
        shannon = []

        if toFile:
            w = open(toFile, 'w')

        # precalculate nlog(n), array is already in -logForm
        dp['nlogn'] = dp['logBP'] * 10 ** (-1 * dp['logBP'])

        for nt in range(1, self.length + 1):
            # if col i or j is the nt include the value and sum it
            mask = ((nt == dp['i']) + (nt == dp['j']))
            summed = np.sum(dp['nlogn'][mask])

            # catch rounding errors:
            if summed < 0:
                summed = 0

            # pairing prob, breakpoint
            # x = 10**(-dp['logBP'][mask])
            # print nt, np.sum(x), x[x>0.001]

            # print to stout if desired
            if printOut:
                print(nt, summed)

            # write to file if desired
            if toFile:
                line = '\t'.join(map(str, [nt, summed]))
                w.write(line + '\n')
            shannon.append(summed)
        if toFile:
            w.close()
        return shannon

    def averageSlippedBPs(self, struct=None, predictedOnly=True):
        """
        replaces PlusandMinus script. If a helix in a predicted structure is
        slipped +/-1 nt we need to sum the predicted probabilities otherwise
        the predicted Shannon entropy will be artificially high.

        turning off predicted only will go through all i<j basepairs and merge
        them in preference of liklihood. This is more compuationally intensive
        """
        dotPlot = self.dp
        # dotPlotCopy = {'logBP':copy.deepcopy(dotPlot['logBP'])}

        # this is the value in -log10(prob), 2 = a prob of 0.01  # noqa
        slippedCutoff = 2
        slipped = []

        # if a reference structure is given, merge pairs to it first
        if struct:
            for pair in range(1, len(struct.ct) - 1):
                # define the base pairs
                pair_i = pair + 1
                pair_j = struct.ct[pair]

                # skip non-paired nucleotides
                if pair_j == 0:
                    continue

                # skip pairs i > j so we don't double count
                if pair_j < pair_i:
                    continue

                # create some search filters
                filter_i = dotPlot['i'] == pair_i
                filter_j = dotPlot['j'] == pair_j

                filter_i_before = dotPlot['i'] == pair_i - 1
                filter_i_after = dotPlot['i'] == pair_i + 1

                filter_j_before = dotPlot['j'] == pair_j - 1
                filter_j_after = dotPlot['j'] == pair_j + 1

                # find i,j union before, at, and after
                filterDP = {}
                filterDP['before_j'] = filter_j_before * filter_i
                filterDP['before_i'] = filter_j * filter_i_before

                filterDP['after_j'] = filter_j_after * filter_i
                filterDP['after_i'] = filter_j * filter_i_after

                # define current point
                at = filter_j * filter_i

                # handle slippage, first define base prob
                prob_at = 10 ** (-dotPlot['logBP'][at])

                # cycle through all filter combinations
                for filterPair in list(filterDP.keys()):

                    # shorthand variable for current filter
                    curr = filterDP[filterPair]

                    # if pair exists ...
                    if np.sum(curr) == 1:
                        # add it to predicted pair probability
                        prob_at += 10 ** (-dotPlot['logBP'][curr])

                        # add pair to slipped list if it meets slip criteria
                        if dotPlot['logBP'][curr] < slippedCutoff:
                            slipped.append((pair_i, pair_j))

                        # now set it to a very low probability for zero
                        dotPlot['logBP'][curr] = 50

                # return to -log10 format
                dotPlot['logBP'][at] = -np.log10(prob_at)

        if not predictedOnly:
            # go through all i<j basepair combinations and check to see if
            # there is a slipped base pair
            for i, j in self.pairList():

                # correct for python counting starting at 0
                pair_i, pair_j = int(i + 0), int(j + 0)

                # see if there exists a basepair for this combination
                filter_i = (dotPlot['i'] == pair_i)
                filter_j = (dotPlot['j'] == pair_j)

                filter_union = (filter_i * filter_j)

                # only continue if the pair is reasonably likely
                # (>1% chance of forming)
                if (np.sum(filter_union) == 1 and dotPlot['logBP'][filter_union] < 3):
                    filterList = {}

                    # define the various types of slippage
                    filter_ibefore = (dotPlot['i'] == pair_i - 1)
                    filter_jbefore = (dotPlot['j'] == pair_j - 1)
                    filter_iafter = (dotPlot['i'] == pair_i + 1)
                    filter_jafter = (dotPlot['j'] == pair_j + 1)

                    # index filters to a dict
                    filterList['before_i'] = filter_ibefore * filter_j
                    filterList['before_j'] = filter_i * filter_jbefore
                    filterList['after_i'] = filter_iafter * filter_j
                    filterList['after_j'] = filter_i * filter_jafter

                    # define the prob at in normal normal space
                    prob_at = 10 ** (-dotPlot['logBP'][filter_union])

                    # go through each of the filters
                    for pairFilter in list(filterList.keys()):
                        curr = filterList[pairFilter]

                        # if the current pair exists in the dotplot
                        if np.sum(curr) == 1:
                            # check to see if it's less likely than current
                            # pairs
                            if dotPlot['logBP'][filter_union] < dotPlot['logBP'][
                                curr
                            ]:  # noqa
                                # and add it to the current pair if it is
                                prob_at += 10 ** (-dotPlot['logBP'][curr])

                                # set to a very low probabliity afterwards
                                dotPlot['logBP'][curr] = 50
                    dotPlot['logBP'][filter_union] = -np.log10(prob_at)

        return slipped

    def pairingProb(self):
        """
        returns the pairing probability

        returned array is equal to the length of the RNA
        """

        dp = self.dp
        prob = []

        # convert the -log10 values to untransformed probability
        dp['prob'] = 10 ** (-1 * dp['logBP'])

        # use a numpy array to mask entires containing i or j and sum
        # the probability for each nt in the RNA
        for nt in range(1, self.length + 1):
            mask = ((nt == dp['i']) + (nt == dp['j']))
            summed = np.sum(dp['prob'][mask])

            prob.append(summed)

        return prob
