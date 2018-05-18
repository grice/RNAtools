import numpy as np
import itertools as it
import time, os
from multiprocessing import Process, Queue, Manager
#cimport numpy as np

# commands to speed up cython, code tracks
#cython: boundscheck=False
#cython: wraparound=False

cdef class relArr:
    # convenience obj to keep track of relative base pairing
    # built once per RNA and overwritten during calcuation
    cdef int [:] relpos
    cdef double [:] part

    def __init__(self, length):
        self.relpos, self.part = self.buildArr(length)
    def buildArr(self,length):
        cdef int [:] relpos = np.zeros(length,dtype=np.int32)
        cdef double [:] part = np.zeros(length,dtype=np.float64)
        return relpos,part

cdef class RNA:
    """
    cython class of RNAtools dotplot object
    """
    cdef int length, count
    cdef int [:] i
    cdef int [:] j
    cdef double [:] prob

    def __init__(self, fIN=None):

        if fIN:
            out = self.read_dp(fIN, 0.001)
            self.length = out[0]
            self.i      = out[1]
            self.j      = out[2]
            self.prob   = out[3]
            self.count  = out[4]

    def read_dp(self, fIN, thresh):
        cdef int length   = 0
        cdef int count    = 0
        cdef int arrinit  = 0

        # run through the file and get the
        # number of entries to init the array
        for n,line in enumerate(open(fIN)):
            if n>=2:
                spl = line.rstrip().split()
                pp = 10.0**(-float(spl[2]))
                if pp > thresh:
                    arrinit += 1

        cdef int [:] i_pair   = np.zeros(arrinit,dtype=np.int32)
        cdef int [:] j_pair   = np.zeros(arrinit,dtype=np.int32)
        cdef double [:] pairprob = np.zeros(arrinit,dtype=np.float64)

        for n,line in enumerate(open(fIN)):
            #print n
            if n == 0:
                length = int( line.rstrip() )
            if n >= 2:
                spl =  line.rstrip().split()
                pp = 10.0**(-float(spl[2]))
                if pp > thresh:
                    i_pair[count] = int(spl[0])
                    j_pair[count] = int(spl[1])
                    pairprob[count] = float(pp)
                    count += 1
        return length,i_pair,j_pair,pairprob,count

cdef double scoreDeltaPFS2(int i,int j, RNA t1, RNA t2, relArr t1_rel, relArr t2_rel):
    """
    this is the scoring function that calculates how closely the two structures
    are to each other. If the relative RNA structures are the same the score
    is high
    """
    #init some things
    cdef int k = 0
    cdef int l = 0
    cdef int count = 0
    cdef int pos_i, pos_j

    # sum_t1,t2 track how long each of the arr steps are
    cdef int sum_t1 = 0
    cdef int sum_t2 = 0
    cdef int relDist = 0

    cdef double score = 0.0

    # loop through RNA1 to calculate the relative positioning
    for k in range(t1.count):
        if abs(t1.i[k] - i) < 2 or abs(t1.j[k] - i) <2:
            sum_t1 += 1
            t1_rel.relpos[count] = t1.j[k] - t1.i[k]
            t1_rel.part[count] = t1.prob[k]

            count += 1

    # loop through RNA2, same code as above
    #   future me: wrap this into a function
    count = 0
    for k in range(t2.count):
        if abs(t2.i[k] - j) < 2 or abs(t2.j[k] - j) < 2:
            sum_t2 += 1
            t2_rel.relpos[count] = t2.j[k] - t2.i[k]
            t2_rel.part[count] = t2.prob[k]

            count += 1

    # scoring loop, strength is based on pairing probability between
    # nucleotides, ex: t1_rel.part[k] value
    for k in range(sum_t1):
        for l in range(sum_t2):
            relDist = abs(t1_rel.relpos[k] - t2_rel.relpos[l])
            # if structures match give the highest score
            if relDist < 1:
                score += 1.0 * (t1_rel.part[k] + t2_rel.part[l])

            # if structures are close give a partial score
            elif relDist < 3:
                score += 0.5 * (t1_rel.part[k] + t2_rel.part[l])
            #elif relDist < 5:
            #    score += 0.25 * (t1_rel.part[k] + t2_rel.part[l])
            # if structures don't match, give a penalty
            else:
                score += -0.25 * (t1_rel.part[k] + t2_rel.part[l])
    return score


cdef align_RNA_partition_backtrace(RNA t1, RNA t2):
    """
    Function to align two RNA partition functions using a modified
    smith-waterman algorthm. There is some function duplication between this
    function and align_RNA_partition.

    Input is two RNA objects to aligns
    Output are two 1-D arrays describing the trajectory of the backtrace for RNA1
    and RNA2 along with the maxalignment scoreMat
    """
    cdef double ins_del = -12.0
    cdef double local_max = 0.0
    cdef double global_max = 0.0
    cdef int i,j,
    cdef int start_i = 0
    cdef int start_j = 0
    cdef double [:,:] scoreMat = np.zeros((t1.length+1,t2.length+1),np.float64)

    t1_rel = relArr(t1.length)
    t2_rel = relArr(t2.length)

    # fill in the scoring matrix
    for i in range(1,t1.length+1):
        for j in range(1,t2.length+1):
            local_max = max(scoreMat[i-1,j-1]+scoreDeltaPFS2(i,j,t1,t2,t1_rel,t2_rel),
                            scoreMat[i-1,j] + ins_del,
                            scoreMat[i,j-1]+ins_del,
                            0)
            scoreMat[i,j] += local_max
            if local_max > global_max:
                global_max = float(local_max)
                start_i, start_j = i, j

    # backtrace implementation
    cdef int curr_i = int(start_i)
    cdef int curr_j = int(start_j)
    cdef int next_i = 0
    cdef int next_j = 0
    cdef int z = 0
    cdef int [:] path_i = np.zeros(t1.length+t2.length, np.int32)
    cdef int [:] path_j = np.zeros(t1.length+t2.length, np.int32)
    cdef double smax = 100.0
    cdef int count = 0

    # start position is found in the fill scoreing matrix routine
    #start_i, start_j = np.unravel_index(np.argmax(scoreMat), scoreMat.shape)

    path_i[count] = start_i
    path_j[count] = start_j

    while smax > 0.0:
        count += 1

        # find the next step in the path
        z = np.argmax([scoreMat[curr_i -1, curr_j -1],
                       scoreMat[curr_i - 1, curr_j],
                       scoreMat[curr_i, curr_j - 1]])

        # move the indicies in that direction
        if z == 0:
            next_i, next_j = curr_i - 1, curr_j - 1
        elif z==1:
            next_i, next_j = curr_i - 1, curr_j
        elif z==2:
            next_i, next_j = curr_i, curr_j - 1

        # store the path to the array
        path_i[count] = next_i
        path_j[count] = next_j

        # swap for the next step in the loop
        curr_i, curr_j = next_i, next_j

        # get the score
        smax = scoreMat[curr_i, curr_j]

    # init an array for reversing the backtrace
    if path_i[count] == 0 or path_j[count] == 0:
        count -= 1

    cdef int [:] path_i_r = np.zeros(count+1, np.int32)
    cdef int [:] path_j_r = np.zeros(count+1, np.int32)
    cdef int r_count = 0
    #print count

    # adjust for zero index, score matrix first row and column
    # set to 0 by definition
    #for i in range(count, -1, -1):
    #    #print i, r_count
    #    path_i_r[r_count] = path_i[i] - 1
    #    path_j_r[r_count] = path_j[i] - 1
    #    r_count += 1
    while count > -1:
        curr_i = path_i[count]
        curr_j = path_j[count]

        # code to catch some the meandering paths
        # if path ventures up, then left just call
        # as a match
        if count > 1:
            next_i = path_i[count-2]
            next_j = path_j[count-2]
            #print "##########"
            #print curr_i, next_i - 1
            #print curr_j, next_j - 1

            if next_i - 1  == curr_i and next_j - 1 == curr_j:
                path_i_r[r_count] = curr_i - 1
                path_j_r[r_count] = curr_j - 1
                count -= 2
                r_count += 1
                continue
        path_i_r[r_count] = curr_i - 1
        path_j_r[r_count] = curr_j - 1

        count -= 1
        r_count += 1

    return path_i_r[:r_count], path_j_r[:r_count], global_max

cdef align_RNA_partition(RNA t1, RNA t2):
    """
    function to calculate the max alignment score from the scoreing matrix
    input is two RNA objects as defined in this cython object
    output is a float of the max alignment score between the two structures
    """
    # indel penalty
    cdef double ins_del = -12.0

    # inits for local and global max variables
    cdef double local_max = 0.0
    cdef double global_max = 0.0
    cdef int i,j

    # init the scoring matrix
    cdef double [:,:] scoreMat = np.zeros((t1.length+1,t2.length+1),np.float64)

    # init the relArr objects --relative distance of pairScore
    # this is done once and the object is reused
    # to reduce the number of init function calls and a significant speedup
    # since creating new objects has a significant timing cost
    t1_rel = relArr(t1.length)
    t2_rel = relArr(t2.length)

    # fill in the scoring matrix
    # wiki has a surprisingly nice visualization of how this is done
    # https://en.wikipedia.org/wiki/Smith-Waterman_algorithm
    for i in range(1,t1.length+1):
        for j in range(1,t2.length+1):
            local_max = max(scoreMat[i-1,j-1]+scoreDeltaPFS2(i,j,t1,t2,t1_rel,t2_rel),
                            scoreMat[i-1,j] + ins_del,
                            scoreMat[i,j-1]+ins_del,
                            0)
            scoreMat[i,j] += local_max
            if local_max > global_max:
                global_max = float(local_max)

    return global_max

def thread_job(inQueue, resultQueue, rnaList):
    while True:
        # get the next set of RNAs to align
        work = inQueue.get()
        # if there arent any break out of the loop
        #print(work)
        if work[0] == None: break

        # set the RNAs
        rna1 = rnaList[work[0]]
        rna2 = rnaList[work[1]]

        # align and package the result
        score = align_RNA_partition(rna1, rna2)
        scoreTuple = (work, score)

        # store result in the queue
        resultQueue.put(scoreTuple)
    return

def align(leftDP, rightDP):
    """
    aligns two RNA .dp files

    output is the backtrace and the globalMax score
    """

    RNA_1 = RNA(leftDP)
    RNA_2 = RNA(rightDP)

    pair_i, pair_j, globalMax = align_RNA_partition_backtrace(RNA_1, RNA_2)

    return list(pair_i), list(pair_j), globalMax

def oneToManyAlign(masterDP, children):
    """
    masterDP is the DP file name for the structure to align to
    children is an array of file names

    output is a one dimensional array of length children
    """

    cdef double scoreMax = 0.0
    master = RNA(masterDP)
    child_scores = []
    for dp in children:
        y = RNA(dp)
        scoreMax = align_RNA_partition(master, y)
        child_scores.append(float(scoreMax))
    return child_scores

def pairwiseAlign(children, numProcs = 1):
    """
    runs a pairwise alignment for the array of alignment files
    calculation can be with multiple threads via the numProcs option

    returns a pairwise alignment matrix
    """

    cdef double scoreMax = 0.0

    childObj = []
    pairScore = np.zeros((len(children),len(children)),dtype=np.float64)

    # load structure data into memory
    for dp in children:
        childObj.append(RNA(dp))

    # init some of the multiprocessing stuff
    manager = Manager()
    work = manager.Queue(numProcs)
    result = Queue()

    # init the threads
    threadPool = []
    for i in range(numProcs):
        p = Process(target=thread_job, args=(work, result, childObj))
        p.start()
        threadPool.append(p)
    #print(threadPool[0])

    # supply the workers
    print("supply work")
    expectSize = 0
    for n, pair in enumerate(it.combinations(range(len(childObj)), 2)):
        if n % 1000 == 0: print("loopcount:", n)
        #print(pair)
        work.put(pair)
        expectSize += 1

    # close finished threads
    print("closing threads")
    for i in range(numProcs):
        work.put((None,))

    # wait until the work is picked up
    #while not work.empty(): time.sleep(1)
    while result.qsize() < expectSize: time.sleep(1)
    # get the results
    while True:
        #print(result.empty())
        if result.empty(): break
        pair, scoreMax = result.get()
        i,j = pair
        pairScore[i][j] = float(scoreMax)

    #for i,j in it.combinations(range(len(childObj)), 2):
    #    scoreMax = align_RNA_partition(childObj[i],childObj[j])
    #    pairScore[i][j] = float(scoreMax)

    return pairScore

def run_alignment_test():
    """
    simple alignment test between two miRNA objects
    ('RNA1', 37)
    ('RNA2', 55)
    ('score:', 192.67144492775518)
    """

    cdef double scoreMax = 0.0

    tstart = time.time()

    test1 = "sequence_hsa-mir-380.txt"
    test2 = "sequence_hsa-mir-383.txt"
    f = relArr(40)
    #x = read_dp(test1, 0.01)
    #y = read_dp(test2, 0.01)
    #align_rna(x,y)
    x = RNA(test1)
    print( "RNA1", x.count)
    y = RNA(test2)
    print("RNA2", y.count)
    #scoreDeltaPFS2(5,10,x,y)
    for i in range(101):
        scoreMax = align_RNA_partition(x,y)
    print("score:", scoreMax)

    print("#"*30)
    print("backtrace")
    tmp1, tmp2, scoreMax  = align_RNA_partition_backtrace(x,y)
    print(list(tmp1))
    print(list(tmp2))
    print scoreMax



    print time.time()-tstart

if __name__ == "__main__":
    run_alignment_test()
