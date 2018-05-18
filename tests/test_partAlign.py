import RNAtools.partAlign as m2
import os, math


filepath = os.path.dirname(__file__)


def test_partAlign():
    """
    Tests partition aligner..
    """
    seq1 = f'{filepath}/../data/sequence_hsa-mir-380.txt'
    seq2 = f'{filepath}/../data/sequence_hsa-mir-383.txt'

    x = m2.RNA(seq1)
    y = m2.RNA(seq2)
    assert  (x.count == 37)
    assert  (y.count == 55)

    scoreMax = m2.align_RNA_partition(x,y)
    assert (math.floor(scoreMax) == 192)



#if __name__ == "__main__":
#    test_partAlign()

# results below from run_alignment_test function built in to partAlign.pyx
#('RNA1', 37)
#('RNA2', 55)
#('score:', 192.67144492775518)
###############################
#backtrace
#[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60]
#[6, 7, 8, 9, 10, 11, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 64, 65]
#192.67144492775518
#0.27054667472839355

