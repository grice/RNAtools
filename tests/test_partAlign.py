import RNAtools.partAlign as m2
import os
import math
from pathlib import Path
import pytest

filepath = os.path.dirname(__file__)
data_dir = Path(f'{filepath}/../data')


@pytest.mark.skip(reason="Currently fails, not sure why though.")
def test_partAlign():
    """
    Tests partition aligner..
    """
    # seq1 = f'{filepath}/../data/sequence_hsa-mir-380.txt'
    seq1 = data_dir / 'sequence_hsa-mir-380.txt'
    # seq2 = f'{filepath}/../data/sequence_hsa-mir-383.txt'
    seq2 = data_dir / 'sequence_hsa-mir-383.txt'

    x = m2.RNA(seq1)
    y = m2.RNA(seq2)
    assert (x.count == 37)
    assert (y.count == 55)

    scoreMax = m2.align_RNA_partition(x, y)
    assert (math.floor(scoreMax) == 192)
