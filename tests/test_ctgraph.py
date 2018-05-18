from RNAtools.graph import CTGraph
import os


filepath = os.path.dirname(__file__)


def test_ct_graph():
    """
    Tests that CT graph construction is done correctly.
    """
    ct = CTGraph(f'{filepath}/../data/foo.ct')
    assert len(ct.graph.nodes) == 21
    assert len(ct.graph.edges) == 28

    print(ct)
