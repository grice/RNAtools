from RNAtools.graph import HelixGraph
import pytest
import networkx as nx
import os.path as osp
import matplotlib.pyplot as plt
import itertools as it


@pytest.fixture
def helix_graph():
    parent = osp.abspath('..')
    ct_path = osp.join(parent, 'data', 'hsa-mir-218-2_890.ct')
    # ct_path = '../../data/foo.ct'
    hg = HelixGraph(ct_path)
    hg.annotate_basepairs()
    return hg


@pytest.fixture
def helix_199():
    parent = osp.abspath('..')
    ct_path = osp.join(parent, 'data', 'hsa-mir-199b.ct')
    hg = HelixGraph(ct_path)
    return hg


@pytest.mark.development
def test_helix_graph(helix_graph):
    assert len(helix_graph.graph) == 10
    # Get nucleotides involved in basepairs attribute
    test_helix_position = 1
    test_helix = helix_graph.graph.nodes(data=True)[test_helix_position]
    print(test_helix)
    basepairs = test_helix['basepairs']
    nucleotides_in_basepairs = sorted(it.chain.from_iterable(basepairs))
    nucleotides = test_helix['nucleotides']
    assert nucleotides_in_basepairs == nucleotides


@pytest.mark.development
def test_helix_199(helix_199):
    print('hsa-mir-199b')
    print(helix_199.graph.edges())


if __name__ == '__main__':

    hg = helix_graph()
    print(hg.hlxDict)
    hg.annotate_helices()
    nx.draw_kamada_kawai(hg.graph, with_labels=True)
    plt.show()

    g = hg.ctgraph.graph
    labeldict = {n: f'{n} {d["helix_id"]}' for n, d in hg.ctgraph.graph.nodes(data=True)}
    print(hg.graph.edges(data=True))
    nx.draw_kamada_kawai(hg.ctgraph.graph, with_labels=True, labels=labeldict)
    plt.show()
