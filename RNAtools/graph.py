"""
Tools to convert RNA 2 structures to a graph representation.
"""

import itertools as it
from collections import Counter
from itertools import product

import numpy as np
import pandas as pd
from sklearn.preprocessing import LabelBinarizer

import networkx as nx
import RNAtools

lb = LabelBinarizer()
letters = ['A', 'U', 'G', 'C']  # fixed. Update as necessary.
edge_kinds = ['basepair', 'backbone']  # fixed. Update as necessary.


class CTGraph(object):
    """
    CTGraph object definition.

    Accepts a .ct file as input, and returns a CTGraph object that houses the
    graph and .ct file data (as a pandas DataFrame).
    """

    def __init__(self, ctfile, data=None, name='', file=None, **attr):
        # Read in the CT file as a pandas DataFrame.
        self.ctfile = ctfile
        edge_df = pd.read_fwf(ctfile, skiprows=1, header=None)
        edge_df.columns = ['position', 'letter', 'pos-1', 'pos+1', 'pair_pos',
                           'position_repeated']
        edge_df['letter'] = \
            edge_df['letter'].apply(lambda x: x.upper().replace("T", "U"))
        edge_df['kind'] = 'basepair'
        del edge_df['position_repeated']

        # Construct graph
        self.graph = nx.from_pandas_edgelist(edge_df,
                                             source='position',
                                             target='pair_pos',
                                             edge_attr='kind')
        # Annotate pairing partners on nodes.
        for n1, n2, d in self.graph.edges(data=True):
            self.graph.node[n1]['pairing_partner'] = n2
            self.graph.node[n2]['pairing_partner'] = n1

        # Add in backbone edges
        for n1, n2 in zip(sorted(self.graph.nodes()),
                          sorted(self.graph.nodes())[1:]):
            self.graph.add_edge(n1, n2, kind='backbone')

        # Add in node metadata
        for r, d in edge_df.iterrows():
            self.graph.node[d['position']]['letter'] = d['letter']

        # Remove the unnecessary node zero.
        self.graph.remove_node(0)

        # Annotate the graph with vectorized features.
        self.annotate()

        self._edge_df = edge_df

        # miRNA name
        self.mirna_name = ctfile.split('/')[-1].split('_')[0]

    def annotate(self):
        """
        Convenience function to annotate nodes and edges.
        """
        self.annotate_nodes()
        self.annotate_edges()

    def annotate_edges(self):
        """
        Annotate the graph edges with features.

        Due to the NetworkX API, it is much more convenient to annotate with
        key-value pairs such that each value is just a non-vector value.

        Thus, rather than having:

        - 'kind': [0, 1]

        We have:

        - 'is_basepair': 1
        - 'is_backbone': 0
        """
        lb.fit(edge_kinds)
        for n1, n2, d in self.graph.edges(data=True):
            if d['kind'] == 'basepair':
                self.graph.edges[n1, n2]['is_basepair'] = 1
                self.graph.edges[n1, n2]['is_backbone'] = 0
            elif d['kind'] == 'backbone':
                self.graph.edges[n1, n2]['is_backbone'] = 1
                self.graph.edges[n1, n2]['is_basepair'] = 0
            self.graph.edges[n1, n2]['kind_enc'] = \
                lb.transform([d['kind']]).ravel()

    def annotate_nodes(self):
        """
        Annotate the graph nodes with vectorized features.
        """
        lb.fit(letters)
        for n, d in self.graph.nodes(data=True):
            self.graph.node[n]['letter_enc'] = \
                lb.transform([d['letter']]).ravel()
            self.graph.node[n]['has_basepair'] = has_basepair(self.graph, n)
            # Initialize bp_prob node feature. Note that this is an
            # initialization, not the actual probability. This must be
            # annotated after-the-fact by calling
            # self.annotate_basepair_probabilities().
            self.graph.node[n]['bp_prob'] = np.array(0)

            # Initialize the fractional position to the termini. This metric is
            # zero-indexed, therefore, the first nucleotide will have a score
            # of 0, while the last nucleotide will have the score (n-1) / n.
            self.graph.node[n]['5p_frac_dist'] = np.array(n / (len(self.graph) + 1))  # noqa: E501
            self.graph.node[n]['3p_frac_dist'] = np.array(1 - self.graph.node[n]['5p_frac_dist'])  # noqa: E501

        # Annotate neighbors identity. This can only be done once all nodes are
        # annotated with `letter_enc`, therefore a second for-loop is
        # necessary.
        for n, d in self.graph.nodes(data=True):
            self.annotate_neighbors_identity(n, d)

    def annotate_neighbors_identity(self, n, d):
        """
        Annotate the neighbors identity. We specifically encode the following:

        - `fp_enc`: 5' nucleotide
        - `tp_enc`: 3' nucleotide
        - `bp_enc`: base paired nucleotide

        Each of these keys' values are numpy arrays.
        """
        # We take advantage of the fact that nodes are 5'->3' ordered by index.
        # Therefore, 5' node is always one less than the given nucleotide n,
        # and 3' node is always one greater than the given nucleotide n.
        fpnode = n - 1  # 5' node.
        tpnode = n + 1  # 3' node.
        bpnode = get_bp_partner(self.graph, n)  # basepair node

        # We check that the node exists in the graph. This is a boundary
        # condition check: the first nucleotide should not have any nucleotides
        # 5' to it, and the last nucleotide should not have any nucleotides
        # 3' to it.
        if not self.graph.has_node(fpnode):
            fpnode = None
        if not self.graph.has_node(tpnode):
            tpnode = None

        # Now, we get the encoding vector of the neighbors.
        if fpnode:
            fpvect = self.graph.node[fpnode]['letter_enc']
        else:
            fpvect = np.array([0, 0, 0, 0])

        if tpnode:
            tpvect = self.graph.node[tpnode]['letter_enc']
        else:
            tpvect = np.array([0, 0, 0, 0])

        if bpnode:
            bpvect = self.graph.node[bpnode]['letter_enc']
        else:
            bpvect = np.array([0, 0, 0, 0])

        # Finally, we encode the data on the graph using the keys `fp_enc`,
        # `tp_enc`, and `bp_enc`.
        self.graph.node[n]['fp_enc'] = fpvect
        self.graph.node[n]['tp_enc'] = tpvect
        self.graph.node[n]['bp_enc'] = bpvect

    def annotate_extrinsic_node_data(self, data):
        """
        Annotates nodes with extrinsic data (i.e. stuff that can't be computed
        directly off the graph data.)

        :param data: A dictionary of key-value pairs. The keys are the string
            label that goes on each node, and the values are a sequence (list,
            tuple, or ndarray) of values that are attached to each node,
            ordered from the 5' position to the 3' position.
        """
        for key, values in data[self.mirna_name].items():
            for n, v in zip(sorted(self.graph.nodes), values):
                if isinstance(v, int) or isinstance(v, float):
                    v = np.array([v])
                self.graph.node[n][key] = v

    def annotate_basepair_probabilities(self, mat):
        """
        Annotates the base pair matrix on the edges and on the nodes.

        :param mat: The base pairing matrix.
        """
        annotate_weights(self, mat)


def has_basepair(G, n):
    has_bp = False
    for u, v, d in G.edges(n, data=True):
        if d['kind'] == 'basepair':
            has_bp = True
            break
    return np.array([has_bp])


class HelixGraph(object):
    """
    Constructs a weighted directed graph of helical connectivity along the
    secondary structure. An edge is drawn between two helices that are
    connected through a junction, the direction of connectivity is 5'->3'
    and weight is the length of the bulge.

    :param ctfile: Path to CT file.
    :param extrinsic_feats: A dictionary of extrinsic features to add to the
        nodes.
    """

    def __init__(self, ctfile, extrinsic_feats):

        # Initialize ctfile object
        self.ct = RNAtools.CT(ctfile)

        # Initialize ctgraph
        self.ctgraph = CTGraph(ctfile)
        self.ctgraph.annotate_extrinsic_node_data(extrinsic_feats)
        # Construct Helix dictionary. Keys are helix IDs, values are 2-tuples
        # of (nt_pos1, nt_pos2).
        self.hlxDict = self.ct.extractHelices()
        # self.hlxDict = hlxDict
        # Swap helix dictionary keys and values, such that keys are nucleotide
        # positions (nt_pos1, nt_pos2 ...) and values are helix indices.
        self.positions = range(
            max(
                list(
                    it.chain.from_iterable(
                        it.chain.from_iterable(self.hlxDict.values())
                    )
                )
            )
        )

        # Helix has the nucleotide positions as keys and helix index as values.
        self.invHlxDict = {}

        # Initialize directed graph and add helix indices as nodes.
        self.graph = nx.DiGraph()
        # self.graph.add_nodes_from(self.hlxDict.keys())
        self.cycles = {}
        self.juncTypes = []
        self.juncTypes_node = np.ones(shape=(len(self.positions), 3)) * -1
        self.scaffold_helices = []
        self.annotate_helices()
        self.generate_helix_graph()
        self.annotate_nodes()
        self.annotate_edges()
        # self.getJuncType()

    def __str__(self):
        """
        Print method for the class
        """

        output_string = 'CT helices:\n'
        for hlx, bp in self.hlxDict.items():
            output_string += "{0}\t:\t{1}\n".format(hlx, str(bp))

        for edge in list(self.graph.edges()):

            output_string += "{0} -> {1}; ".format(edge[0], edge[1])

        output_string += "\n"
        output_string += "self.cycles: " + str(self.cycles) + "\n"
        output_string += "self.scaffold_helices: " + str(self.scaffold_helices) + "\n"  # noqa: E501
        output_string += "self.juncTypes: " + str(self.juncTypes) + "\n"
        output_string += str(self.juncTypes_node)
        return output_string

    def annotate_helices(self):
        """
        Annotates the CTGraph nodes with helix IDs.

        Credit to Varun Shivashankar for making this happen!!!!
        """
        h = 0  # helix ID
        b = 0  # bulge ID
        visited_nodes = list()
        g = self.ctgraph.graph
        prev_fp_node = 0
        prev_tp_node = len(g.nodes()) + 1
        prev_in_helix = False
        for fp_node in sorted(g.nodes()):
            if fp_node not in visited_nodes:
                tp_node = g.node[fp_node]['pairing_partner']
                if tp_node:
                    # The only time we increment a helix is if we were
                    # previously not in a helix, or if we skipped a node on the
                    # 3' end.
                    if (not prev_in_helix) or (tp_node != prev_tp_node - 1):
                        h += 1

                    g.node[fp_node]['helix_id'] = h
                    g.node[fp_node]['bulge_id'] = None
                    g.node[tp_node]['helix_id'] = h
                    g.node[tp_node]['bulge_id'] = None
                    visited_nodes.append(tp_node)
                    prev_fp_node = fp_node  # noqa: F841
                    prev_tp_node = tp_node
                else:
                    # The only time we increment a bulge is if we were
                    # previously in a helix.
                    if prev_in_helix:
                        b += 1
                    g.node[fp_node]['helix_id'] = None
                    g.node[fp_node]['bulge_id'] = b

            prev_in_helix = g.node[fp_node]['has_basepair']
            visited_nodes.append(fp_node)

    def generate_helix_graph(self):
        """
        The function performs a 5' to 3' walk on the CTGraph.graph object and
        annotates the nodes and edges.

        Add those pairs that belong to different helices, as edges.

        The weight of the graph is the length of the junction.
        """
        self.graph.add_node('start', nucleotides=[])
        self.graph.add_node('end', nucleotides=[])
        prev_helix_id = 'start'
        edge_nucleotides = []
        prev_in_helix = False
        for n, d in sorted(self.ctgraph.graph.nodes(data=True)):
            helix_id = d['helix_id']
            # If we are in a helix...
            if helix_id is not None:

                # Add helix ID to graph if it doesn't already exist.
                if not self.graph.has_node(helix_id):
                    self.graph.add_node(helix_id)
                    # We make sure that the helix graph node is annotated with
                    # a 'nucleotides' list.
                    self.graph.node[helix_id]['nucleotides'] = []
                # If the graph already has the helix ID as a node, then append
                # the nucleotide position to its 'nucleotides' list.
                if self.graph.has_node(helix_id):
                    self.graph.node[helix_id]['nucleotides'].append(n)

                # Check to see if we are switching states or not. If we
                # previously was not in a helix, then since we are now in a
                # helix, it is time to add an edge between the previous helix
                # and the current helix. This captures the "self-loop" state.
                # This is the (not prev_in_helix) conditional.
                #
                # The other criteria for adding an edge is if we are moving
                # from one helix to a new one. (i.e. prev_helix_id != helix_id)
                if not prev_in_helix or prev_helix_id != helix_id:
                    self.graph.add_edge(prev_helix_id, helix_id, nucleotides=edge_nucleotides)  # noqa: E501
                    edge_nucleotides = []
                    prev_helix_id = helix_id
                prev_in_helix = True

            # If we are not in a helix, then we append the current node to the
            # list of nucleotides spanning an edge. We also make sure to reset
            # the prev_in_helix state to False.
            else:
                edge_nucleotides.append(n)
                prev_in_helix = False

        # Finally, add in the "end" node.
        self.graph.add_edge(prev_helix_id, 'end', nucleotides=edge_nucleotides)

    def annotate_nodes(self):
        """
        Annotate helix graph nodes with features, in the same style as the
        CTGraph object.

        Features that are annotated are:

        - `num_nodes`: Number of nucleotides in the helix
        - `base_count`: The percentage of each base in the helix
        - `bp_percentage`: The percentage of each base pair type in the helix
        """
        # Annotate node-intrinsic metadata.
        for helix, data in self.graph.nodes(data=True):
            # Annotate number of nodes.
            self.annotate_num_node_nucleotides(helix, data)

            # Annotate the percentage of each nucleotide.
            annotate_base_count(self, helix, data)

            # Annotate base pairs
            self.annotate_basepairs(helix, data)

            # Annotate the counts of each type of base pair.
            self.annoate_basepair_composition(helix, data)

            # TODO: Annotate stacking energy of helix.

            # TODO: Annotate sum of shannon entropy

            # TODO: Annotate sum of base pair probabilities
            self.annotate_shannon(helix, data)

            # Add slots on the graph to annotate the node's edge summary stats.
            self.graph.node[helix]['num_in_edges'] = np.array(0)
            self.graph.node[helix]['num_out_edges'] = np.array(0)
            self.graph.node[helix]['num_self_loops'] = np.array(0)

            self.graph.node[helix]['sum_in_edge_lengths'] = np.array(0)
            self.graph.node[helix]['sum_out_edge_lengths'] = np.array(0)
            self.graph.node[helix]['sum_self_loops_lengths'] = np.array(0)
            # Base composition is stored in this order: ['A', 'U', 'G', 'C']
            self.graph.node[helix]['base_comp_in_edges'] = np.array([0, 0, 0, 0])  # noqa: E501
            self.graph.node[helix]['base_comp_out_edges'] = np.array([0, 0, 0, 0])  # noqa: E501
            self.graph.node[helix]['base_comp_self_loops'] = np.array([0, 0, 0, 0])  # noqa: E501

        # Annotate node's edge summary statistics.
        for h1, h2, data in self.graph.edges(data=True):
            # Deal with self-loops separately.
            if h1 == h2:
                bases = [self.ctgraph.graph.node[n]['letter'] for n in data['nucleotides']]  # noqa: E501
                bases = Counter(bases)

                self.graph.node[h1]['num_self_loops'] += 1
                self.graph.node[h1]['sum_self_loops_lengths'] += len(data['nucleotides'])  # noqa: E501
                self.graph.node[h1]['base_comp_self_loops'] += np.array([bases['A'], bases['U'], bases['C'], bases['G']])  # noqa: E501

            # Deal with regular edges.
            else:  # here, h1 != h2, and edge goes from h1 -> h2
                # Grab base composition first.
                bases = [self.ctgraph.graph.node[n]['letter'] for n in data['nucleotides']]  # noqa: E501
                bases = Counter(bases)

                # Annotate in-edge (for node h2)
                self.graph.node[h2]['num_in_edges'] += 1
                self.graph.node[h2]['sum_in_edge_lengths'] += len(data['nucleotides'])  # noqa: E501
                self.graph.node[h2]['base_comp_in_edges'] += np.array([bases['A'], bases['U'], bases['C'], bases['G']])  # noqa:

                # Annotate out-edges (for node h1)
                self.graph.node[h1]['num_out_edges'] += 1
                self.graph.node[h1]['sum_out_edge_lengths'] += len(data['nucleotides'])   # noqa: E501
                self.graph.node[h1]['base_comp_out_edges'] += np.array([bases['A'], bases['U'], bases['C'], bases['G']])  # noqa:

    def annotate_shannon(self, helix, data):
        """
        Helper method that annotates the sum over all shannon entropies of
        nucleotides in the helix and the shannon entropy normalized by number
        of nodes.
        """
        sum_shannon = np.array(0.0)
        for nuc in data['nucleotides']:
            shannon = self.ctgraph.graph.node[nuc]['shannon']
            sum_shannon += shannon[0]
        self.graph.node[helix]['sum_shannon'] = sum_shannon

        avg_shannon = np.array(sum_shannon / len(data['nucleotides']))
        if np.isnan(avg_shannon):
            avg_shannon = np.array(0)
        self.graph.node[helix]['avg_shannon'] = avg_shannon

    def annotate_num_node_nucleotides(self, helix, data):
        """
        Helper method that annotates the number of nodes in a helix.

        :param helix: A node in the HelixGraph.
        :param data: The HelixGraph node's attribute dictionary.
        """
        self.graph.node[helix]['num_nodes'] = np.array(len(data['nucleotides']))  # noqa: E501

    def annotate_basepairs(self, helix, data):
        """
        Helper method that annotates helix graph nodes with basepairs
        contained within them.

        Attribute data structure is:

            basepairs: [(nt_pos, pair_nt_pos)]

        Note, we annotate the paired nucleotide position, and NOT the
        nucleotide identity. That has to be computed from the underlying
        `self.CTGraph.graph` data.

        We also do not record a single base pair two ways. If we record
        (nt_pos, pair_nt_pos), then we do not record (pair_nt_pos, nt_pos).

        Pairing information is taken from the CTGraph node attribute
        'pairing_partner'.

        :param helix: A node in the HelixGraph.
        :param data: The HelixGraph node's attribute dictionary.
        """
        # Record visited nucleotides and basepairs per helix node
        visited = []
        basepairs = []
        for nuc in data['nucleotides']:
            # Check if each nucleotide in the 'nucleotides' attr
            # has been visited
            if nuc not in visited:
                # append nuc and CTGraph 'pairing_partner' to basepairs
                # mark nuc as visited
                pp = self.ctgraph.graph.nodes[nuc]['pairing_partner']
                basepairs.append(tuple([nuc, pp]))
                visited.append(nuc)
                visited.append(pp)
        # Annotate helix graph nodes with attr basepairs
        self.graph.node[helix]['basepairs'] = basepairs

    def annoate_basepair_composition(self, helix, data):
        """
        Annotates the base pair composition.

        We read the 'basepairs' attribute of each node, which gives us the
        nucleotide positions of base pairs. We then go back into the CTGraph
        object and read the exact nucleotide identity from there.
        """
        nt_pairings = [f'{l1}{l2}' for l1, l2 in product('AUGC', repeat=2)]  # AA, AU, AG, AC, ...  # noqa: E501
        bp_count = np.zeros(len(nt_pairings))
        for (pos1, pos2) in data['basepairs']:
            nuc1 = self.ctgraph.graph.node[pos1]['letter']
            nuc2 = self.ctgraph.graph.node[pos2]['letter']

            # Through the next three lines of code, "AU" is treated the same
            # as "UA".
            bp = f'{nuc1}{nuc2}'
            if bp not in nt_pairings:
                bp = f'{nuc2}{nuc1}'
            # print(bp_count[nt_pairings.index(bp)])

            bp_count[nt_pairings.index(bp)] += 1
        self.graph.node[helix]['bp_count'] = bp_count

    def annotate_edges(self):
        """
        Annotate edges with metadata.
        """

        for u, v, d in self.graph.edges(data=True):
            # Annotate how many nucleotides exist on the edge.
            self.annotate_num_edge_nucleotides(u, v, d)
            # Annotate how many backbone bonds exist on the edge. This is a
            # measure of the path length between two nodes.
            self.annotate_num_edge_backbones(u, v, d)
            # Annotate edge weight. This is taken to be 1 / number of backbone
            # bonds.
            self.annotate_edge_weight(u, v, d)

    def annotate_num_edge_nucleotides(self, u, v, d):
        """
        Helper method that simply annotates on each edge how many nucleotides
        exist on the edge.
        """
        self.graph.edges[u, v]['num_nucleotides'] = len(d['nucleotides'])

    def annotate_num_edge_backbones(self, u, v, d):
        """
        Annotates the number of nodes associated with an edge.
        """
        self.graph.edges[u, v]['num_backbone_bonds'] = len(d['nucleotides']) + 1  # noqa: E501

    def annotate_edge_weight(self, u, v, d):
        """
        Annotates the edge weights. Defined as 1 / number of edges.
        """
        self.graph.edges[u, v]['weight'] = 1 / (len(d['nucleotides']) + 1)


class HelixBulgeGraph(object):
    """
    Helix Bulge Graph.

    This graph is defined as follows:
    - Nodes are one of helices or unpaired regions (bulges).
    - Edges are unweighted and undirected, and connect helices to bulges and
      bulges to helices.

    To construct the HelixBulgeGraph, the easiest way is to leverage the
    HelixGraph. See the `generate_graph()` class method for details.
    """

    def __init__(self, ctfile, extrinsic_feats):
        super(HelixBulgeGraph, self).__init__()
        self.helix_graph = HelixGraph(ctfile, extrinsic_feats=extrinsic_feats)
        self.graph = nx.Graph()

        self.generate_graph()
        self.annotate_nodes()

    def generate_graph(self):
        """
        Leverages the HelixGraph to generate the HelixBulgeGraph. The key here
        is to iterate over all edges, and add a "bulge" node.
        """

        for i, (h1, h2, d) in enumerate(self.helix_graph.graph.edges(data=True)):  # noqa: E501
            # It is easier to just add in all the nodes first, and then remove
            # the 'start' and 'end' dummy nodes at the end.
            self.graph.add_node(f'H{h1}',
                                nucleotides=self.helix_graph.graph.node[h2]['nucleotides'],  # noqa: E501
                                bipartite='helix')
            self.graph.add_node(f'H{h2}',
                                nucleotides=self.helix_graph.graph.node[h2]['nucleotides'],  # noqa: E501
                                bipartite='helix')

            # Add in the bulge nodes. These are what were the edges in the
            # helix_graph.
            self.graph.add_node(f'B{i}',
                                nucleotides=d['nucleotides'],
                                bipartite='bulge')

            self.graph.add_edge(f'H{h1}', f'B{i}')
            self.graph.add_edge(f'H{h2}', f'B{i}')

        self.graph.remove_node(f'Hstart')
        self.graph.remove_node(f'Hend')

    def annotate_nodes(self):
        """
        Annotate nodes.
        """
        for n, d in self.graph.nodes(data=True):
            annotate_base_count(self, n, d)
            self.annotate_partition(n, d)
            self.annotate_size(n, d)
            self.annotate_partition_composition(n, d)

    def annotate_partition(self, node, data):
        """
        Annotate node partitions. It is named "partitions" because this graph
        is a bipartite graph, in which one partition is the "helix" partition,
        and the other is the "bulge" partition.
        """
        self.graph.node[node]['is_helix'] = \
            np.array(int(data['bipartite'] == 'helix'))  # noqa: E501

    def annotate_size(self, node, data):
        """
        Annotate the size of the node. This is the number of nucleotides
        captured in the nodes.
        """
        if data['bipartite'] == 'helix':
            self.graph.node[node]['helix_size'] = np.array(len(data['nucleotides']))  # noqa: E501
            self.graph.node[node]['bulge_size'] = np.array(0)
        elif data['bipartite'] == 'bulge':
            self.graph.node[node]['helix_size'] = np.array(0)
            self.graph.node[node]['bulge_size'] = np.array(len(data['nucleotides']))  # noqa: E501

    def annotate_partition_composition(self, node, data):
        """
        Takes advantage of the base percentages already annotated to annotate
        the partition base composition.
        """
        if data['bipartite'] == 'helix':
            self.graph.node[node]['helix_base_count'] = data['base_count']
            self.graph.node[node]['bulge_base_count'] = np.zeros(4)
        elif data['bipartite'] == 'bulge':
            self.graph.node[node]['bulge_base_count'] = data['base_count']
            self.graph.node[node]['helix_base_count'] = np.zeros(4)


def annotate_base_count(obj, helix, data):
    """
    Helper method that annotates the percentage of each nucleotide on
    each node.

    :param obj: Should be one of HelixGraph, HelixBulgeGraph.
    :param helix: A node in the HelixGraph.
    :param data: The HelixGraph node's attribute dictionary.
    """
    if isinstance(obj, HelixGraph):
        nts = [obj.ctgraph.graph.node[n]['letter']
               for n
               in data['nucleotides']
               ]
    elif isinstance(obj, HelixBulgeGraph):
        nts = [obj.helix_graph.ctgraph.graph.node[n]['letter']
               for n
               in data['nucleotides']
               ]
    nt_counts = Counter(nts)
    # total_nts = sum(nt_counts.values())
    base_count = []
    for letter in 'AUGC':
        try:
            base_count.append(nt_counts[letter])
        except ZeroDivisionError:
            base_count.append(0)
    obj.graph.node[helix]['base_count'] = np.array(base_count)


def get_bp_partner(g, n):
    """
    Returns the base pairing partner of n.

    If node is not base paired, returns None. Else, returns other node.
    """
    if not has_basepair(g, n):
        return None
    else:
        for u, v, d in g.edges(n, data=True):
            if d['kind'] == 'basepair':
                if u == n:
                    return v
                else:
                    return u


def annotate_weights(ct, bp_matrix):
    """
    Annotate weights on the graph based on base pairing probability.

    Both the backbone and base pairing edges get weights.

    Note: the backbones automatically get a weight of 1, while the base
    pairing edges get a weight that corresponds to their base pairing
    probability.

    In order to not interfere with graph layout algorithms, which
    automatically respect the 'weight' parameter, the probability (i.e.
    weight) is stored with the `prob_edge` key.

    In addition to storing the probability on the edge, we also store the edge
    probability on the node. We store under the key bp_prob as follows:

        0.0 < p <= 1.0 for probability,
        0.0 for unpaired.

    :param ct: A CTGraph object.
    :param bp_matrix: A square matrix. Values are base pairing probabilities.
        Note: The CTGraph is 1-indexed, but the bp matrix is zero-indexed.
    :returns: A CTGraph object with weights annotated.
    """
    for u, v, d in ct.graph.edges(data=True):
        if d["kind"] == "backbone":
            w = 1
        elif d["kind"] == "basepair":
            # Annotate weight from basepair matrix.
            # we use u-1 and v-1 because matrix is zero-indexed but graph is
            # 1-indexed.
            w = bp_matrix[u - 1, v - 1]

            # Annotate the node with the basepair probability.
            ct.graph.node[u]['bp_prob'] = np.array(w)
            ct.graph.node[v]['bp_prob'] = np.array(w)
        ct.graph.edges[u, v]["prob_edge"] = w
