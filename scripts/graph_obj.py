#!/usr/bin/env python
"""

Class of the directed graph produced by minigraph

Consist of three classess:

    Graph : collection of nodes and edges
    Node  : node with attributes
    Edge  : edge with attributes

"""

from dataclasses import dataclass
from collections import defaultdict


class Graph():
    """

    Store graph properties

    Attributes
    nodes : collection of nodes object
    edges: collection of edges object

    Methods
    add_node : add new node in the graph
    add_edge: add new edge in the graph
    graph_len: total pangenome length
    nonref_len: total non-ref sequences in the graph
    conv_adjlist: convert edges into an adjacency list representation



    """

    def __init__(self, nodes=None, edges=None):
        if nodes:
            self.nodes = nodes
        else:
            self.nodes = []

        if edges:
            self.edges = edges
        else:
            self.edges = []

    def add_nodes(self, node):
        """

        add a node to the graph

        Params:
        A node object 

        Returns: None, but will add node to the nodes

        """

        self.nodes.append(node)

    def add_edges(self, edge):
        """

        add an edge to graph

        Params: an edges object

        Returns: None, but will add edge to the edges
        """

        self.edges.append(edge)

    def graph_len(self):
        """

        Calculate the length of graph

        Params: None
        Output: length of graphs (bp)

        """
        totlen = 0
        for node in self.nodes:
            totlen += node.nodelen
        return totlen

    def nonref_len(self):
        """
        Calculate length of the non-ref sequence in the graph

        Params: None
        Output: Length non-ref sequences (bp)

        """
        nonref = 0
        for node in self.nodes:
            if node.noderank > 0:
                nonref += node.nodelen
        return nonref

    def conv_adjlist(self):
        """

        Convert edges representation into a adjacency list

        Params: None
        Output: Dict with parent as key and list containing all childs

        """

        adj_list = defaultdict(list)

        for edge in self.edges:
            adj_list[edge.parent].append(edge.child)


@dataclass
class Node():
    """

    Store node properties

    Attributes:

    nodeid : node identifier (s1,s2..)
    nodelen : length of sequences in the node (bp)
    nodeseq : nucleotide sequences in the node
    nodechr : chromosomal location node
    nodecoord : coordinate of the node
    nodecol : node label
    noderank : order of inclusion of node (0 reference)
    """

    nodeid: str
    nodelen: int
    nodeseq: str = ""
    nodechr: str
    nodecoord: int
    nodecol: str = ""
    noderank: int


@dataclass
class Edge():
    """

    Store edge properties
    Edges always directed

    parent ---> child

    Strand denotes orientation (+/-)

    """

    parent: str
    child: str
    strand_parent: str
    strand_child: str
