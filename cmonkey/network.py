"""network.py - cMonkey network module

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""


class NetworkEdge:
    """class to represent a network edge"""

    def __init__(self, source, target, score):
        """create an edge instance"""
        self.__source = source
        self.__target = target
        self.__score = score

    def source(self):
        """returns the source node"""
        return self.__source

    def target(self):
        """returns the target node"""
        return self.__target

    def score(self):
        """returns the edge score"""
        return self.__score

    def set_score(self, score):
        """sets a new score for this edge"""
        self.__score = score

    def __str__(self):
        """returns string representation"""
        return "%s -> %s w = %s" % (self.__source, self.__target,
                                    str(self.__score))

    def __repr__(self):
        """returns string representation"""
        return str(self)


class Network:
    """class to represent a network graph"""

    def __init__(self, name, edges):
        """creates a network from a list of edges"""
        self.__name = name
        self.__edges = edges

    def name(self):
        """returns the name of the network"""
        return self.__name

    def edges(self):
        """returns the list of edges"""
        return self.__edges

    def num_edges(self):
        """returns the number of edges in this graph"""
        return len(self.__edges)

    def total_score(self):
        """returns the sum of edge scores"""
        total = 0.0
        for edge in self.__edges:
            total += edge.score()
        return total

    def normalize_scores_to(self, score):
        """normalizes all edge scores so that they sum up to
        the specified score"""
        total = self.total_score()
        if score != total:
            # score_e / score_total * score == score_e * (score_total / score)
            # we use this to save a division per loop iteration
            scale = float(score) / float(total)
            for edge in self.__edges:
                edge.set_score(edge.score() * scale)

    def __repr__(self):
        return "Network: %s\n# edges: %d\n" % (self.__name,
                                               len(self.__edges))

    @classmethod
    def create(cls, name, edges):
        """standard Factory method"""
        added = {}
        network_edges = []

        def add_edge(edge):
            """adds an edge to network_edges"""
            key = "%s:%s" % (edge.source(), edge.target())
            added[key] = True
            network_edges.append(edge)

        for edge in edges:
            key = "%s:%s" % (edge.source(), edge.target())
            key_rev = "%s:%s" % (edge.target(), edge.source())
            if key not in added:
                add_edge(edge)
            if key_rev not in added:
                add_edge(NetworkEdge(edge.target(), edge.source(),
                                     edge.score()))
        return Network(name, network_edges)