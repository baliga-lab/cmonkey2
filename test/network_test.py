"""network_test.py - unit tests for network module

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import unittest
from network import NetworkEdge, Network


class NetworkEdgeTest(unittest.TestCase):  # pylint: disable-msg=R0904
    """Test class for NetworkEdge"""

    def test_create(self):
        """tests creating an edge"""
        edge = NetworkEdge('from', 'to', 3.14)
        self.assertEquals('from', edge.source())
        self.assertEquals('to', edge.target())
        self.assertEquals(3.14, edge.score())

    def test_set_score(self):
        """tests setting an edge score"""
        edge = NetworkEdge('from', 'to', 3.14)
        edge.set_score(2.5)
        self.assertEquals(2.5, edge.score())


class NetworkTest(unittest.TestCase):  # pylint: disable-msg=R0904
    """Test class for Network"""

    def test_construct(self):
        """tests creating a network using the constructor"""
        edge1 = NetworkEdge('n1', 'n2', 123)
        edge2 = NetworkEdge('n3', 'n2', 234)
        network = Network('network', [edge1, edge2])
        self.assertEquals('network', network.name())
        self.assertEquals(2, network.num_edges())
        self.assertTrue(edge1 in network.edges())
        self.assertTrue(edge2 in network.edges())
        self.assertEquals(357, network.total_score())

    def test_create_unique_edges(self):
        """tests creating a network using the standard factory method"""
        edge1 = NetworkEdge('n1', 'n2', 123)
        edge2 = NetworkEdge('n3', 'n2', 234)
        network = Network.create('network', [edge1, edge2])
        self.assertEquals('network', network.name())
        self.assertEquals(4, network.num_edges())
        self.assertEquals(714, network.total_score())

    def test_create_overlapping_edges(self):
        """tests creating a network using the standard factory method
        no duplicate edge will be generated"""
        edge1 = NetworkEdge('n1', 'n2', 123)
        edge2 = NetworkEdge('n3', 'n2', 234)
        edge3 = NetworkEdge('n2', 'n3', 234)
        network = Network.create('network', [edge1, edge2, edge3])
        self.assertEquals(4, network.num_edges())
        self.assertEquals(714, network.total_score())

    def test_normalize_to_same_score(self):
        """tests creating a network"""
        edge1 = NetworkEdge('n1', 'n2', 123)
        edge2 = NetworkEdge('n3', 'n2', 234)
        network = Network('network', [edge1, edge2])
        network.normalize_scores_to(357)
        self.assertEquals(357, network.total_score())

    def test_normalize_to_different_score(self):
        """tests creating a network"""
        edge1 = NetworkEdge('n1', 'n2', 123)
        edge2 = NetworkEdge('n3', 'n2', 234)
        network = Network('network', [edge1, edge2])
        network.normalize_scores_to(400)
        self.assertEquals(400, network.total_score())
