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

    def test_create(self):
        """tests creating a network"""
        edge1 = NetworkEdge('n1', 'n2', 123)
        edge2 = NetworkEdge('n3', 'n2', 234)
        network = Network([edge1, edge2])
        self.assertEquals(2, network.num_edges())
        self.assertEquals(357, network.total_score())
