"""network_test.py - unit tests for network module

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import unittest
import network as nw


class NetworkTest(unittest.TestCase):  # pylint: disable-msg=R0904
    """Test class for Network"""

    def test_construct(self):
        """tests creating a network using the constructor"""
        edge1 = ('n1', 'n2', 123)
        edge2 = ('n3', 'n2', 234)
        network = nw.Network.create('network', [edge1, edge2], 123, check_size=False)
        self.assertEquals('network', network.name)
        self.assertEquals(2, network.num_edges())
        self.assertTrue(edge1 in network.edges)
        self.assertTrue(edge2 in network.edges)
        self.assertEquals(714, network.total_score())
        self.assertEquals(123, network.weight)

    def test_create_unique_edges(self):
        """tests creating a network using the standard factory method"""
        edge1 = ('n1', 'n2', 123)
        edge2 = ('n3', 'n2', 234)
        network = nw.Network.create('network', [edge1, edge2], 234, check_size=False)
        self.assertEquals('network', network.name)
        self.assertEquals(2, network.num_edges())
        self.assertEquals(714, network.total_score())
        self.assertEquals(234, network.weight)

    def test_create_overlapping_edges(self):
        """tests creating a network using the standard factory method
        no duplicate edge will be generated"""
        edge1 = ('n1', 'n2', 123)
        edge2 = ('n3', 'n2', 234)
        edge3 = ('n2', 'n3', 234)
        network = nw.Network.create('network', [edge1, edge2, edge3], 456, check_size=False)
        self.assertEquals(2, network.num_edges())
        self.assertEquals(714, network.total_score())

    def test_normalize_to_same_score(self):
        """tests creating a network"""
        edge1 = ('n1', 'n2', 123)
        edge2 = ('n3', 'n2', 234)
        network = nw.Network.create('network', [edge1, edge2], 123, check_size=False)
        network.normalize_scores_to(357)
        self.assertEquals(357, network.total_score())

    def test_normalize_to_different_score(self):
        """tests creating a network"""
        edge1 = ('n1', 'n2', 123)
        edge2 = ('n3', 'n2', 234)
        network = nw.Network.create('network', [edge1, edge2], 234, check_size=False)
        network.normalize_scores_to(400)
        self.assertEquals(400, network.total_score())

    def test_edges_with_source_in(self):
        """tests the edges that have a source in the input"""
        edge1 = ('n1', 'n2', 123)
        edge2 = ('n3', 'n2', 234)
        edge3 = ('n4', 'n2', 234)
        edge4 = ('n4', 'n1', 123)
        network = nw.Network.create('network', [edge1, edge2, edge3, edge4], 42,
                                    check_size=False)
        res_edges = network.edges_with_node('n3')
        self.assertEquals(1, len(res_edges))
        self.assertTrue(edge2 in res_edges)
        
