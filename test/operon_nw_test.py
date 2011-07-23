"""util_test.py - test classes for operon module

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import unittest
from microbes_online import get_network_factory


class MockMicrobesOnline:
    """mock class for MicrobesOnline"""

    def __init__(self, filename):
        """creates a mock instance"""
        with open(filename) as inputfile:
            self.content = inputfile.read()

    def get_operon_predictions_for(self, _):
        """mocked MicrobesOnline operon prediction result"""
        return self.content


class MockOrganism:
    """mock Organism class"""
    def __init__(self, taxonomy_id):
        self.__taxonomy_id = taxonomy_id

    def taxonomy_id(self):
        return self.__taxonomy_id


class ReadOperonNetworkTest(unittest.TestCase):  # pylint: disable-msg=R0904
    """Test class for read_from_microbes_online"""

    def test_get_network_factory(self):
        """test happy path"""
        microbes_online = MockMicrobesOnline('testdata/gnc64091.named')
        network = get_network_factory(microbes_online)(MockOrganism('64091'))
        self.assertEquals(14, network.num_edges())
        self.assertEquals(14000, network.total_score())
