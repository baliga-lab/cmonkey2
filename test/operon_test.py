"""util_test.py - test classes for operon module

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import unittest
from operons import read_from_microbes_online


class MockMicrobesOnline:
    """mock class for MicrobesOnline"""

    def __init__(self, filename):
        """creates a mock instance"""
        with open(filename) as inputfile:
            self.content = inputfile.read()

    def get_operon_predictions_for(self, _):
        """mocked MicrobesOnline operon prediction result"""
        return self.content


class ReadMicrobesOnlineTest(unittest.TestCase):  # pylint: disable-msg=R0904
    """Test class for read_from_microbes_online"""

    def test_read_success(self):
        """test happy path"""
        microbes_online = MockMicrobesOnline('testdata/gnc64091.named')
        preds = read_from_microbes_online(microbes_online, '64091')
        self.assertEquals(7, len(preds))
