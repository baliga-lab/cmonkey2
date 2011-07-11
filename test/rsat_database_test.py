"""rsat_database_test.py - unit tests for RsatDatabase class

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import unittest
from organism import RsatDatabase


class RsatDatabaseTest(unittest.TestCase):  # pylint: disable-msg=R0904
    """Test class for RsatDatabase. These tests interacts with a real web
    service and will therefore run relatively slowly. Should be run
    as part of an integration test suite and not of the unit test suite"""

    def setUp(self):  # pylint: disable-msg=C0103
        """test fixture"""
        self.database = RsatDatabase('http://rsat.ccb.sickkids.ca')

    def test_get_directory_html(self):
        """test get_directory_html method"""
        html = self.database.get_directory_html()
        self.assertIsNotNone(html)

    def test_get_organism_file(self):
        """test get_directory_html method"""
        text = self.database.get_organism_file('Helicobacter_pylori_26695')
        self.assertIsNotNone(text)
