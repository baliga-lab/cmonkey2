"""rsat_test.py - unit tests for RsatDatabase class

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import unittest
import os
import shutil
import rsat
import util


class RsatDatabaseTest(unittest.TestCase):  # pylint: disable-msg=R0904
    """Test class for RsatDatabase. These tests interacts with a real web
    service and will therefore run relatively slowly. Should be run
    as part of an integration test suite and not of the unit test suite.
    There is no real attempt to actually check the contents of the files,
    it is mainly a check for link availability.
    """

    def setUp(self):  # pylint: disable-msg=C0103
        """test fixture"""
        if not os.path.exists('testcache'):
            os.mkdir('testcache')
        self.database = rsat.RsatDatabase('http://rsat.ccb.sickkids.ca',
                                          'testcache')

    def tearDown(self):  # pylint: disable-msg=C0103
        """test cleanup"""
        if os.path.exists('testcache'):
            shutil.rmtree('testcache')

    def test_get_directory(self):
        """test get_directory method"""
        html = self.database.get_directory()
        self.assertIsNotNone(html)

    def test_get_organism(self):
        """test get_organism method"""
        text = self.database.get_organism('Helicobacter_pylori_26695')
        self.assertIsNotNone(text)

    def test_get_organism_names(self):
        """test get_organism_names method"""
        text = self.database.get_organism_names(
            'Helicobacter_pylori_26695')
        self.assertIsNotNone(text)

    def test_get_ensembl_organism_names(self):
        """test get_ensembl_organism_names method"""
        text = self.database.get_ensembl_organism_names(
            'Saccharomyces_cerevisiae')
        self.assertIsNotNone(text)

    def test_get_ensembl_organism_names_404(self):
        """test get_ensembl_organism_names method with not found"""
        self.assertRaises(util.DocumentNotFound,
                          self.database.get_ensembl_organism_names,
                          'nonexist')

    def test_get_features(self):
        """test get_features method"""
        text = self.database.get_features('Helicobacter_pylori_26695')
        self.assertIsNotNone(text)

    def test_get_feature_names(self):
        """test get_feature_names method"""
        text = self.database.get_feature_names('Helicobacter_pylori_26695')
        self.assertIsNotNone(text)
