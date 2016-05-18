"""rsat_test.py - unit tests for RsatDatabase class

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import unittest
import os
import shutil
import cmonkey.rsat as rsat
import cmonkey.util as util


RSAT_BASE_URL = 'http://rsat01.biologie.ens.fr/rsat'


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
        self.database = rsat.RsatDatabase(RSAT_BASE_URL, 'testcache',
                                          'Helicobacter_pylori_26695_uid57787', 85962)

    def tearDown(self):  # pylint: disable-msg=C0103
        """test cleanup"""
        if os.path.exists('testcache'):
            shutil.rmtree('testcache')

    def test_get_organism_names(self):
        """test get_organism_names method"""
        self.assertEquals("85962", self.database.get_taxonomy_id('Helicobacter_pylori_26695_uid57787'))

    def test_get_features(self):
        """test get_features method"""
        text = self.database.get_features('Helicobacter_pylori_26695_uid57787')
        self.assertIsNotNone(text)

    def test_get_feature_names(self):
        """test get_feature_names method"""
        text = self.database.get_feature_names('Helicobacter_pylori_26695_uid57787')
        self.assertIsNotNone(text)
