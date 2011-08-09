"""rsat_test.py - unit tests for RsatDatabase class

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import unittest
import os
import shutil
import microbes_online as mo


class MicrobesOnlineTest(unittest.TestCase):  # pylint: disable-msg=R0904
    """Test class for accessing the actual Microbes Online Service
    Don't use in TDD unit tests"""

    def setUp(self):  # pylint: disable-msg=C0103
        """test fixture"""
        if not os.path.exists('testcache'):
            os.mkdir('testcache')
        self.service = mo.MicrobesOnline(mo.MICROBES_ONLINE_BASE_URL, 'testcache')

    def tearDown(self):  # pylint: disable-msg=C0103
        """test cleanup"""
        if os.path.exists('testcache'):
            shutil.rmtree('testcache')

    def test_get_predictions_success(self):
        """doing a successful retrieval"""
        preds = self.service.get_operon_predictions_for('64091')
        self.assertIsNotNone(preds)
