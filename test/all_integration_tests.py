"""all_integration_tests.py - run all integration tests in the project

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import unittest
from rsat_test import RsatDatabaseTest
from microbes_online_test import MicrobesOnlineTest

if __name__ == '__main__':
    SUITE = []
    SUITE.append(unittest.TestLoader().loadTestsFromTestCase(
            RsatDatabaseTest))
    SUITE.append(unittest.TestLoader().loadTestsFromTestCase(
            MicrobesOnlineTest))

    unittest.TextTestRunner(verbosity=2).run(unittest.TestSuite(SUITE))
