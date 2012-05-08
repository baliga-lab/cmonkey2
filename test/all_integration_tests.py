"""all_integration_tests.py - run all integration tests in the project

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import unittest
import rsat_test
import meme430_test
import halo_genes_test
import microbes_online_test as mo_test


if __name__ == '__main__':
    SUITE = []
    SUITE.append(unittest.TestLoader().loadTestsFromTestCase(
            rsat_test.RsatDatabaseTest))
    SUITE.append(unittest.TestLoader().loadTestsFromTestCase(
            mo_test.MicrobesOnlineTest))
    SUITE.append(unittest.TestLoader().loadTestsFromTestCase(
            meme430_test.Meme430Test))
    SUITE.append(unittest.TestLoader().loadTestsFromTestCase(
            halo_genes_test.HaloGeneTest))

    unittest.TextTestRunner(verbosity=2).run(unittest.TestSuite(SUITE))
