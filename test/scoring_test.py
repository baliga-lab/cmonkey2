"""scoring_test.py - test classes for scoring module

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import unittest
import scoring
import numpy as np


class DefaultScalingTest(unittest.TestCase):  # pylint: disable-msg=R0904
    """Test class for default scaling"""

    def test_get_default_motif_scaling(self):
        scaling_fun = scoring.get_default_motif_scaling(100)
        self.assertAlmostEqual(0.0, scaling_fun(1))
        self.assertAlmostEqual(0.02702703, scaling_fun(3))
        self.assertAlmostEqual(0.98648649, scaling_fun(74))
        self.assertAlmostEqual(1.0, scaling_fun(75))
        self.assertAlmostEqual(1.0, scaling_fun(76))
        self.assertAlmostEqual(1.0, scaling_fun(100))

    def test_get_default_network_scaling(self):
        scaling_fun = scoring.get_default_network_scaling(100)
        self.assertAlmostEqual(0.0, scaling_fun(1))
        self.assertAlmostEqual(0.013513514, scaling_fun(3))
        self.assertAlmostEqual(0.493243243, scaling_fun(74))
        self.assertAlmostEqual(0.5, scaling_fun(75))
        self.assertAlmostEqual(0.5, scaling_fun(76))
        self.assertAlmostEqual(0.5, scaling_fun(100))


if __name__ == '__main__':
    SUITE = []
    SUITE.append(unittest.TestLoader().loadTestsFromTestCase(DefaultScalingTest))
