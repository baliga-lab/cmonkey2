"""all_tests.py - run all unit tests in the project

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import unittest
from datamatrix_test import DataMatrixTest, DataMatrixCollectionTest
from datamatrix_test import DataMatrixFactoryTest, NoChangeFilterTest
from datamatrix_test import CenterScaleFilterTest
from cmonkey_test import CMonkeyTest, MembershipTest
from util_test import DelimitedFileTest, UtilsTest
from util_test import LevenshteinDistanceTest, BestMatchingLinksTest
from organism_test import KeggOrganismCodeMapperTest, RsatOrganismMapperTest
from organism_test import GoTaxonomyMapperTest, OrganismTest
from thesaurus_test import DelimitedFileFactoryTest

# pylint: disable-msg=C0301
if __name__ == '__main__':
    SUITE = []
    SUITE.append(unittest.TestLoader().loadTestsFromTestCase(DataMatrixTest))
    SUITE.append(unittest.TestLoader().loadTestsFromTestCase(DataMatrixCollectionTest))
    SUITE.append(unittest.TestLoader().loadTestsFromTestCase(DataMatrixFactoryTest))
    SUITE.append(unittest.TestLoader().loadTestsFromTestCase(NoChangeFilterTest))
    SUITE.append(unittest.TestLoader().loadTestsFromTestCase(CenterScaleFilterTest))

    SUITE.append(unittest.TestLoader().loadTestsFromTestCase(CMonkeyTest))
    SUITE.append(unittest.TestLoader().loadTestsFromTestCase(MembershipTest))

    SUITE.append(unittest.TestLoader().loadTestsFromTestCase(DelimitedFileTest))
    SUITE.append(unittest.TestLoader().loadTestsFromTestCase(UtilsTest))
    SUITE.append(unittest.TestLoader().loadTestsFromTestCase(LevenshteinDistanceTest))
    SUITE.append(unittest.TestLoader().loadTestsFromTestCase(BestMatchingLinksTest))

    SUITE.append(unittest.TestLoader().loadTestsFromTestCase(KeggOrganismCodeMapperTest))
    SUITE.append(unittest.TestLoader().loadTestsFromTestCase(RsatOrganismMapperTest))
    SUITE.append(unittest.TestLoader().loadTestsFromTestCase(GoTaxonomyMapperTest))
    SUITE.append(unittest.TestLoader().loadTestsFromTestCase(OrganismTest))

    SUITE.append(unittest.TestLoader().loadTestsFromTestCase(DelimitedFileFactoryTest))

    unittest.TextTestRunner(verbosity=2).run(unittest.TestSuite(SUITE))
