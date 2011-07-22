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
from organism_test import GoTaxonomyMapperTest, OrganismFactoryTest
from organism_test import OrganismTest
from seqtools_test import SeqtoolsTest

from thesaurus_test import DelimitedFileFactoryTest
from operon_test import ReadMicrobesOnlineTest
from network_test import NetworkEdgeTest, NetworkTest

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
    SUITE.append(unittest.TestLoader().loadTestsFromTestCase(OrganismFactoryTest))
    SUITE.append(unittest.TestLoader().loadTestsFromTestCase(OrganismTest))

    SUITE.append(unittest.TestLoader().loadTestsFromTestCase(SeqtoolsTest))

    SUITE.append(unittest.TestLoader().loadTestsFromTestCase(DelimitedFileFactoryTest))

    SUITE.append(unittest.TestLoader().loadTestsFromTestCase(ReadMicrobesOnlineTest))

    SUITE.append(unittest.TestLoader().loadTestsFromTestCase(NetworkEdgeTest))
    SUITE.append(unittest.TestLoader().loadTestsFromTestCase(NetworkTest))

    unittest.TextTestRunner(verbosity=2).run(unittest.TestSuite(SUITE))
