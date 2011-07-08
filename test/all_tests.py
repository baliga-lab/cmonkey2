"""Run all tests in the project"""
import unittest
from datatypes_test import DataMatrixTest, DataMatrixCollectionTest
from cmonkey_test import CMonkeyTest, MembershipTest
from util_test import DelimitedFileTest

if __name__ == '__main__':
    SUITE = []
    SUITE.append(unittest.TestLoader().loadTestsFromTestCase(DataMatrixTest))
    SUITE.append(unittest.TestLoader().loadTestsFromTestCase(CMonkeyTest))
    SUITE.append(unittest.TestLoader().loadTestsFromTestCase(MembershipTest))
    SUITE.append(unittest.TestLoader().loadTestsFromTestCase(
            DataMatrixCollectionTest))
    SUITE.append(unittest.TestLoader().loadTestsFromTestCase(
            DelimitedFileTest))
    unittest.TextTestRunner(verbosity=2).run(unittest.TestSuite(SUITE))
