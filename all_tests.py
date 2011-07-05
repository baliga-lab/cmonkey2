"""Run all tests in the project"""
import unittest
from datatypes_test import DataMatrixTest
from cmonkey_test import CMonkeyTest, MembershipTest

if __name__ == '__main__':
    SUITE = []
    SUITE.append(unittest.TestLoader().loadTestsFromTestCase(DataMatrixTest))
    SUITE.append(unittest.TestLoader().loadTestsFromTestCase(CMonkeyTest))
    SUITE.append(unittest.TestLoader().loadTestsFromTestCase(MembershipTest))
    unittest.TextTestRunner(verbosity=2).run(unittest.TestSuite(SUITE))

