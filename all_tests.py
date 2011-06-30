import unittest
from datatypes_test import *
from cmonkey_test import *

if __name__ == '__main__':
  suite = []
  suite.append(unittest.TestLoader().loadTestsFromTestCase(DataMatrixTest))
  suite.append(unittest.TestLoader().loadTestsFromTestCase(cMonkeyTest))
  unittest.TextTestRunner(verbosity=2).run(unittest.TestSuite(suite))

