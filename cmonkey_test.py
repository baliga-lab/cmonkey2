import unittest
from cmonkey import *

class cMonkeyTest(unittest.TestCase):

    def test_create_cmonkey(self):
        cm = cMonkey([])
        self.assertFalse(cm.run_finished)
        self.assertEquals('hpy', cm.organism)

    def test_create_cmonkey_with_organism(self):
        cm = cMonkey([], 'homo sapiens')
        self.assertFalse(cm.run_finished)
        self.assertEquals('homo sapiens', cm.organism)

if __name__ == '__main__':
    unittest.main()
