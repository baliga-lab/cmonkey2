"""unit test module for cmonkey"""
import unittest
from cmonkey import CMonkey

class CMonkeyTest(unittest.TestCase):
    """Test class for CMonkey"""

    def test_create_cmonkey(self):
        """create a CMonkey object"""
        cmonkey = CMonkey([])
        self.assertFalse(cmonkey.run_finished)
        self.assertEquals('hpy', cmonkey.organism)

    def test_create_cmonkey_with_organism(self):
        """create CMonkey object, specifying an organism"""
        cmonkey = CMonkey([], 'homo sapiens')
        self.assertFalse(cmonkey.run_finished)
        self.assertEquals('homo sapiens', cmonkey.organism)

if __name__ == '__main__':
    unittest.main()
