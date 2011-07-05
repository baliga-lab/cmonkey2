"""unit test module for cmonkey"""
import unittest
from cmonkey import CMonkey, Membership

class CMonkeyTest(unittest.TestCase):
    """Test class for CMonkey"""

    def test_create_cmonkey(self):
        """create a CMonkey object"""
        cmonkey = CMonkey([])
        self.assertFalse(cmonkey.run_finished)
        self.assertEquals('hpy', cmonkey.configuration['organism'])

    def test_create_cmonkey_with_config(self):
        """create CMonkey object, specifying an organism"""
        config = {'organism': 'homo sapiens'}
        cmonkey = CMonkey([], config)
        self.assertFalse(cmonkey.run_finished)
        self.assertEquals('homo sapiens', cmonkey.configuration['organism'])

class MembershipTest(unittest.TestCase):
    """Test class for Membership"""

    def test_map_to_is_member_matrix(self):
        in_matrix = [[1, 2],[2, 3]]
        out = Membership.map_to_is_member_matrix(in_matrix, 3)
        self.assertEquals([[True, False], [True, True], [False, True]], out)

if __name__ == '__main__':
    unittest.main()
