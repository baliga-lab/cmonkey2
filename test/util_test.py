"""Test classes for util module"""
import unittest
from util import DelimitedFile

class DelimitedFileTest(unittest.TestCase):
    """Test class for DelimitedFile"""

    def test_read_with_tabs(self):
        """Reads a tab delimited file"""
        dfile = DelimitedFile.read("testdata/simple.tsv")
        lines = dfile.get_lines()
        self.assertEquals(["value11", "value12"], lines[0])
        self.assertEquals(["value21", "value22"], lines[1])
        self.assertIsNone(dfile.get_header())

    def test_read_with_tabs_and_header(self):
        """Reads a tab delimited file with a header"""
        dfile = DelimitedFile.read("testdata/simple.tsv", has_header=True)
        lines = dfile.get_lines()
        self.assertEquals(1, len(lines))
        self.assertEquals(["value11", "value12"], dfile.get_header())

    def test_read_with_semicolon_header_and_comments(self):
        """Reads a semicolon delimited file with a header and comments"""
        dfile = DelimitedFile.read("testdata/withcomments.ssv", sep=';',
                                    has_header=True, comment='#')
        lines = dfile.get_lines()
        self.assertEquals(2, len(lines))
        self.assertEquals(["header1", "header2"], dfile.get_header())

    def test_read_with_quotes(self):
        """Reads a semicolon delimited file with quotes"""
        dfile = DelimitedFile.read("testdata/withquotes.ssv", sep=';',
                                   has_header=False, comment='#', quote='"')
        lines = dfile.get_lines()
        self.assertEquals(["value11", "value12"], lines[0])
        self.assertEquals(["value21", "value22"], lines[1])
