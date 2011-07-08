"""Test classes for util module"""
import unittest
from util import DelimitedFile, levenshtein_distance, best_matching_links


class DelimitedFileTest(unittest.TestCase):  # pylint: disable-msg=R0904
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

    def test_read_with_empty_lines(self):
        """Reads a semicolon delimited file containing emptylines"""
        dfile = DelimitedFile.read("testdata/withemptylines.ssv", sep=';',
                                   has_header=True, comment='#', quote='"')
        lines = dfile.get_lines()
        self.assertEquals(["header1", "header2"], dfile.get_header())
        self.assertEquals(2, len(lines))
        self.assertEquals(["value11", "value12"], lines[0])
        self.assertEquals(["value21", "value22"], lines[1])


class LevenshteinDistanceTest(unittest.TestCase):  # pylint: disable-msg=R0904
    """Test class for levenshtein_distance"""

    def test_kitten_sitting(self):
        """compare kitten with sitting"""
        self.assertEquals(3, levenshtein_distance('sitting', 'kitten'))

    def test_saturday_sunday(self):
        """compare Saturday with Sunday"""
        self.assertEquals(3, levenshtein_distance('Sunday', 'Saturday'))


RSAT_LIST_FILE_PATH = "testdata/RSAT_genomes_listing.txt"


class BestMatchingLinksTest(unittest.TestCase):  # pylint: disable-msg=R0904
    """Test class for best_matching_links"""

    def test_best_rsat_matches(self):
        """test the best_matching_links function"""
        matches = best_matching_links('Halobacterium', RSAT_LIST_FILE_PATH)
        self.assertEquals("Halobacterium_sp/", matches[0])
