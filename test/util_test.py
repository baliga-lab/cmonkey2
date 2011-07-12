"""util_test.py - test classes for util module

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import unittest
from util import DelimitedFile, levenshtein_distance, best_matching_links
from util import quantile, make_matrix, r_stddev


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
        with open(RSAT_LIST_FILE_PATH) as inputfile:
            html = inputfile.read()
        matches = best_matching_links('Halobacterium', html)
        self.assertEquals("Halobacterium_sp/", matches[0])


class UtilsTest(unittest.TestCase):  # pylint: disable-msg=R0904
    """Test class for utility functions"""

    def test_make_matrix(self):
        """tests the make_matrix function"""
        matrix = make_matrix(["row1", "row2"], 3)
        self.assertEquals(2, len(matrix))
        self.assertEquals(3, len(matrix['row1']))
        self.assertEquals(0, matrix['row1'][0])

    def test_quantile(self):
        """tests the quantile function"""
        data = [1, 2, 3, 4, 5]
        self.assertEquals(1, quantile(data, 0))
        self.assertEquals(1.8, quantile(data, 0.2))
        self.assertEquals(2, quantile(data, 0.25))
        self.assertEquals(3, quantile(data, 0.5))
        self.assertEquals(4, quantile(data, 0.75))
        self.assertEquals(5, quantile(data, 1))

    def test_r_stddev(self):
        """tests the standard deviation function"""
        self.assertEquals(0.1, r_stddev([0.1, 0.2, 0.3]))
