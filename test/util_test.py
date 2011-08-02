"""util_test.py - test classes for util module

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import unittest
from util import DelimitedFile, levenshtein_distance, best_matching_links
from util import quantile, make_matrix, r_stddev, r_variance_columns
from util import column_means, row_means


class DelimitedFileTest(unittest.TestCase):  # pylint: disable-msg=R0904
    """Test class for DelimitedFile"""

    def test_read_with_tabs(self):
        """Reads a tab delimited file"""
        dfile = DelimitedFile.read("testdata/simple.tsv")
        lines = dfile.lines()
        self.assertEquals(["value11", "value12"], lines[0])
        self.assertEquals(["value21", "value22"], lines[1])
        self.assertIsNone(dfile.header())

    def test_read_with_tabs_and_header(self):
        """Reads a tab delimited file with a header"""
        dfile = DelimitedFile.read("testdata/simple.tsv", has_header=True)
        lines = dfile.lines()
        self.assertEquals(1, len(lines))
        self.assertEquals(["value11", "value12"], dfile.header())

    def test_read_with_semicolon_header_and_comments(self):
        """Reads a semicolon delimited file with a header and comments"""
        dfile = DelimitedFile.read("testdata/withcomments.ssv", sep=';',
                                    has_header=True, comment='#')
        lines = dfile.lines()
        self.assertEquals(2, len(lines))
        self.assertEquals(["header1", "header2"], dfile.header())

    def test_read_with_quotes(self):
        """Reads a semicolon delimited file with quotes"""
        dfile = DelimitedFile.read("testdata/withquotes.ssv", sep=';',
                                   has_header=False, comment='#', quote='"')
        lines = dfile.lines()
        self.assertEquals(["value11", "value12"], lines[0])
        self.assertEquals(["value21", "value22"], lines[1])

    def test_read_with_empty_lines(self):
        """Reads a semicolon delimited file containing emptylines"""
        dfile = DelimitedFile.read("testdata/withemptylines.ssv", sep=';',
                                   has_header=True, comment='#', quote='"')
        lines = dfile.lines()
        self.assertEquals(["header1", "header2"], dfile.header())
        self.assertEquals(2, len(lines))
        self.assertEquals(["value11", "value12"], lines[0])
        self.assertEquals(["value21", "value22"], lines[1])

    def test_create_from_text(self):
        """Reads a tab delimited file from a text"""
        dfile = DelimitedFile.create_from_text(
            "value11\tvalue12\nvalue21\tvalue22")
        lines = dfile.lines()
        self.assertEquals(["value11", "value12"], lines[0])
        self.assertEquals(["value21", "value22"], lines[1])
        self.assertIsNone(dfile.header())

    def test_create_from_text_empty_line_at_end(self):
        """Reads a tab delimited file from a text"""
        dfile = DelimitedFile.create_from_text(
            "value11\tvalue12\nvalue21\tvalue22\n")
        lines = dfile.lines()
        self.assertEquals(2, len(lines))
        self.assertEquals(["value11", "value12"], lines[0])
        self.assertEquals(["value21", "value22"], lines[1])
        self.assertIsNone(dfile.header())


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

    def test_r_variance_columns(self):
        """tests the column variance function"""
        matrix = [[0.0010, 0.1234, 0.21370, 0.0342],
                  [0.2123, -0.2135, -0.99980, -0.0213],
                  [-0.4534, 0.5546, 0.79123, 0.00312321]]
        result = r_variance_columns(matrix)
        self.assertAlmostEqual(0.1157139233, result[0])
        self.assertAlmostEqual(0.1482354433, result[1])
        self.assertAlmostEqual(0.8356519353, result[2])
        self.assertAlmostEqual(0.0007737517, result[3])

    def test_column_means(self):
        """tests the column_means() function"""
        matrix = [[0.0010, 0.1234, 0.21370, 0.0342],
                  [0.2123, -0.2135, -0.99980, -0.0213],
                  [-0.4534, 0.5546, 0.79123, 0.00312321]]
        result = column_means(matrix)
        self.assertAlmostEqual(-0.08003333, result[0])
        self.assertAlmostEqual(0.15483333, result[1])
        self.assertAlmostEqual(0.00171, result[2])
        self.assertAlmostEqual(0.00534107, result[3])

    def test_row_means(self):
        """tests the row_means() function"""
        matrix = [[0.0010, 0.1234, 0.21370, 0.0342],
                  [0.2123, -0.2135, -0.99980, -0.0213],
                  [-0.4534, 0.5546, 0.79123, 0.00312321]]
        result = row_means(matrix)
        self.assertAlmostEqual(0.0930750, result[0])
        self.assertAlmostEqual(-0.255575, result[1])
        self.assertAlmostEqual(0.2238883025, result[2])
