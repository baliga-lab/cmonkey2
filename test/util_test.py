"""util_test.py - test classes for util module

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import unittest
import util
import numpy as np


class DelimitedFileTest(unittest.TestCase):  # pylint: disable-msg=R0904
    """Test class for DelimitedFile"""

    def test_read_with_tabs(self):
        """Reads a tab delimited file"""
        dfile = util.DelimitedFile.read("testdata/simple.tsv")
        lines = dfile.lines()
        self.assertEquals(["value11", "value12"], lines[0])
        self.assertEquals(["value21", "value22"], lines[1])
        self.assertIsNone(dfile.header())

    def test_read_with_tabs_and_header(self):
        """Reads a tab delimited file with a header"""
        dfile = util.DelimitedFile.read("testdata/simple.tsv", has_header=True)
        lines = dfile.lines()
        self.assertEquals(1, len(lines))
        self.assertEquals(["value11", "value12"], dfile.header())

    def test_read_with_semicolon_header_and_comments(self):
        """Reads a semicolon delimited file with a header and comments"""
        dfile = util.DelimitedFile.read("testdata/withcomments.ssv", sep=';',
                                        has_header=True, comment='#')
        lines = dfile.lines()
        self.assertEquals(2, len(lines))
        self.assertEquals(["header1", "header2"], dfile.header())

    def test_read_with_quotes(self):
        """Reads a semicolon delimited file with quotes"""
        dfile = util.DelimitedFile.read("testdata/withquotes.ssv", sep=';',
                                        has_header=False, comment='#', quote='"')
        lines = dfile.lines()
        self.assertEquals(["value11", "value12"], lines[0])
        self.assertEquals(["value21", "value22"], lines[1])

    def test_read_with_empty_lines(self):
        """Reads a semicolon delimited file containing emptylines"""
        dfile = util.DelimitedFile.read("testdata/withemptylines.ssv", sep=';',
                                        has_header=True, comment='#', quote='"')
        lines = dfile.lines()
        self.assertEquals(["header1", "header2"], dfile.header())
        self.assertEquals(2, len(lines))
        self.assertEquals(["value11", "value12"], lines[0])
        self.assertEquals(["value21", "value22"], lines[1])

    def test_create_from_text(self):
        """Reads a tab delimited file from a text"""
        dfile = util.DelimitedFile.create_from_text(
            "value11\tvalue12\nvalue21\tvalue22")
        lines = dfile.lines()
        self.assertEquals(["value11", "value12"], lines[0])
        self.assertEquals(["value21", "value22"], lines[1])
        self.assertIsNone(dfile.header())

    def test_create_from_text_empty_line_at_end(self):
        """Reads a tab delimited file from a text"""
        dfile = util.DelimitedFile.create_from_text(
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
        self.assertEquals(3, util.levenshtein_distance('sitting', 'kitten'))

    def test_saturday_sunday(self):
        """compare Saturday with Sunday"""
        self.assertEquals(3, util.levenshtein_distance('Sunday', 'Saturday'))


RSAT_LIST_FILE_PATH = "testdata/RSAT_genomes_listing.txt"


class BestMatchingLinksTest(unittest.TestCase):  # pylint: disable-msg=R0904
    """Test class for best_matching_links"""

    def test_best_rsat_matches(self):
        """test the best_matching_links function"""
        with open(RSAT_LIST_FILE_PATH) as inputfile:
            html = inputfile.read()
        matches = util.best_matching_links('Halobacterium', html)
        self.assertEquals("Halobacterium_sp/", matches[0])


class UtilsTest(unittest.TestCase):  # pylint: disable-msg=R0904
    """Test class for utility functions"""

    def test_make_matrix(self):
        """tests the make_matrix function"""
        matrix = util.make_matrix(["row1", "row2"], 3)
        self.assertEquals(2, len(matrix))
        self.assertEquals(3, len(matrix['row1']))
        self.assertEquals(0, matrix['row1'][0])

    def test_quantile(self):
        """tests the quantile function"""
        data = [1, 2, 3, 4, 5]
        self.assertEquals(1, util.quantile(data, 0))
        self.assertEquals(1.8, util.quantile(data, 0.2))
        self.assertEquals(2, util.quantile(data, 0.25))
        self.assertEquals(3, util.quantile(data, 0.5))
        self.assertEquals(4, util.quantile(data, 0.75))
        self.assertEquals(5, util.quantile(data, 1))

    def test_r_stddev(self):
        """tests the standard deviation function"""
        self.assertEquals(0.1, util.r_stddev([0.1, 0.2, 0.3]))

    def test_r_variance_columns(self):
        """tests the column variance function"""
        matrix = [[0.0010, 0.1234, 0.21370, 0.0342],
                  [0.2123, -0.2135, -0.99980, -0.0213],
                  [-0.4534, 0.5546, 0.79123, 0.00312321]]
        result = util.r_variance_columns(matrix)
        self.assertAlmostEqual(0.1157139233, result[0])
        self.assertAlmostEqual(0.1482354433, result[1])
        self.assertAlmostEqual(0.8356519353, result[2])
        self.assertAlmostEqual(0.0007737517, result[3])

    def test_column_means(self):
        """tests the column_means() function"""
        matrix = [[0.0010, 0.1234, 0.21370, 0.0342],
                  [0.2123, -0.2135, -0.99980, -0.0213],
                  [-0.4534, 0.5546, 0.79123, 0.00312321]]
        result = util.column_means(matrix)
        self.assertAlmostEqual(-0.08003333, result[0])
        self.assertAlmostEqual(0.15483333, result[1])
        self.assertAlmostEqual(0.00171, result[2])
        self.assertAlmostEqual(0.00534107, result[3])

    def test_column_means_with_nans(self):
        """tests the column_means() function, containing NaNs"""
        matrix = [[0.0010, 0.1234, 0.21370, np.nan],
                  [0.2123, np.nan, -0.99980, -0.0213],
                  [np.nan, 0.5546, 0.79123, 0.00312321]]
        result = util.column_means(matrix)
        self.assertAlmostEqual(0.10664999, result[0])
        self.assertAlmostEqual(0.33899999, result[1])
        self.assertAlmostEqual(0.00171, result[2])
        self.assertAlmostEqual(-0.00908839499, result[3])

    def test_row_means(self):
        """tests the row_means() function"""
        matrix = [[0.0010, 0.1234, 0.21370, 0.0342],
                  [0.2123, -0.2135, -0.99980, -0.0213],
                  [-0.4534, 0.5546, 0.79123, 0.00312321]]
        result = util.row_means(matrix)
        self.assertAlmostEqual(0.0930750, result[0])
        self.assertAlmostEqual(-0.255575, result[1])
        self.assertAlmostEqual(0.2238883025, result[2])

    def test_row_means_with_nans(self):
        """tests the row_means() function"""
        matrix = [[0.0010, np.nan, 0.21370, 0.0342],
                  [0.2123, -0.2135, -0.99980, -0.0213],
                  [-0.4534, 0.5546, 0.79123, np.nan]]
        result = util.row_means(matrix)
        self.assertAlmostEqual(0.08296666, result[0])
        self.assertAlmostEqual(-0.255575, result[1])
        self.assertAlmostEqual(0.297476666, result[2])

    def test_trim_mean_nonmedian(self):
        self.assertAlmostEqual(
            40.625,
            util.trim_mean([2, 4, 6, 7, 11, 21, 81, 90, 105, 121], 0.1))

    def test_trim_mean_median(self):
        self.assertAlmostEqual(3.5, util.trim_mean([.1, .2, 3, 4, 5, 6], 0.5))

    def test_trim_mean_no_values(self):
        self.assertEqual(0, util.trim_mean([], 0.05))

    def test_trim_mean_real(self):
        values = [0.0, 0.0, -8.7520618359684352, -8.7520618359684352, 0.0, 0.0,
                  0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        self.assertAlmostEqual(-1.4586770, util.trim_mean(values, 0.05))

    def test_mean_with_nans(self):
        """tests the mean() function"""
        array = np.array([2.0, 3.0, np.nan, 1.0])
        result = util.mean(array)
        self.assertAlmostEqual(2.0, result)

    def test_density(self):
        kvalues = [3.4268700450682301, 3.3655160468930152, -8.0654569044842539,
                   2.0762815314005487, 4.8537715329554203, 1.2374476248622075]
        cluster_values = [-3.5923001345962162, 0.77069901513184735,
                           -4.942909785931378, -3.1580950032999096] 
        bandwidth = 2.69474878768
        dmin = -13.8848342423
        dmax = 12.6744452247
        result = util.density(kvalues, cluster_values, bandwidth, dmin, dmax)
        self.assertAlmostEquals(0.08663036966690765, result[0])
        self.assertAlmostEquals(0.08809242907902183, result[1])
        self.assertAlmostEquals(0.49712338305039777, result[2])
        self.assertAlmostEquals(0.12248549621579163, result[3])
        self.assertAlmostEquals(0.05708884005243133, result[4])
        self.assertAlmostEquals(0.14857948193544993, result[5])

    def test_sd_rnorm(self):
        result = util.sd_rnorm([1.3, 1.6, 1.2, 1.05], 9, 0.748951)
        # the results are fairly random, make sure we have the right
        # number of values
        self.assertEquals(9, len(result))
        

class Order2StringTest(unittest.TestCase):  # pylint: disable-msg=R09042
    """Test class for order2string"""

    def test_order2string(self):
        self.assertEquals("1st", util.order2string(1))
        self.assertEquals("2nd", util.order2string(2))
        self.assertEquals("3rd", util.order2string(3))
        self.assertEquals("4th", util.order2string(4))
        self.assertEquals("11th", util.order2string(11))
        self.assertEquals("12th", util.order2string(12))
        self.assertEquals("21st", util.order2string(21))
        self.assertEquals("22nd", util.order2string(22))
        self.assertEquals("23rd", util.order2string(23))
