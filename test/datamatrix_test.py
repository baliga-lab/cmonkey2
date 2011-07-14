"""datatypes_test.py - test classes for datatypes module

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import unittest
from datamatrix import DataMatrix, DataMatrixCollection, DataMatrixFactory
from datamatrix import nochange_filter, center_scale_filter
from copy import deepcopy


class DataMatrixTest(unittest.TestCase):  # pylint: disable-msg=R0904
    """Test class for DataMatrix"""

    def test_create_with_0_row_size(self):
        """create DataMatrix with a 0 row size"""
        matrix = DataMatrix(0, 3)
        self.assertEquals(0, matrix.num_rows())
        self.assertEquals(0, matrix.num_columns())

    def test_create_without_names(self):
        """create DataMatrix without row and column names"""
        matrix = DataMatrix(3, 4)
        self.assertEquals(3, matrix.num_rows())
        self.assertEquals(4, matrix.num_columns())
        self.assertEquals(0.0, matrix.value_at(0, 0))
        self.assertEquals("Row 0", matrix.row_name(0))
        self.assertEquals("Row 1", matrix.row_name(1))
        self.assertEquals("Column 0", matrix.column_name(0))
        self.assertEquals("Column 1", matrix.column_name(1))

    def test_create_with_names(self):
        """create DataMatrix with row and column names"""
        matrix = DataMatrix(3, 2, ["MyRow1", "MyRow2", "MyRow3"],
                            ["MyCol1", "MyCol2"])
        self.assertEquals(3, matrix.num_rows())
        self.assertEquals(2, matrix.num_columns())
        self.assertEquals(0.0, matrix.value_at(0, 0))
        self.assertEquals("MyRow1", matrix.row_name(0))
        self.assertEquals("MyRow2", matrix.row_name(1))
        self.assertEquals("MyCol1", matrix.column_name(0))
        self.assertEquals("MyCol2", matrix.column_name(1))
        self.assertIsNotNone(str(matrix))

    def test_create_with_wrong_row_name_count(self):
        """create DataMatrix, providing the wrong number of row names"""
        self.assertRaises(ValueError, DataMatrix,
                          3, 2, row_names=["MyRow1", "MyRow2"])

    def test_create_with_wrong_column_name_count(self):
        """create DataMatrix, providing the wrong number of column names"""
        self.assertRaises(ValueError, DataMatrix,
                          3, 2, col_names=["MyCol1"])

    def test_set_value(self):
        """set a value in the matrix"""
        matrix = DataMatrix(3, 4)
        matrix.set_value_at(0, 1, 42.0)
        self.assertEquals(42.0, matrix.value_at(0, 1))

    def test_set_values_none(self):
        """invoke set_values() with None should result in ValueError"""
        matrix = DataMatrix(2, 2)
        self.assertRaises(ValueError, matrix.set_values, None)

    def test_set_values_wrong_row_count(self):
        """invoke set_values() with wrong row count should result in
        ValueError"""
        matrix = DataMatrix(2, 2)
        self.assertRaises(ValueError, matrix.set_values, [[1, 2]])

    def test_set_values_wrong_column_count(self):
        """invoke set_values() with wrong row count should result in
        ValueError"""
        matrix = DataMatrix(2, 2)
        self.assertRaises(ValueError, matrix.set_values,
                          [[1, 2], [2, 3, 4]])

    def test_set_values_ok(self):
        """invoke set_values() with wrong row count should result in
        ValueError"""
        matrix = DataMatrix(2, 2)
        matrix.set_values([[1, 2], [2, 3]])
        self.assertEquals([[1, 2], [2, 3]], matrix.values())

    def test_get_row_values(self):
        """tests the get_row_values"""
        matrix = DataMatrix(2, 2)
        matrix.set_values([[1, 2], [2, 3]])
        self.assertEquals([1, 2], matrix.get_row_values(0))
        self.assertEquals([2, 3], matrix.get_row_values(1))

    def test_get_column_values(self):
        """tests the get_column_values"""
        matrix = DataMatrix(2, 2)
        matrix.set_values([[3, 7], [8, 9]])
        self.assertEquals([3, 8], matrix.get_column_values(0))
        self.assertEquals([7, 9], matrix.get_column_values(1))


class DataMatrixCollectionTest(unittest.TestCase):  # pylint: disable-msg=R0904
    """Test class for MatrixCollection"""

    def test_create_with_one(self):
        """creates a DataMatrixCollection with one matrix"""
        matrix = DataMatrix(2, 3, ["row0", "row1"], ["col0", "col1", "col2"])
        coll = DataMatrixCollection([matrix])
        self.assertEquals(["row0", "row1"], coll.unique_row_names)
        self.assertEquals(["col0", "col1", "col2"],
                          coll.unique_column_names)
        self.assertEquals(2, coll.num_unique_rows())
        self.assertEquals(3, coll.num_unique_columns())


class MockDelimitedFile:  # pylint: disable-msg=R0903
    """Mock DelimitedFile"""

    def __init__(self, header, lines):
        """create a mock instance"""
        self.header = header
        self.lines = lines


class DataMatrixFactoryTest(unittest.TestCase):  # pylint: disable-msg=R0904
    """Test class for DataMatrixFactory"""

    def setUp(self):  # pylint: disable-msg=C0103
        """text fixture"""
        self.dfile = MockDelimitedFile(["H1", "H2", "H3"],
                                       [["R1", 1, 2], ["R2", 3, 4]])
        self.dfile_with_na = MockDelimitedFile(["H1", "H2", "H3"],
                                               [["R1", 'NA', 2],
                                                ["R2", 'NA', 4]])

    def test_no_filters(self):
        """test a factory without filters"""
        factory = DataMatrixFactory([])
        matrix = factory.create_from(self.dfile)
        self.assertEquals(2, matrix.num_rows())
        self.assertEquals(2, matrix.num_columns())
        self.assertEquals(["H2", "H3"], matrix.column_names())
        self.assertEquals(["R1", "R2"], matrix.row_names())
        self.assertEquals([1, 2], matrix.values()[0])
        self.assertEquals([3, 4], matrix.values()[1])

    def test_with_na_values(self):
        """test a factory with a DelimitedFile containing NA values"""
        factory = DataMatrixFactory([])
        matrix = factory.create_from(self.dfile_with_na)
        self.assertEquals(2, matrix.num_rows())
        self.assertEquals(2, matrix.num_columns())
        self.assertEquals(["H2", "H3"], matrix.column_names())
        self.assertEquals(["R1", "R2"], matrix.row_names())
        self.assertEquals([None, 2], matrix.values()[0])
        self.assertEquals([None, 4], matrix.values()[1])

    def test_simple_filter(self):
        """test a factory using a single filter"""
        factory = DataMatrixFactory([times2])
        matrix = factory.create_from(self.dfile)
        self.assertEquals(2, matrix.num_rows())
        self.assertEquals(2, matrix.num_columns())
        self.assertEquals(["H2", "H3"], matrix.column_names())
        self.assertEquals(["R1", "R2"], matrix.row_names())
        self.assertEquals([2, 4], matrix.values()[0])
        self.assertEquals([6, 8], matrix.values()[1])


def times2(matrix):
    """a simple filter that multiplies all values in the matrix by 2"""
    result = deepcopy(matrix)
    for row in range(matrix.num_rows()):
        for col in range(matrix.num_columns()):
            result.set_value_at(row, col, matrix.value_at(row, col) * 2)
    return result


class NoChangeFilterTest(unittest.TestCase):  # pylint: disable-msg=R0904
    """Test class for nochange_filter"""

    def test_simple(self):
        """simplest test case: everything kept"""
        matrix = DataMatrix(2, 2, ['R1', 'R2'], ['C1', 'C2'])
        matrix.set_values([[0.24, -0.35], [-0.42, 0.42]])
        filtered = nochange_filter(matrix)
        self.assertEquals(2, filtered.num_rows())
        self.assertEquals(2, filtered.num_columns())
        self.assertEquals([[0.24, -0.35], [-0.42, 0.42]], filtered.values())

    def test_remove_row(self):
        """remove one row"""
        matrix = DataMatrix(2, 2, ['R1', 'R2'], ['C1', 'C2'])
        matrix.set_values([[0.24, -0.35], [-0.001, None]])
        filtered = nochange_filter(matrix)
        self.assertEquals(1, filtered.num_rows())
        self.assertEquals(2, filtered.num_columns())
        self.assertEquals([[0.24, -0.35]], filtered.values())

    def test_remove_column(self):
        """remove one column"""
        matrix = DataMatrix(2, 2, ['R1', 'R2'], ['C1', 'C2'])
        matrix.set_values([[0.001, -0.35], [0, 0.42]])
        filtered = nochange_filter(matrix)
        self.assertEquals(2, filtered.num_rows())
        self.assertEquals(1, filtered.num_columns())
        self.assertEquals([[-0.35], [0.42]], filtered.values())


class CenterScaleFilterTest(unittest.TestCase):  # pylint: disable-msg=R0904
    """Test class for center_median_filter"""

    def test_filter(self):
        """test the centering"""
        matrix = DataMatrix(2, 2, ['R1', 'R2'], ['C1', 'C2'])
        matrix.set_values([[2, 3], [3, 4]])
        filtered = center_scale_filter(matrix)
        self.assertAlmostEqual(-0.70710678237309499, filtered.value_at(0, 0))
        self.assertAlmostEqual(0.70710678237309499, filtered.value_at(0, 1))
        self.assertAlmostEqual(-0.70710678237309499, filtered.value_at(1, 0))
        self.assertAlmostEqual(0.70710678237309499, filtered.value_at(1, 1))
