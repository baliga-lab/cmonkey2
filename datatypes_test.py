"""Test classes for datatypes module"""
import unittest
from datatypes import DataMatrix

class DataMatrixTest(unittest.TestCase):
    """Test class for DataMatrix"""

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
        matrix = DataMatrix(3, 2, ["MyRow1", "MyRow2","MyRow3"],
                            ["MyCol1", "MyCol2"])
        self.assertEquals(3, matrix.num_rows())
        self.assertEquals(2, matrix.num_columns())
        self.assertEquals(0.0, matrix.value_at(0, 0))
        self.assertEquals("MyRow1", matrix.row_name(0))
        self.assertEquals("MyRow2", matrix.row_name(1))
        self.assertEquals("MyCol1", matrix.column_name(0))
        self.assertEquals("MyCol2", matrix.column_name(1))

    def test_create_with_wrong_row_name_count(self):
        """create DataMatrix, providing the wrong number of row names"""
        self.assertRaises(ValueError, DataMatrix, 
                          3, 2, row_names = ["MyRow1", "MyRow2"])

    def test_create_with_wrong_column_name_count(self):
        """create DataMatrix, providing the wrong number of column names"""
        self.assertRaises(ValueError, DataMatrix, 
                          3, 2, column_names = ["MyCol1"])

    def test_set_value(self):
        """set a value in the matrix"""
        matrix = DataMatrix(3, 4)
        matrix.set_value_at(0, 1, 42.0)
        self.assertEquals(42.0, matrix.value_at(0, 1))
