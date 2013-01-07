"""datamatrix_test.py - test classes for datamatrix module

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import unittest
import copy
import datamatrix as dm
import numpy as np


class DataMatrixTest(unittest.TestCase):  # pylint: disable-msg=R0904
    """Test class for DataMatrix"""

    def test_create_with_0_row_size(self):
        """create DataMatrix with a 0 row size"""
        matrix = dm.DataMatrix(0, 3)
        self.assertEquals(0, matrix.num_rows)
        self.assertEquals(0, matrix.num_columns)

    def test_create_without_names(self):
        """create DataMatrix without row and column names"""
        matrix = dm.DataMatrix(3, 4)
        self.assertEquals(3, matrix.num_rows)
        self.assertEquals(4, matrix.num_columns)
        self.assertEquals(0.0, matrix.values[0][0])
        self.assertEquals("Row 0", matrix.row_names[0])
        self.assertEquals("Row 1", matrix.row_names[1])
        self.assertEquals("Col 0", matrix.column_names[0])
        self.assertEquals("Col 1", matrix.column_names[1])

    def test_create_with_names(self):
        """create DataMatrix with row and column names"""
        matrix = dm.DataMatrix(3, 2, ["MyRow1", "MyRow2", "MyRow3"],
                               ["MyCol1", "MyCol2"])
        self.assertEquals(3, matrix.num_rows)
        self.assertEquals(2, matrix.num_columns)
        self.assertEquals(0.0, matrix.values[0][0])
        self.assertEquals("MyRow1", matrix.row_names[0])
        self.assertEquals("MyRow2", matrix.row_names[1])
        self.assertEquals("MyCol1", matrix.column_names[0])
        self.assertEquals("MyCol2", matrix.column_names[1])
        self.assertIsNotNone(str(matrix))

    def test_create_with_wrong_row_name_count(self):
        """create DataMatrix, providing the wrong number of row names"""
        self.assertRaises(ValueError, dm.DataMatrix,
                          3, 2, row_names=["MyRow1", "MyRow2"])

    def test_create_with_wrong_column_name_count(self):
        """create DataMatrix, providing the wrong number of column names"""
        self.assertRaises(ValueError, dm.DataMatrix,
                          3, 2, col_names=["MyCol1"])

    def test_create_with_init_value(self):
        """create DataMatrix with an initialization value"""
        matrix = dm.DataMatrix(2, 2, init_value=42.0)
        self.assertTrue((matrix.values == [[42.0, 42.0], [42.0, 42.0]]).all())

    def test_init_with_values_none(self):
        """invoke set_values() with None should result in ValueError"""
        matrix = dm.DataMatrix(2, 2)
        self.assertTrue((matrix.values == [[0, 0], [0, 0]]).all())

    def test_init_with_values_wrong_row_count(self):
        """invoke set_values() with wrong row count should result in
        ValueError"""
        self.assertRaises(ValueError, dm.DataMatrix, 2, 2, values=[[1, 2]])

    def test_init_with_values_wrong_column_count(self):
        """invoke set_values() with wrong row count should result in
        ValueError"""
        self.assertRaises(ValueError, dm.DataMatrix, 2, 2,
                          values=[[1, 2], [2, 3, 4]])

    def test_init_with_values_ok(self):
        """invoke set_values() with wrong row count should result in
        ValueError"""
        matrix = dm.DataMatrix(2, 2, values=[[1, 2], [2, 3]])
        self.assertTrue((matrix.values == [[1, 2], [2, 3]]).all())

    def test_submatrix_by_name_rows_only(self):
        """test creating sub matrices by row/column names"""
        matrix = dm.DataMatrix(4, 4,
                               row_names=['R0', 'R1', 'R2', 'R3'],
                               col_names=['C0', 'C1', 'C2', 'C3'],
                               values=[[1, 2, 3, 4],
                                       [4, 5, 6, 7],
                                       [8, 9, 10, 11],
                                       [12, 13, 14, 15]])
        submatrix = matrix.submatrix_by_name(row_names=['R0', 'R2'])
        self.assertEquals(submatrix.row_names, ['R0', 'R2'])
        self.assertTrue((submatrix.values == [[1, 2, 3, 4],
                                              [8, 9, 10, 11]]).all())

    def test_submatrix_by_name_columns_only(self):
        """test creating sub matrices by row/column names"""
        matrix = dm.DataMatrix(4, 4,
                               row_names=['R0', 'R1', 'R2', 'R3'],
                               col_names=['C0', 'C1', 'C2', 'C3'],
                               values=[[1, 2, 3, 4],
                                       [4, 5, 6, 7],
                                       [8, 9, 10, 11],
                                       [12, 13, 14, 15]])
        submatrix = matrix.submatrix_by_name(row_names=None,
                                             column_names=['C1', 'C3'])
        self.assertEquals(submatrix.column_names, ['C1', 'C3'])
        self.assertTrue((submatrix.values == [[2, 4],
                                              [5, 7],
                                              [9, 11],
                                              [13, 15]]).all())

    def test_submatrix_by_name_rows_and_cols(self):
        """test creating sub matrices by row/column name selection"""
        matrix = dm.DataMatrix(4, 4,
                               row_names=['R0', 'R1', 'R2', 'R3'],
                               col_names=['C0', 'C1', 'C2', 'C3'],
                               values=[[1, 2, 3, 4],
                                       [4, 5, 6, 7],
                                       [8, 9, 10, 11],
                                       [12, 13, 14, 15]])
        submatrix = matrix.submatrix_by_name(row_names=['R0', 'R2'],
                                             column_names=['C1', 'C3'])
        self.assertEquals(submatrix.row_names, ['R0', 'R2'])
        self.assertEquals(submatrix.column_names, ['C1', 'C3'])
        self.assertTrue((submatrix.values == [[2, 4],
                                              [9, 11]]).all())

    def test_submatrix_by_name_rows_and_cols_with_nonexisting(self):
        """test creating sub matrices by row/column name selection
        using non-existing names"""
        matrix = dm.DataMatrix(4, 4,
                               row_names=['R0', 'R1', 'R2', 'R3'],
                               col_names=['C0', 'C1', 'C2', 'C3'],
                               values=[[1, 2, 3, 4],
                                       [4, 5, 6, 7],
                                       [8, 9, 10, 11],
                                       [12, 13, 14, 15]])
        submatrix = matrix.submatrix_by_name(row_names=['R0', 'R2', 'R5'],
                                             column_names=['C1', 'C3', 'C5'])
        self.assertEquals(submatrix.row_names, ['R0', 'R2'])
        self.assertEquals(submatrix.column_names, ['C1', 'C3'])
        self.assertTrue((submatrix.values == [[2, 4],
                                              [9, 11]]).all())

    def test_submatrix_by_rows(self):
        """test creating sub matrices by providing row indexes"""
        matrix = dm.DataMatrix(4, 2,
                               row_names=['R0', 'R1', 'R2', 'R3'],
                               col_names=['C0', 'C1'],
                               values=[[1, 2],
                                       [3, 4],
                                       [5, 6],
                                       [7, 8]])
        submatrix = matrix.submatrix_by_rows([1, 3])
        self.assertEquals(submatrix.row_names, ['R1', 'R3'])
        self.assertEquals(submatrix.column_names, ['C0', 'C1'])
        self.assertTrue((submatrix.values == [[3, 4], [7, 8]]).all())

    def test_sorted_by_rowname(self):
        matrix = dm.DataMatrix(3, 3,
                               row_names=['R0', 'R2', 'R1'],
                               col_names=['C0', 'C1', 'C2'],
                               values=[[1, 2, 3],
                                       [4, 5, 6],
                                       [8, 9, 10]])
        sorted_matrix = matrix.sorted_by_row_name()
        self.assertEquals(sorted_matrix.row_names, ['R0', 'R1', 'R2'])
        self.assertTrue((sorted_matrix.values == [[1, 2, 3],
                                                  [8, 9, 10],
                                                  [4, 5, 6]]).all())

    def test_sorted_by_rowname_duplicate_row_names(self):
        matrix = dm.DataMatrix(4, 3,
                               row_names=['R0', 'R2', 'R1', 'R1'],
                               col_names=['C0', 'C1', 'C2'],
                               values=[[1, 2, 3],
                                       [4, 5, 6],
                                       [8, 9, 10],
                                       [11, 12, 13]])
        sorted_matrix = matrix.sorted_by_row_name()
        self.assertEquals(sorted_matrix.row_names, ['R0', 'R1', 'R1', 'R2'])
        self.assertTrue((sorted_matrix.values == [[1, 2, 3],
                                                  [8, 9, 10],
                                                  [11, 12, 13],
                                                  [4, 5, 6]]).all())

    def test_multiply_column_by(self):
        """tests the multiply_column_by method"""
        matrix = dm.DataMatrix(2, 2,
                               row_names=['R0', 'R1'],
                               col_names=['C0', 'C1'],
                               values=[[1, 2],
                                       [3, 4]])
        multiplied = matrix.multiply_column_by(1, 2)
        self.assertEquals(multiplied.row_names, ['R0', 'R1'])
        self.assertEquals(multiplied.column_names, ['C0', 'C1'])
        self.assertEquals(matrix, multiplied)
        self.assertTrue((multiplied.values == [[1, 4], [3, 8]]).all())


    def test_max(self):
        """tests the max() method"""
        matrix = dm.DataMatrix(2, 2,
                                row_names=['R0', 'R1'],
                                col_names=['C0', 'C1'],
                                values=[[1, np.inf],
                                        [np.nan, 4]])
        self.assertEquals(4, matrix.max())

    def test_min(self):
        """tests the min() method"""
        matrix = dm.DataMatrix(2, 2,
                                row_names=['R0', 'R1'],
                                col_names=['C0', 'C1'],
                                values=[[1, -np.inf],
                                        [np.nan, 4]])
        self.assertEquals(1, matrix.min())

    def test_row_values(self):
        """tests the row_values() method"""
        matrix = dm.DataMatrix(2, 3,
                               values=[[1.0, 2.0, 3.0],
                                       [4.0, 5.0, 6.0]])
        rowvals = matrix.row_values(0)
        self.assertTrue((rowvals == [1.0, 2.0, 3.0]).all())
        rowvals[0] = 42.0
        self.assertTrue((rowvals == [42.0, 2.0, 3.0]).all())
        self.assertTrue((matrix.values == [[1.0, 2.0, 3.0],
                                           [4.0, 5.0, 6.0]]).all())

    def test_column_values(self):
        """tests the column_values() method"""
        matrix = dm.DataMatrix(2, 3,
                               values=[[1.0, 2.0, 3.0],
                                       [4.0, 5.0, 6.0]])
        colvals = matrix.column_values(1)
        self.assertTrue((colvals == [2.0, 5.0]).all())
        colvals[1] = 42.0
        self.assertTrue((colvals == [2.0, 42.0]).all())
        self.assertTrue((matrix.values == [[1.0, 2.0, 3.0],
                                           [4.0, 5.0, 6.0]]).all())

    def test_residual(self):
        """tests the residual() method"""
        matrix = dm.DataMatrix(2, 2,
                               values=[[4.0, 7.0],
                                       [5.0, 8.0]])
        self.assertEquals(0, matrix.residual())

    def test_residual2(self):
        """tests the residual() method"""
        matrix = dm.DataMatrix(3, 3,
                               values=[[1000, -4000, 7000],
                                       [-2000, 5000, -8000],
                                       [3000, -6000, 9000]])
        self.assertAlmostEqual(4049.38271604938, matrix.residual())

    def test_residual_var_normalize(self):
        """tests the residual() method. Note that this method
        seems to make rounding errors in the 5th place"""
        matrix = dm.DataMatrix(3, 3,
                               values=[[1000, -4000, 7000],
                                       [-2000, 5000, -8000],
                                       [3000, -6000, 9000]])
        max_row_var = matrix.row_variance()
        self.assertAlmostEqual(0.000105128205128205,
                               matrix.residual(max_row_variance=max_row_var), places=4)

    def test_fix_extreme_values(self):
        """tests the adjustment function"""
        matrix = dm.DataMatrix(3, 2,
                               row_names=['R0', 'R1', 'R3'],
                               col_names=['C0', 'C1'],
                               values=[[-1.01, np.nan],
                                       [np.inf, -22.0],
                                       [-19.9, -25.3]])
        matrix.fix_extreme_values()
        self.assertTrue((matrix.values == [[-1.01, -1.01],
                                           [-1.01, -19.9],
                                           [-19.9, -19.9]]).all())
        

class DataMatrixCollectionTest(unittest.TestCase):  # pylint: disable-msg=R0904
    """Test class for MatrixCollection"""

    def test_create_with_one(self):
        """creates a DataMatrixCollection with one matrix"""
        matrix = dm.DataMatrix(2, 3, ["row0", "row1"], ["col0", "col1", "col2"])
        coll = dm.DataMatrixCollection([matrix])
        self.assertEquals(["row0", "row1"], coll.unique_row_names())
        self.assertEquals(["col0", "col1", "col2"],
                          coll.unique_column_names())
        self.assertEquals(2, coll.num_unique_rows())
        self.assertEquals(3, coll.num_unique_columns())


class MockDelimitedFile:  # pylint: disable-msg=R0903
    """Mock DelimitedFile"""

    def __init__(self, header, lines):
        """create a mock instance"""
        self.__header = header
        self.__lines = lines

    def header(self):
        """returns the header"""
        return self.__header

    def lines(self):
        """returns thes lines"""
        return self.__lines


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
        factory = dm.DataMatrixFactory([])
        matrix = factory.create_from(self.dfile)
        self.assertEquals(2, matrix.num_rows)
        self.assertEquals(2, matrix.num_columns)
        self.assertEquals(matrix.column_names, ["H2", "H3"])
        self.assertEquals(matrix.row_names, ["R1", "R2"])
        self.assertTrue((matrix.values[0] == [1, 2]).all())
        self.assertTrue((matrix.values[1] == [3, 4]).all())

    def test_with_na_values(self):
        """test a factory with a DelimitedFile containing NA values"""
        factory = dm.DataMatrixFactory([])
        matrix = factory.create_from(self.dfile_with_na)
        self.assertEquals(2, matrix.num_rows)
        self.assertEquals(2, matrix.num_columns)
        self.assertEquals(matrix.column_names, ["H2", "H3"])
        self.assertEquals(matrix.row_names, ["R1", "R2"])
        self.assertTrue(np.isnan(matrix.values[0][0]))
        self.assertEquals(2.0, matrix.values[0][1])
        self.assertTrue(np.isnan(matrix.values[1][0]))
        self.assertEquals(4.0, matrix.values[1][1])

    def test_simple_filter(self):
        """test a factory using a single filter"""
        factory = dm.DataMatrixFactory([times2])
        matrix = factory.create_from(self.dfile)
        self.assertEquals(2, matrix.num_rows)
        self.assertEquals(2, matrix.num_columns)
        self.assertEquals(matrix.column_names, ["H2", "H3"])
        self.assertEquals(matrix.row_names, ["R1", "R2"])
        self.assertTrue((matrix.values[0] == [2, 4]).all())
        self.assertTrue((matrix.values[1] == [6, 8]).all())


def times2(matrix):
    """a simple filter that multiplies all values in the matrix by 2"""
    result = copy.deepcopy(matrix)
    for row in range(matrix.num_rows):
        for col in range(matrix.num_columns):
            result.values[row][col] = matrix.values[row][col] * 2
    return result


class NoChangeFilterTest(unittest.TestCase):  # pylint: disable-msg=R0904
    """Test class for nochange_filter"""

    def test_simple(self):
        """simplest test case: everything kept"""
        matrix = dm.DataMatrix(2, 2, ['R1', 'R2'], ['C1', 'C2'],
                               values=[[0.24, -0.35], [-0.42, 0.42]])
        filtered = dm.nochange_filter(matrix)
        self.assertEquals(2, filtered.num_rows)
        self.assertEquals(2, filtered.num_columns)
        self.assertTrue((filtered.values == [[0.24, -0.35],
                                               [-0.42, 0.42]]).all())

    def test_remove_row(self):
        """remove one row"""
        matrix = dm.DataMatrix(2, 2, ['R1', 'R2'], ['C1', 'C2'],
                               values=[[0.24, -0.35], [-0.001, np.nan]])
        filtered = dm.nochange_filter(matrix)
        self.assertEquals(1, filtered.num_rows)
        self.assertEquals(2, filtered.num_columns)
        self.assertTrue((filtered.values == [[0.24, -0.35]]).all())

    def test_remove_column(self):
        """remove one column"""
        matrix = dm.DataMatrix(2, 2, ['R1', 'R2'], ['C1', 'C2'],
                               values=[[0.001, -0.35], [np.nan, 0.42]])
        filtered = dm.nochange_filter(matrix)
        self.assertEquals(2, filtered.num_rows)
        self.assertEquals(1, filtered.num_columns)
        self.assertTrue((filtered.values == [[-0.35], [0.42]]).all())


class CenterScaleFilterTest(unittest.TestCase):  # pylint: disable-msg=R0904
    """Test class for center_median_filter"""

    def test_filter(self):
        """test the centering"""
        matrix = dm.DataMatrix(2, 2, ['R1', 'R2'], ['C1', 'C2'],
                               values=[[2, 3], [3, 4]])
        filtered = dm.center_scale_filter(matrix).values
        self.assertAlmostEqual(-0.70710678237309499, filtered[0][0])
        self.assertAlmostEqual(0.70710678237309499, filtered[0][1])
        self.assertAlmostEqual(-0.70710678237309499, filtered[1][0])
        self.assertAlmostEqual(0.70710678237309499, filtered[1][1])


def as_sorted_flat_values(matrices):
    """this method is now inlined into quantile_normalize_scores
    we keep it here with its tests to demonstrate how it works
    """
    return np.transpose(np.asarray([np.sort(matrix.values.flatten())
                                    for matrix in matrices]))


class QuantileNormalizeTest(unittest.TestCase): # pylint: disable-msg=R0904
    """Test cases for quantile normalization. This is the
    central algorithm for combining lists of scoring matrices"""

    def test_as_flat_values_empty_input(self):
        """tests as_sorted_flat_values() with no input"""
        self.assertTrue((as_sorted_flat_values([]) == []).all())

    def test_as_sorted_flat_values(self):
        """tests that the flat values of the input matrices are
        all put in one big numpy array"""
        m1 = dm.DataMatrix(2, 2, values=[[2, np.nan], [3, 4]])
        m2 = dm.DataMatrix(2, 2, values=[[6, 5], [4, 3]])
        flat_values = as_sorted_flat_values([m1, m2])
        self.assertEquals(4, len(flat_values))
        self.assertTrue((flat_values[0] == [2, 3]).all())
        self.assertTrue((flat_values[1] == [3, 4]).all())
        self.assertTrue((flat_values[2] == [4, 5]).all())
        self.assertTrue(np.isnan(flat_values[3][0]))
        self.assertEquals(6, flat_values[3][1])

    def test_ranks(self):
        ranked = dm.ranks(np.array([3.0, 1.0, 2.0]))
        self.assertTrue((ranked ==  [2, 0, 1]).all())

    def test_qm_result_matrices(self):
        m1 = dm.DataMatrix(2, 2, values=[[2, 1], [3, 4]])
        m2 = dm.DataMatrix(2, 2, values=[[6, 5], [4, 3]])
        tmp_mean = np.array([1.0, 2.0, 3.0, 4.0])
        result = dm.qm_result_matrices([m1, m2], tmp_mean)
        self.assertEquals(2, len(result))
        qm1 = result[0]
        qm2 = result[1]
        self.assertTrue((qm1.values == [[2, 1], [3, 4]]).all())
        self.assertTrue((qm2.values == [[4, 3], [2, 1]]).all())

    def test_quantile_normalize_scores_with_all_defined_weights(self):
        """happy path for quantile normalization"""
        m1 = dm.DataMatrix(2, 2, values=[[1, 3], [2, 4]])
        m2 = dm.DataMatrix(2, 2, values=[[2.3, 2.5], [2.1, 2.31]])
        result = dm.quantile_normalize_scores([m1, m2], [6.0, 1.0])

        outmatrix1 = result[0].values
        self.assertAlmostEqual(0.5785714, outmatrix1[0][0])
        self.assertAlmostEqual(1.45071428, outmatrix1[0][1])
        self.assertAlmostEqual(1.02142857, outmatrix1[1][0])
        self.assertAlmostEqual(1.89285714, outmatrix1[1][1])

        outmatrix2 = result[1].values
        self.assertAlmostEqual(1.02142857, outmatrix2[0][0])
        self.assertAlmostEqual(1.89285714, outmatrix2[0][1])
        self.assertAlmostEqual(0.5785714, outmatrix2[1][0])
        self.assertAlmostEqual(1.45071428, outmatrix2[1][1])

    def test_quantile_normalize_scores_with_no_weights(self):
        """no weights -> fall back to row means"""
        m1 = dm.DataMatrix(2, 2, values=[[1, 3], [2, 4]])
        m2 = dm.DataMatrix(2, 2, values=[[2.3, 2.5], [2.1, 2.31]])
        result = dm.quantile_normalize_scores([m1, m2], None)

        outmatrix1 = result[0].values
        self.assertAlmostEqual(1.55, outmatrix1[0][0])
        self.assertAlmostEqual(2.655, outmatrix1[0][1])
        self.assertAlmostEqual(2.15, outmatrix1[1][0])
        self.assertAlmostEqual(3.25, outmatrix1[1][1])

        outmatrix2 = result[1].values
        self.assertAlmostEqual(2.15, outmatrix2[0][0])
        self.assertAlmostEqual(3.25, outmatrix2[0][1])
        self.assertAlmostEqual(1.55, outmatrix2[1][0])
        self.assertAlmostEqual(2.655, outmatrix2[1][1])

    def test_quantile_normalize_scores_with_undefined_weight(self):
        """one undefined weight"""
        m1 = dm.DataMatrix(2, 2, values=[[1, 3], [2, 4]])
        m2 = dm.DataMatrix(2, 2, values=[[2.3, 2.5], [2.1, 2.31]])
        result = dm.quantile_normalize_scores([m1, m2], [6.0, np.nan])

        outmatrix1 = result[0].values
        self.assertAlmostEqual(1.0, outmatrix1[0][0])
        self.assertAlmostEqual(3.0, outmatrix1[0][1])
        self.assertAlmostEqual(2.0, outmatrix1[1][0])
        self.assertAlmostEqual(4.0, outmatrix1[1][1])

        outmatrix2 = result[1].values
        self.assertAlmostEqual(2.0, outmatrix2[0][0])
        self.assertAlmostEqual(4.0, outmatrix2[0][1])
        self.assertAlmostEqual(1.0, outmatrix2[1][0])
        self.assertAlmostEqual(3.0, outmatrix2[1][1])


if __name__ == '__main__':
    SUITE = []
    SUITE.append(unittest.TestLoader().loadTestsFromTestCase(QuantileNormalizeTest))
    SUITE.append(unittest.TestLoader().loadTestsFromTestCase(DataMatrixTest))
    SUITE.append(unittest.TestLoader().loadTestsFromTestCase(DataMatrixFactoryTest))
    SUITE.append(unittest.TestLoader().loadTestsFromTestCase(CenterScaleFilterTest))
    unittest.TextTestRunner(verbosity=2).run(unittest.TestSuite(SUITE))
