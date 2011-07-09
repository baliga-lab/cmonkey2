"""datatypes.py - data types for cMonkey

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""


class DataMatrix:
    """
    A 2 dimensional data matrix class, with optional row and column names
    The matrix is initialized with a fixed dimension size and can not
    be resized after initialization.
    Names and values of a matrix instance can be modified
    DataMatrix does not make any assumptions about its value or name types
    """
    def __init__(self, nrows, ncols, row_names=None, column_names=None):
        """create a DataMatrix instance"""
        self.values = [[0.0 for _ in range(ncols)] for _ in range(nrows)]
        if not row_names:
            self.row_names = ["Row " + str(i) for i in range(nrows)]
        else:
            if len(row_names) != nrows:
                raise ValueError("number of row names should be %d" % nrows)
            self.row_names = row_names

        if not column_names:
            self.column_names = ["Column " + str(i) for i in range(ncols)]
        else:
            if len(column_names) != ncols:
                raise ValueError("number of column names should be %d" % ncols)
            self.column_names = column_names

    def num_rows(self):
        """returns the number of rows"""
        return len(self.values)

    def num_columns(self):
        """returns the number of columns"""
        if self.num_rows() == 0:
            return 0
        else:
            return len(self.values[0])

    def value_at(self, row, column):
        """retrieve the value at the specified position"""
        return self.values[row][column]

    def set_value_at(self, row, column, value):
        """set the value at the specified position"""
        self.values[row][column] = value

    def row_name(self, row):
        """retrieve the name for the specified row"""
        return self.row_names[row]

    def column_name(self, row):
        """retrieve the name for the specified column"""
        return self.column_names[row]


class DataMatrixCollection:
    """A collection of DataMatrix objects containing gene expression values
    It also offers functionality to combine the comtained matrices
    """
    def __init__(self, matrices):
        self.matrices = matrices
        self.unique_row_names = self.make_unique_names(
            lambda matrix: matrix.row_names)
        self.unique_column_names = self.make_unique_names(
            lambda matrix: matrix.column_names)

    def num_unique_rows(self):
        """returns the number of unique rows"""
        return len(self.unique_row_names)

    def num_unique_columns(self):
        """returns the number of unique columns"""
        return len(self.unique_column_names)

    def make_unique_names(self, name_extract_fun):
        """helper method to create a unique name list
        name_extract_fun is a function to return a list of names from
        a matrix"""
        result = []
        for matrix in self.matrices:
            names = name_extract_fun(matrix)
            for name in names:
                result.append(name)
        result.sort()
        return result

__all__ = ['DataMatrix', 'DataMatrixCollection']
