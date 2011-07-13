"""datatypes.py - data types for cMonkey

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
from scipy import median
from util import r_stddev


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

    def get_row_values(self, row):
        """returns the specified row values"""
        return self.values[row]

    def get_column_values(self, column):
        """returns the specified column values"""
        result = []
        for row_index in range(self.num_rows()):
            result.append(self.value_at(row_index, column))
        return result

    def value_at(self, row, column):
        """retrieve the value at the specified position"""
        return self.values[row][column]

    def set_value_at(self, row, column, value):
        """set the value at the specified position"""
        self.values[row][column] = value

    def set_values(self, values):
        """Sets values from a two-dimensional list"""
        if values == None:
            raise ValueError("values should be a two-dimensional list")
        if len(values) != self.num_rows():
            raise ValueError("number of rows should be %d" % self.num_rows())
        for row_index in range(self.num_rows()):
            inrow = values[row_index]
            if len(inrow) != self.num_columns():
                raise ValueError("number of columns should be %d" %
                                 self.num_columns())
            for col_index in range(self.num_columns()):
                self.set_value_at(row_index, col_index, inrow[col_index])

    def row_name(self, row):
        """retrieve the name for the specified row"""
        return self.row_names[row]

    def column_name(self, row):
        """retrieve the name for the specified column"""
        return self.column_names[row]

    def __str__(self):
        """returns a string representation of this matrix"""
        result = 'Gene\t' + '\t'.join(self.column_names) + '\n'
        for row_index in range(self.num_rows()):
            result += self.row_names[row_index] + '\t'
            result += '\t'.join([str(value)
                                 for value in self.get_row_values(row_index)])
            result += '\n'
        return result


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


class DataMatrixFactory:
    """Reader class for creating a DataMatrix from a delimited file,
    applying all supplied filters. Currently, the assumption is
    that all input ratio files have a header row and the first column
    denotes the gene names.
    (To be moved to filter class comments):
    There are a couple of constraints and special things to consider:
    - Row names are unique, throw out rows with names that are already
      in the matrix, first come, first in. This is to throw out the second
      probe of a gene in certain cases
    """

    def __init__(self, filters):
        """create a reader instance with the specified filters"""
        self.filters = filters

    def create_from(self, delimited_file):
        """creates and returns an initialized, filtered DataMatrix instance"""
        lines = delimited_file.lines
        header = delimited_file.header
        nrows = len(lines)
        ncols = len(header) - 1
        colnames = header[1:len(header)]
        rownames = []
        for row in range(nrows):
            rownames.append(delimited_file.lines[row][0])
        data_matrix = DataMatrix(nrows, ncols, rownames, colnames)
        for row in range(nrows):
            for col in range(ncols):
                strval = lines[row][col + 1]
                if strval == 'NA':
                    value = None
                else:
                    value = float(strval)
                data_matrix.set_value_at(row, col, value)

        for matrix_filter in self.filters:
            data_matrix = matrix_filter(data_matrix)
        return data_matrix


FILTER_THRESHOLD = 0.98
ROW_THRESHOLD = 0.17
COLUMN_THRESHOLD = 0.1


def nochange_filter(matrix):
    """returns a new filtered DataMatrix containing only the columns and
    rows that have large enough measurements"""
    rows_to_keep = nochange_filter_rows(matrix)
    cols_to_keep = nochange_filter_columns(matrix)
    colnames = [matrix.column_name(col) for col in cols_to_keep]
    rownames = [matrix.row_name(row) for row in rows_to_keep]
    numrows = len(rows_to_keep)
    numcols = len(cols_to_keep)

    result = DataMatrix(numrows, numcols, rownames, colnames)
    for row_index in range(numrows):
        for col_index in range(numcols):
            value = matrix.value_at(rows_to_keep[row_index],
                                    cols_to_keep[col_index])
            result.set_value_at(row_index, col_index, value)
    return result


def nochange_filter_rows(data_matrix):
    """subfunction of nochange_filter to filter row-wise"""
    keep = []
    for row_index in range(data_matrix.num_rows()):
        count = 0
        for col_index in range(data_matrix.num_columns()):
            value = data_matrix.value_at(row_index, col_index)
            if value == None or abs(value) <= ROW_THRESHOLD:
                count += 1
        mean = float(count) / data_matrix.num_columns()
        if mean < FILTER_THRESHOLD:
            keep.append(row_index)
    return keep


def nochange_filter_columns(data_matrix):
    """subfunction of nochange_filter to filter column-wise"""
    keep = []
    for col_index in range(data_matrix.num_columns()):
        count = 0
        for row_index in range(data_matrix.num_rows()):
            value = data_matrix.value_at(row_index, col_index)
            if value == None or abs(value) <= COLUMN_THRESHOLD:
                count += 1
        mean = float(count) / data_matrix.num_rows()
        if mean < FILTER_THRESHOLD:
            keep.append(col_index)
    return keep


def row_filter(matrix, fun):
    """generalize a matrix filter that is applying a function for each row"""
    values = []
    for row_index in range(matrix.num_rows()):
        row = [value for value in matrix.get_row_values(row_index)
               if value != None]
        values.append(fun(row))
    result = DataMatrix(matrix.num_rows(), matrix.num_columns(),
                        matrix.row_names, matrix.column_names)
    result.set_values(values)
    return result


def center_scale_filter(matrix):
    """center the values of each row around their median and scale
    by their standard deviation"""

    def center_scale(row):
        """centers the provided row around the median"""
        center = median(row)
        scale = r_stddev(row)
        return [((value - center) / scale) for value in row]

    return row_filter(matrix, center_scale)

__all__ = ['DataMatrix', 'DataMatrixCollection', 'DataMatrixFactory',
           'nochange_filter', 'center_scale_filter']
