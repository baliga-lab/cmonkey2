"""datatypes.py - data types for cMonkey

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import scipy
import numpy
from util import r_stddev


class DataMatrix:
    """
    A two-dimensional data matrix class, with optional row and column names
    The matrix is initialized with a fixed dimension size and can not
    be resized after initialization.
    Names and values of a matrix instance can be modified.
    The values themselves are implemented as a two-dimensional numpy array
    and returned values are all based on numpy arrays.
    """
    def __init__(self, nrows, ncols, row_names=None, col_names=None,
                 values=None):
        """create a DataMatrix instance"""
        def check_values():
            """Sets values from a two-dimensional list"""
            if len(values) != nrows:
                raise ValueError("number of rows should be %d" % nrows)
            for row_index in range(nrows):
                inrow = values[row_index]
                if len(inrow) != ncols:
                    raise ValueError("number of columns should be %d" % ncols)

        if not row_names:
            self.__row_names = ["Row " + str(i) for i in range(nrows)]
        else:
            if len(row_names) != nrows:
                raise ValueError("number of row names should be %d" % nrows)
            self.__row_names = row_names

        if not col_names:
            self.__column_names = ["Column " + str(i) for i in range(ncols)]
        else:
            if len(col_names) != ncols:
                raise ValueError("number of column names should be %d" % ncols)
            self.__column_names = col_names

        if values == None:
            self.__values = numpy.zeros((nrows, ncols))
        else:
            check_values()
            self.__values = numpy.array(values)

    def num_rows(self):
        """returns the number of rows"""
        return len(self.__row_names)

    def num_columns(self):
        """returns the number of columns"""
        if self.num_rows() == 0:
            return 0
        else:
            return len(self.__column_names)

    def row_names(self):
        """return the row names"""
        return self.__row_names

    def column_names(self):
        """return the column names"""
        return self.__column_names

    def values(self):
        """returns this matrix's values"""
        return self.__values

    def __getitem__(self, row_index):
        """return the row at the specified position"""
        return self.__values[row_index]

    def row_name(self, row):
        """retrieve the name for the specified row"""
        return self.__row_names[row]

    def column_name(self, row):
        """retrieve the name for the specified column"""
        return self.__column_names[row]

    def submatrix(self, row_indexes=None, column_indexes=None):
        """extract a submatrix with the specified rows and columns"""
        def select_names(names, indexes):
            """from an array of strings, select the ones given in indexes"""
            result = []
            for index in indexes:
                result.append(names[index])
            return result

        def make_values(row_indexes, column_indexes):
            """creates an array from the selected rows and columns"""
            result = []
            for row_index in row_indexes:
                row = []
                result.append(row)
                for col_index in column_indexes:
                    row.append(self.__values[row_index][col_index])
            return result

        if row_indexes == None:
            row_indexes = range(self.num_rows())
        if column_indexes == None:
            column_indexes = range(self.num_rows())
        new_row_names = select_names(self.__row_names, row_indexes)
        new_col_names = select_names(self.__column_names, column_indexes)
        new_values = make_values(row_indexes, column_indexes)
        return DataMatrix(len(row_indexes), len(column_indexes),
                          row_names=new_row_names,
                          col_names=new_col_names,
                          values=new_values)

    def __str__(self):
        """returns a string representation of this matrix"""
        result = 'Gene\t' + '\t'.join(self.__column_names) + '\n'
        for row_index in range(self.num_rows()):
            result += self.__row_names[row_index] + '\t'
            result += '\t'.join([str(value)
                                 for value in self.__values[row_index]])
            result += '\n'
        return result


class DataMatrixCollection:
    """A collection of DataMatrix objects containing gene expression values
    It also offers functionality to combine the comtained matrices
    """
    def __init__(self, matrices):
        self.__matrices = matrices
        self.__unique_row_names = self.__make_unique_names(
            lambda matrix: matrix.row_names())
        self.__unique_column_names = self.__make_unique_names(
            lambda matrix: matrix.column_names())

    def __make_unique_names(self, name_extract_fun):
        """helper method to create a unique name list
        name_extract_fun is a function to return a list of names from
        a matrix"""
        result = []
        for matrix in self.__matrices:
            names = name_extract_fun(matrix)
            for name in names:
                result.append(name)
        result.sort()
        return result

    def __getitem__(self, index):
        return self.__matrices[index]

    def num_unique_rows(self):
        """returns the number of unique rows"""
        return len(self.__unique_row_names)

    def num_unique_columns(self):
        """returns the number of unique columns"""
        return len(self.__unique_column_names)

    def unique_row_names(self):
        """returns the unique row names"""
        return self.__unique_row_names

    def unique_column_names(self):
        """returns the unique column names"""
        return self.__unique_column_names


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
        lines = delimited_file.lines()
        header = delimited_file.header()
        nrows = len(lines)
        ncols = len(header) - 1
        colnames = header[1:len(header)]
        rownames = []
        for row in range(nrows):
            rownames.append(lines[row][0])
        data_matrix = DataMatrix(nrows, ncols, rownames, colnames)
        for row in range(nrows):
            for col in range(ncols):
                strval = lines[row][col + 1]
                if strval == 'NA':
                    value = numpy.nan
                else:
                    value = float(strval)
                data_matrix[row][col] = value

        for matrix_filter in self.filters:
            data_matrix = matrix_filter(data_matrix)
        return data_matrix


FILTER_THRESHOLD = 0.98
ROW_THRESHOLD = 0.17
COLUMN_THRESHOLD = 0.1


def nochange_filter(matrix):
    """returns a new filtered DataMatrix containing only the columns and
    rows that have large enough measurements"""

    def nochange_filter_rows(data_matrix):
        """subfunction of nochange_filter to filter row-wise"""
        keep = []
        for row_index in range(data_matrix.num_rows()):
            count = 0
            for col_index in range(data_matrix.num_columns()):
                value = data_matrix[row_index][col_index]
                if numpy.isnan(value) or abs(value) <= ROW_THRESHOLD:
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
                value = data_matrix[row_index][col_index]
                if numpy.isnan(value) or abs(value) <= COLUMN_THRESHOLD:
                    count += 1
            mean = float(count) / data_matrix.num_rows()
            if mean < FILTER_THRESHOLD:
                keep.append(col_index)
        return keep

    rows_to_keep = nochange_filter_rows(matrix)
    cols_to_keep = nochange_filter_columns(matrix)
    colnames = [matrix.column_name(col) for col in cols_to_keep]
    rownames = [matrix.row_name(row) for row in rows_to_keep]
    numrows = len(rows_to_keep)
    numcols = len(cols_to_keep)

    result = DataMatrix(numrows, numcols, rownames, colnames)
    for row_index in range(numrows):
        for col_index in range(numcols):
            value = matrix[rows_to_keep[row_index]][cols_to_keep[col_index]]
            result[row_index][col_index] = value
    return result


def row_filter(matrix, fun):
    """generalize a matrix filter that is applying a function for each row"""
    values = []
    for row_index in range(matrix.num_rows()):
        row = [value for value in matrix[row_index]
               if not numpy.isnan(value)]
        values.append(fun(row))
    result = DataMatrix(matrix.num_rows(), matrix.num_columns(),
                        matrix.row_names(), matrix.column_names(),
                        values=values)
    return result


def center_scale_filter(matrix):
    """center the values of each row around their median and scale
    by their standard deviation"""

    def center_scale(row):
        """centers the provided row around the median"""
        center = scipy.median(row)
        scale = r_stddev(row)
        return [((value - center) / scale) for value in row]

    return row_filter(matrix, center_scale)

__all__ = ['DataMatrix', 'DataMatrixCollection', 'DataMatrixFactory',
           'nochange_filter', 'center_scale_filter']
