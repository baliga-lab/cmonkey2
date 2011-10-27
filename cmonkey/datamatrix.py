"""datatypes.py - data types for cMonkey

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import scipy
import numpy
import util


class DataMatrix:
    """
    A two-dimensional data matrix class, with optional row and column names
    The matrix is initialized with a fixed dimension size and can not
    be resized after initialization.
    Names and values of a matrix instance can be modified.
    The values themselves are implemented as a two-dimensional numpy array
    and returned values are all based on numpy arrays.
    """

    # pylint: disable-msg=R0913
    def __init__(self, nrows, ncols, row_names=None, col_names=None,
                 values=None, init_value=None):
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
            self.__column_names = ["Col " + str(i) for i in range(ncols)]
        else:
            if len(col_names) != ncols:
                raise ValueError("number of column names should be %d" % ncols)
            self.__column_names = col_names

        if values != None:
            check_values()
            self.__values = numpy.array(values, dtype=numpy.float64)
        else:
            self.__values = numpy.zeros((nrows, ncols))
            if init_value != None:
                self.__values.fill(init_value)

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

    def flat_values(self):
        """returns all values as a single sequence"""
        return self.__values.flatten()

    def row_values(self, row):
        """returns the values in the specified row"""
        result = []
        for column in range(self.num_columns()):
            result.append(self.__values[row][column])
        return result

    def column_values(self, column):
        """returns the values in the specified column"""
        result = []
        for row in range(self.num_rows()):
            result.append(self.__values[row][column])
        return result

    def __getitem__(self, row_index):
        """return the row at the specified position"""
        return self.__values[row_index]

    def row_name(self, row):
        """retrieve the name for the specified row"""
        return self.__row_names[row]

    def column_name(self, row):
        """retrieve the name for the specified column"""
        return self.__column_names[row]

    def submatrix_by_rows(self, row_indexes):
        """extract a submatrix with the specified rows.
        row_indexes needs to be sorted"""
        new_values = []
        for row_index in row_indexes:
            new_values.append([value for value in self.__values[row_index]])
        return DataMatrix(len(row_indexes), self.num_columns(),
                          row_names=[self.__row_names[index]
                                     for index in row_indexes],
                          col_names=self.__column_names,
                          values=new_values)

    def submatrix_by_name(self, row_names=None, column_names=None):
        """extract a submatrix with the specified rows and columns
        Selecting by name is more common than selecting by index
        in cMonkey, because submatrices are often selected based
        on memberships.
        Note: Currently, no duplicate row names or column names are
        supported"""
        def find_indexes(search_names, my_names):
            """returns the indexes for the specified search names"""
            return [my_names.index(name) for name in search_names
                    if name in my_names]

        def make_values(row_indexes, column_indexes):
            """creates an array from the selected rows and columns"""
            result = []
            for row_index in row_indexes:
                row = []
                result.append(row)
                for col_index in column_indexes:
                    row.append(self.__values[row_index][col_index])
            return result

        if row_names == None:
            row_names = self.row_names()
        if column_names == None:
            column_names = self.column_names()
        row_indexes = find_indexes(row_names, self.row_names())
        col_indexes = find_indexes(column_names, self.column_names())
        new_values = make_values(row_indexes, col_indexes)
        return DataMatrix(len(row_indexes), len(col_indexes),
                          row_names=[name for name in row_names
                                     if name in self.__row_names],
                          col_names=[name for name in column_names
                                     if name in self.__column_names],
                          values=new_values)

    def sorted_by_row_name(self):
        """returns a version of this table, sorted by row name"""
        row_pairs = []
        for row_index in range(len(self.__row_names)):
            row_pairs.append((self.__row_names[row_index], row_index))
        row_pairs.sort()
        new_rows = []
        new_row_names = []
        for row_pair in row_pairs:
            new_row_names.append(row_pair[0])
            new_rows.append(self.__values[row_pair[1]])
        return DataMatrix(self.num_rows(), self.num_columns(),
                          new_row_names, self.column_names(),
                          values=new_rows)

    def column_means(self):
        """Returns a new DataMatrix containing the column means"""
        return DataMatrix(1, self.num_columns(), row_names=["Column Means"],
                          col_names=self.column_names(),
                          values=[util.column_means(self.values())])

    ######################################################################
    #### Operations on the matrix values
    ######################################################################

    def multiply_column_by(self, column, factor):
        """Mulitplies the specified column by a certain factor"""
        for row in range(self.num_rows()):
            self.__values[row][column] *= factor
        return self

    def max(self):
        """return the maximum value in this matrix"""
        return numpy.amax(self.values())

    def quantile(self, probability):
        """returns the result of the quantile function over all contained
        values"""
        return util.quantile(self.__values.ravel(), probability)

    def min(self):
        """return the minimum value in this matrix"""
        return numpy.amin(self.values())

    def __neg__(self):
        """returns a new DataMatrix with the values in the matrix negated"""
        return DataMatrix(self.num_rows(), self.num_columns(),
                          self.row_names(), self.column_names(),
                          -self.values())

    def __add__(self, value):
        """adding a value to this matrix can be either a scalar or a matrix"""
        if isinstance(value, DataMatrix):
            return DataMatrix(self.num_rows(), self.num_columns(),
                              self.row_names(), self.column_names(),
                              self.__values + value.values())
        else:
            return DataMatrix(self.num_rows(), self.num_columns(),
                              self.row_names(), self.column_names(),
                              self.__values + value)

    def __sub__(self, value):
        """subtract a value from the matrix"""
        if value != 0.0:
            return DataMatrix(self.num_rows(), self.num_columns(),
                              self.row_names(), self.column_names(),
                              self.__values - value)
        else:
            return self

    def __mul__(self, factor):
        """returns a new DataMatrix with the values in the matrix negated"""
        return DataMatrix(self.num_rows(), self.num_columns(),
                          self.row_names(), self.column_names(),
                          self.values() * factor)

    def __repr__(self):
        """returns a string representation of this matrix"""
        return str(self)

    def __str__(self):
        """returns a string representation of this matrix"""
        result = "%10s" % 'Row'
        result += ' '.join([("%10s" % name)
                            for name in self.__column_names]) + '\n'
        for row_index in range(self.num_rows()):
            result += ("%10s" % self.__row_names[row_index]) + ' '
            result += ' '.join([("%10f" % value)
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
        scale = util.r_stddev(row)
        return [((value - center) / scale) for value in row]

    return row_filter(matrix, center_scale)


def quantile_normalize_scores(matrices, weights=None):
    """quantile normalize scores against each other"""
    def extract_flat_values():
        """optimization: flat_values are only retrieved once"""
        result = []
        for index in range(len(matrices)):
            result.append(matrices[index].flat_values())
        return result

    def compute_tmp_mean_weighted(flat_values):
        """compute weighted row means"""
        in_values = []
        for index in range(len(matrices)):
            values = flat_values[index] * weights[index]
            values.sort()
            in_values.append(values)
        result = []
        num_values = len(in_values[0])
        scale = sum(weights)
        for row in range(len(in_values[0])):
            result.append(numpy.mean([inarray[row]
                                      for inarray in in_values]) / scale)
        return result

    def compute_tmp_mean_unweighted(flat_values):
        """compute unweighted row means"""
        in_values = []
        for index in range(len(matrices)):
            values = sorted(flat_values)
            in_values.append(values)
        result = []
        for row in range(len(in_values[0])):
            result.append(numpy.mean([inarray[row] for inarray in in_values]))
        return result

    def build_result_matrices(tmp_mean):
        """builds the resulting matrices by looking at the rank of their
        original values and retrieving the means at the specified position"""
        result = []
        for index in range(len(matrices)):
            sorted_values = sorted(flat_values[index])
            ranks = [sorted_values.index(value)
                     for value in flat_values[index]]
            outmatrix = DataMatrix(matrices[index].num_rows(),
                                   matrices[index].num_columns(),
                                   matrices[index].row_names(),
                                   matrices[index].column_names())
            for row in range(outmatrix.num_rows()):
                for col in range(outmatrix.num_columns()):
                    seqindex = row * outmatrix.num_columns() + col
                    outmatrix[row][col] = tmp_mean[ranks[seqindex]]
            result.append(outmatrix)
        return result

    flat_values = extract_flat_values()
    if weights != None:
        tmp_mean = compute_tmp_mean_weighted(flat_values)
    else:
        tmp_mean = compute_tmp_mean_unweighted(flat_values)
    return build_result_matrices(tmp_mean)


__all__ = ['DataMatrix', 'DataMatrixCollection', 'DataMatrixFactory',
           'nochange_filter', 'center_scale_filter']
