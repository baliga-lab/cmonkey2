# vi: sw=4 ts=4 et:
"""datatypes.py - data types for cMonkey

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import scipy
import numpy as np
import operator
import util
import logging
import rpy2.robjects as robj


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
            for row_index in xrange(nrows):
                inrow = values[row_index]
                if len(inrow) != ncols:
                    raise ValueError(("row %d: number of columns should be " +
                                      "%d (was %d)")
                                     % (row_index, ncols, len(inrow)))

        if row_names == None:
            self.row_names = ["Row " + str(i) for i in range(nrows)]
        else:
            if len(row_names) != nrows:
                raise ValueError("number of row names should be %d" % nrows)
            self.row_names = list(row_names)

        if col_names == None:
            self.column_names = ["Col " + str(i) for i in xrange(ncols)]
        else:
            if len(col_names) != ncols:
                raise ValueError("number of column names should be %d" % ncols)
            self.column_names = list(col_names)

        if values != None:
            check_values()
            self.values = np.array(values, dtype=np.float64)
        else:
            self.values = np.zeros((nrows, ncols))
            if init_value != None:
                self.values.fill(init_value)
        self.__row_variance = None

    def num_rows(self):
        """returns the number of rows"""
        return len(self.row_names)

    def num_columns(self):
        """returns the number of columns"""
        if self.num_rows() == 0:
            return 0
        else:
            return len(self.column_names)

    def row_indexes(self, row_names):
        """returns the row indexes with the matching names"""
        return self.__find_indexes(self.row_names, row_names)

    def column_indexes(self, column_names):
        """returns the column indexes with the matching names"""
        return self.__find_indexes(self.column_names, column_names)

    def __find_indexes(self, names, search_names):
        """generic finder method to search name indexes in a numpy array"""
        result = []
        for name in search_names:
            result.append(names.index(name))
        return result

    def flat_values(self):
        """returns all values as a single sequence"""
        return self.values.flatten()

    def row_values(self, row):
        """returns the values in the specified row"""
        return self.values[[row]][0]

    def column_values(self, column):
        """returns the values in the specified column"""
        return self.values[:, [column]].flatten()

    def __getitem__(self, row_index):
        """return the row at the specified position"""
        return self.values[row_index]

    def row_name(self, row):
        """retrieve the name for the specified row"""
        return self.row_names[row]

    def column_name(self, row):
        """retrieve the name for the specified column"""
        return self.column_names[row]

    def submatrix_by_rows(self, row_indexes):
        """extract a submatrix with the specified rows.
        row_indexes needs to be sorted"""
        new_values = self.values[[row_indexes]]
        return DataMatrix(len(row_indexes), self.num_columns(),
                          row_names=[self.row_names[index]
                                     for index in row_indexes],
                          col_names=self.column_names,
                          values=new_values)

    def submatrix_by_name(self, row_names=None, column_names=None):
        """extract a submatrix with the specified rows and columns
        Selecting by name is more common than selecting by index
        in cMonkey, because submatrices are often selected based
        on memberships.
        Note: Currently, no duplicate row names or column names are
        supported. Furthermore, the submatrices potentially share
        the original matrix's values and so, writing to the submatrix
        will change the original matrix, too. Recommended to use
        submatrices read-only
        """
        def make_values(row_indexes, column_indexes):
            """creates an array from the selected rows and columns"""
            if row_indexes == None and column_indexes == None:
                return self.values
            elif row_indexes == None:
                return self.values[:, column_indexes]
            elif column_indexes == None:
                return self.values[row_indexes]
            else:
                return self.values[row_indexes][:, column_indexes]

        if row_names == None:
            row_names = self.row_names
            row_indexes = None
        else:
            row_names = [name for name in row_names
                         if name in self.row_names]
            row_indexes = self.row_indexes(row_names)

        if column_names == None:
            column_names = self.column_names
            col_indexes = None
        else:
            column_names = [name for name in column_names
                            if name in self.column_names]
            col_indexes = self.column_indexes(column_names)

        new_values = make_values(row_indexes, col_indexes)
        return DataMatrix(len(row_names), len(column_names), row_names,
                          column_names, values=new_values)

    def sorted_by_row_name(self):
        """returns a version of this table, sorted by row name"""
        row_pairs = []
        for row_index in xrange(len(self.row_names)):
            row_pairs.append((self.row_names[row_index], row_index))
        row_pairs.sort()
        new_rows = []
        new_row_names = []
        for row_pair in row_pairs:
            new_row_names.append(row_pair[0])
            new_rows.append(self.values[row_pair[1]])
        return DataMatrix(self.num_rows(), self.num_columns(),
                          new_row_names, self.column_names,
                          values=new_rows)

    def column_means(self):
        """Returns a numpy array, containing the column means"""
        return util.column_means(self.values)

    def row_means(self):
        """Returns a numpy array, containing the column means"""
        return util.row_means(self.values)

    ######################################################################
    #### Operations on the matrix values
    ######################################################################

    def multiply_column_by(self, column, factor):
        """Mulitplies the specified column by a certain factor"""
        self.values[:, column] *= factor
        return self

    def max(self):
        """return the maximum value in this matrix"""
        return np.amax(self.values[np.isfinite(self.values)])

    def quantile(self, probability):
        """returns the result of the quantile function over all contained
        values"""
        return util.quantile(self.values.ravel(), probability)

    def min(self):
        """return the minimum value in this matrix"""
        return np.amin(self.values[np.isfinite(self.values)])

    def replace_nan_with(self, value):
        """replaces NaN with the specified value"""
        self.values[np.isnan(self.values)] = value

    def apply_log(self):
        """applies np.log to all values"""
        self.values[self.values != 0.0] = np.log(self.values[self.values != 0.0])

    def __neg__(self):
        """returns a new DataMatrix with the values in the matrix negated"""
        return DataMatrix(self.num_rows(), self.num_columns(),
                          self.row_names, self.column_names,
                          -self.values)

    def __add__(self, value):
        """adding a value to this matrix can be either a scalar or a matrix"""
        if isinstance(value, DataMatrix):
            return DataMatrix(self.num_rows(), self.num_columns(),
                              self.row_names, self.column_names,
                              self.values + value.values)
        else:
            return DataMatrix(self.num_rows(), self.num_columns(),
                              self.row_names, self.column_names,
                              self.values + value)

    def __sub__(self, value):
        """subtract a value from the matrix"""
        if value != 0.0:
            return DataMatrix(self.num_rows(), self.num_columns(),
                              self.row_names, self.column_names,
                              self.values - value)
        else:
            return self

    def __mul__(self, factor):
        """returns a new DataMatrix with the values in the matrix negated"""
        return DataMatrix(self.num_rows(), self.num_columns(),
                          self.row_names, self.column_names,
                          self.values * factor)

    def mean(self):
        """returns the mean value"""
        return util.mean(self.values)

    def row_variance(self):
        if self.__row_variance == None:
            self.__row_variance = util.max_row_var(self.values)
        return self.__row_variance

    def residual(self, max_row_variance=None):
        """computes the residual for this matrix, if max_row_variance is given,
        result is normalized by the row variance"""
        d_rows = self.row_means()
        d_cols = self.column_means()
        d_all = util.mean(d_rows)
        tmp = self.values + d_all - util.r_outer(d_rows, d_cols, operator.add)
        average = util.mean(np.abs(tmp))
        if max_row_variance != None:
            row_var = self.row_variance()
            if np.isnan(row_var) or row_var > max_row_variance:
                row_var = max_row_variance
            average = average / row_var
        return average

    def fix_extreme_values(self, min_value=-20.0):
        """replaces values < -20 with the smallest value that is >= -20
        replaces all NA/Inf values with the maximum value in the matrix
        """
        masked = self.values[np.isfinite(self.values)]
        minval = np.min(masked[masked >= min_value])
        maxval = np.max(masked)
        self.values[self.values < -20.0] = np.min(masked[masked >= -20.0])
        self.values[np.isinf(self.values)] = maxval
        self.values[np.isnan(self.values)] = maxval

    def __repr__(self):
        """returns a string representation of this matrix"""
        return str(self)

    def __str__(self):
        """returns a string representation of this matrix"""
        result = "%10s" % 'Row'
        result += ' '.join([("%10s" % name)
                            for name in self.column_names]) + '\n'
        for row_index in xrange(self.num_rows()):
            result += ("%10s" % self.row_names[row_index]) + ' '
            result += ' '.join([("%10f" % value)
                                 for value in self.values[row_index]])
            result += '\n'
        return result

    def write_tsv_file(self, path):
        """writes this matrix to tab-separated file"""
        with open(path, 'w') as outfile:
            title = ['GENE']
            title.extend(self.column_names)
            titlerow = '\t'.join(title)
            outfile.write(titlerow + '\n')
            for row_index in range(len(self.row_names)):
                row = [self.row_names[row_index]]
                row.extend([('%f' % value) for value in self.values[row_index]])
                outfile.write('\t'.join(row) + '\n')
            outfile.flush()

class DataMatrixCollection:
    """A collection of DataMatrix objects containing gene expression values
    It also offers functionality to combine the comtained matrices
    """
    def __init__(self, matrices):
        self.__matrices = matrices
        self.__unique_row_names = self.__make_unique_names(
            lambda matrix: matrix.row_names)
        self.__unique_column_names = self.__make_unique_names(
            lambda matrix: matrix.column_names)

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
        for row in xrange(nrows):
            rownames.append(lines[row][0])
        values = np.empty([nrows, ncols])
        for row in xrange(nrows):
            for col in xrange(ncols):
                strval = lines[row][col + 1]
                if strval == 'NA':
                    value = np.nan
                else:
                    value = float(strval)
                values[row][col] = value
        data_matrix = DataMatrix(nrows, ncols, rownames, colnames,
                                 values=values)

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
        for row_index in xrange(data_matrix.num_rows()):
            count = 0
            for col_index in xrange(data_matrix.num_columns()):
                value = data_matrix[row_index][col_index]
                if np.isnan(value) or abs(value) <= ROW_THRESHOLD:
                    count += 1
            mean = float(count) / data_matrix.num_columns()
            if mean < FILTER_THRESHOLD:
                keep.append(row_index)
        return keep

    def nochange_filter_columns(data_matrix):
        """subfunction of nochange_filter to filter column-wise"""
        keep = []
        for col_index in xrange(data_matrix.num_columns()):
            count = 0
            for row_index in xrange(data_matrix.num_rows()):
                value = data_matrix[row_index][col_index]
                if np.isnan(value) or abs(value) <= COLUMN_THRESHOLD:
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
    for row_index in xrange(numrows):
        for col_index in xrange(numcols):
            value = matrix[rows_to_keep[row_index]][cols_to_keep[col_index]]
            result[row_index][col_index] = value
    return result


def row_filter(matrix, fun):
    """generalize a matrix filter that is applying a function for each row"""
    values = []
    for row_index in xrange(matrix.num_rows()):
        values.append(fun(matrix[row_index]))
    result = DataMatrix(matrix.num_rows(), matrix.num_columns(),
                        matrix.row_names, matrix.column_names,
                        values=values)
    return result


def center_scale_filter(matrix):
    """center the values of each row around their median and scale
    by their standard deviation"""

    def center_scale(row):
        """centers the provided row around the median"""
        filtered = row[np.isfinite(row)]
        center = scipy.median(filtered)
        scale = util.r_stddev(filtered)
        nurow = [((value - center) / scale)
                if not np.isnan(value) else value for value in row]
        return nurow

    return row_filter(matrix, center_scale)


def quantile_normalize_scores(matrices, weights=None):
    """quantile normalize scores against each other"""

    flat_values = as_sorted_flat_values(matrices)
    #logging.info("COMPUTING WEIGHTED MEANS...")
    #start_time = util.current_millis()
    if weights != None:
        tmp_mean = weighted_row_means(flat_values, weights)
    else:
        tmp_mean = util.row_means(flat_values)
    #elapsed = util.current_millis() - start_time
    #logging.info("weighted means in %f s.", elapsed / 1000.0)
    #start_time = util.current_millis()
    result = qm_result_matrices(matrices, tmp_mean)
    #elapsed = util.current_millis() - start_time
    #logging.info("result matrices built in %f s.", elapsed / 1000.0)
    return result


def as_sorted_flat_values(matrices):
    """rearranges the scores in the input matrices into a matrix
    with |matrices| columns where the columns contain the values
    of each matrix in sorted order"""
    return np.transpose(np.asarray([np.sort(matrix.flat_values())
                                    for matrix in matrices]))


def weighted_row_means(matrix, weights):
    """compute weighted row means"""
    #start_time = util.current_millis()
    # multiply each column of matrix with each component of the
    # weight vector: Using matrix multiplication resulted in speedup
    # from 125 s. to 0.125 seconds over apply_along_axis() (1000x faster)!
    scaled = weights * matrix
    #elapsed = util.current_millis() - start_time
    #logging.info("APPLIED WEIGHTS TO COLUMNS in %f s.", elapsed / 1000.0)
    scale = np.sum(np.ma.masked_array(weights, np.isnan(weights)))
    return util.row_means(scaled) / scale


def ranks(values):
    """optimization: write a map from value to first index in
    sorted_values"""
    values = values.argsort()
    ranks = np.empty(len(values), int)
    ranks[values] = np.arange(len(values))
    return ranks


def qm_result_matrices(matrices, tmp_mean):
    """builds the resulting matrices by looking at the rank of their
    original values and retrieving the means at the specified position"""
    result = []
    for i in range(len(matrices)):
        matrix = matrices[i]
        values = matrix.values
        num_rows, num_cols = values.shape
        xvec = robj.FloatVector(values.reshape(values.size))
        xr = robj.r.t(robj.r.matrix(xvec, nrow=num_rows, ncol=num_cols, byrow=True))
        rankvals = [value - 1 for value in util.rrank(xr)]  # adjust to be 0-based
        values = np.reshape(tmp_mean[rankvals], (matrix.num_rows(),
                                                 matrix.num_columns()))
        outmatrix = DataMatrix(matrix.num_rows(),
                               matrix.num_columns(),
                               matrix.row_names,
                               matrix.column_names,
                               values=values)
        result.append(outmatrix)
    return result


__all__ = ['DataMatrix', 'DataMatrixCollection', 'DataMatrixFactory',
           'nochange_filter', 'center_scale_filter']
