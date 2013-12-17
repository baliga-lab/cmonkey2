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
import gzip
import multiprocessing as mp
import os, os.path, random

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
            self.row_names = row_names

        if col_names == None:
            self.column_names = ["Col " + str(i) for i in xrange(ncols)]
        else:
            if len(col_names) != ncols:
                raise ValueError("number of column names should be %d" % ncols)
            self.column_names = col_names

        if values != None:
            check_values()
            self.values = np.array(values, dtype=np.float64)
        else:
            self.values = np.zeros((nrows, ncols))
            if init_value != None:
                self.values.fill(init_value)

        self.__row_variance = None
        self.row_indexes = None
        self.column_indexes = None

        # store the number of rows/columns as an attribute, so we do not have
        # to recalculate them
        self.num_rows = nrows
        self.num_columns = 0 if nrows == 0 else ncols

    def row_indexes_for(self, row_names):
        """returns the row indexes with the matching names"""
        if self.row_indexes == None:
            self.row_indexes = { row: index for index, row in enumerate(self.row_names) }
        return [self.row_indexes[name] if name in self.row_indexes else -1
                for name in row_names]

    def column_indexes_for(self, column_names):
        """returns the column indexes with the matching names"""
        if self.column_indexes == None:
            self.column_indexes = { col: index for index, col in enumerate(self.column_names) }
        return [self.column_indexes[name] if name in self.column_indexes else -1
                for name in column_names]

    def row_values(self, row):
        """returns the values in the specified row"""
        return self.values[[row]][0]

    def column_values(self, column):
        """returns the values in the specified column"""
        return self.values[:, [column]].flatten()

    def submatrix_by_rows(self, row_indexes):
        """extract a submatrix with the specified rows.
        row_indexes needs to be sorted"""
        new_values = self.values[[row_indexes]]
        return DataMatrix(len(row_indexes), self.num_columns,
                          row_names=[self.row_names[index] for index in row_indexes],
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
            row_names = [name for name in sorted(row_names)
                         if name in set(self.row_names)]
            row_indexes = self.row_indexes_for(row_names)

        if column_names == None:
            column_names = self.column_names
            col_indexes = None
        else:
            column_names = [name for name in sorted(column_names)
                            if name in set(self.column_names)]
            col_indexes = self.column_indexes_for(column_names)

        new_values = make_values(row_indexes, col_indexes)
        return DataMatrix(len(row_names), len(column_names), row_names,
                          column_names, values=new_values)

    def sorted_by_row_name(self):
        """returns a version of this table, sorted by row name"""
        row_names = self.row_names
        row_pairs = [(row_names[row], row) for row in xrange(len(self.row_names))]
        row_pairs.sort()
        new_row_names = [row_pair[0] for row_pair in row_pairs] 
        new_rows = [self.values[row_pair[1]] for row_pair in row_pairs]
        return DataMatrix(self.num_rows, self.num_columns,
                          new_row_names, self.column_names,
                          values=new_rows)

    ######################################################################
    #### Operations on the matrix values
    ######################################################################

    def multiply_column_by(self, column, factor):
        """Mulitplies the specified column by a certain factor"""
        self.values[:, column] *= factor
        return self

    def subtract_with_quantile(self, quantile):
        """subtracts this matrix's values with the specified quantile of its values"""
        self.values - self.quantile(quantile)

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
        d_rows = util.row_means(self.values)
        d_cols = util.column_means(self.values)
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
        for row_index in xrange(self.num_rows):
            result += ("%10s" % self.row_names[row_index]) + ' '
            result += ' '.join([("%10f" % value)
                                 for value in self.values[row_index]])
            result += '\n'
        return result

    def write_tsv_file(self, path, compressed=True):
        """writes this matrix to tab-separated file"""
        def write_data(outfile):
            title = ['GENE']
            title.extend(self.column_names)
            titlerow = '\t'.join(title)
            outfile.write(titlerow + '\n')
            for row_index in range(len(self.row_names)):
                row = [self.row_names[row_index]]
                row.extend([('%f' % value) for value in self.values[row_index]])
                outfile.write('\t'.join(row) + '\n')
        if compressed:
            if not path.endswith('.gz'):
                path = path + '.gz'
            with gzip.open(path, 'wb') as zipfile:
                write_data(zipfile)
        else:
            with open(path, 'w') as outfile:
                write_data(outfile)
                outfile.flush()


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

        # This handles header formats that omit the 0-column
        if len(lines) > 0 and len(lines[0]) - 1 > ncols:
            ncols = len(lines[0]) - 1
            colnames = header
        else:
            colnames = header[1:len(header)]

        rownames = []
        for row in xrange(nrows):
            rownames.append(lines[row][0])
        values = np.empty([nrows, ncols])
        for row in xrange(nrows):
            for col in xrange(ncols):
                strval = lines[row][col + 1]
                value = np.nan if len(strval) == 0 or strval == 'NA' else float(strval)
                values[row][col] = value

        data_matrix = DataMatrix(nrows, ncols, rownames, colnames,
                                 values=values)

        for matrix_filter in self.filters:
            data_matrix = matrix_filter(data_matrix)
        return data_matrix.sorted_by_row_name()


FILTER_THRESHOLD = 0.98
ROW_THRESHOLD = 0.17
COLUMN_THRESHOLD = 0.1


def nochange_filter(matrix):
    """returns a new filtered DataMatrix containing only the columns and
    rows that have large enough measurements"""

    def nochange_filter_rows(data_matrix):
        """subfunction of nochange_filter to filter row-wise"""
        keep = []
        dmvalues = data_matrix.values
        for row_index in xrange(data_matrix.num_rows):
            count = 0
            for col_index in xrange(data_matrix.num_columns):
                value = dmvalues[row_index][col_index]
                if np.isnan(value) or abs(value) <= ROW_THRESHOLD:
                    count += 1
            mean = float(count) / data_matrix.num_columns
            if mean < FILTER_THRESHOLD:
                keep.append(row_index)
        return keep

    def nochange_filter_columns(data_matrix):
        """subfunction of nochange_filter to filter column-wise"""
        keep = []
        dmvalues = data_matrix.values
        for col_index in xrange(data_matrix.num_columns):
            count = 0
            for row_index in xrange(data_matrix.num_rows):
                value = dmvalues[row_index][col_index]
                if np.isnan(value) or abs(value) <= COLUMN_THRESHOLD:
                    count += 1
            mean = float(count) / data_matrix.num_rows
            if mean < FILTER_THRESHOLD:
                keep.append(col_index)
        return keep

    rows_to_keep = nochange_filter_rows(matrix)
    cols_to_keep = nochange_filter_columns(matrix)
    colnames = map(lambda col: matrix.column_names[col], cols_to_keep)
    rownames = map(lambda row: matrix.row_names[row], rows_to_keep)
    numrows = len(rows_to_keep)
    numcols = len(cols_to_keep)

    result = DataMatrix(numrows, numcols, rownames, colnames)
    rvalues = result.values
    mvalues = matrix.values
    for row_index in xrange(numrows):
        for col_index in xrange(numcols):
            value = mvalues[rows_to_keep[row_index]][cols_to_keep[col_index]]
            rvalues[row_index][col_index] = value
    return result


def row_filter(matrix, fun):
    """generalize a matrix filter that is applying a function for each row"""
    values = []
    for row_index in xrange(matrix.num_rows):
        values.append(fun(matrix.values[row_index]))
    result = DataMatrix(matrix.num_rows, matrix.num_columns,
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

    logging.info("COMPUTING WEIGHTED MEANS...")
    start_time = util.current_millis()

    # rearranges the scores in the input matrices into a matrix
    # with |matrices| columns where the columns contain the values
    # of each matrix in sorted order
    flat_values = np.transpose(np.asarray([np.sort(matrix.values.flatten())
                                           for matrix in matrices]))

    elapsed = util.current_millis() - start_time
    logging.info("flattened/sorted score matrices in %f s.", elapsed / 1000.0)

    start_time = util.current_millis()
    if weights != None:
        # multiply each column of matrix with each component of the
        # weight vector: Using matrix multiplication resulted in speedup
        # from 125 s. to 0.125 seconds over apply_along_axis() (1000x faster)!
        scaled = weights * flat_values
        scale = np.sum(np.ma.masked_array(weights, np.isnan(weights)))
        tmp_mean = util.row_means(scaled) / scale
    else:
        tmp_mean = util.row_means(flat_values)
    elapsed = util.current_millis() - start_time
    logging.info("weighted means in %f s.", elapsed / 1000.0)
    start_time = util.current_millis()

    result = qm_result_matrices(matrices, tmp_mean)

    elapsed = util.current_millis() - start_time
    logging.info("result matrices built in %f s.", elapsed / 1000.0)
    return result


def ranks(values):
    """optimization: write a map from value to first index in
    sorted_values
    This used to be the original ranking, unfortunately, it is not close enough
    to the behaviour of the R function, but we keep it around here as a reference
    and base for potentially faster ranking.
    """
    values = values.argsort()
    ranks = np.empty(len(values), int)
    ranks[values] = np.arange(len(values))
    return ranks


def rank_fun(mat_mean):
    """ranking function that is run within Pool.map()"""
    values, row_names, column_names, tmp_mean = mat_mean
    num_rows, num_cols = values.shape
    rankvals = util.rrank_matrix(values)
    values = np.reshape(tmp_mean[rankvals], (num_rows, num_cols))
    return DataMatrix(num_rows, num_cols, row_names, column_names,
                      values=values)


def qm_result_matrices(matrices, tmp_mean, multiprocessing=True):
    """builds the resulting matrices by looking at the rank of their
    original values and retrieving the means at the specified position"""
    if multiprocessing:
        # parallelized ranking
        pool = mp.Pool()
        results = pool.map(rank_fun,
                           [(matrix.values, matrix.row_names, matrix.column_names, tmp_mean)
                            for matrix in matrices])
        pool.close()
        pool.join()
        return results
    else:
        # non-parallelized
        result = []
        for i in range(len(matrices)):
            matrix = matrices[i]
            values = matrix.values
            num_rows, num_cols = values.shape
            rankvals = util.rrank_matrix(values)
            values = np.reshape(tmp_mean[rankvals], (num_rows, num_cols))
            outmatrix = DataMatrix(num_rows,
                                   num_cols,
                                   matrix.row_names,
                                   matrix.column_names,
                                   values=values)
            result.append(outmatrix)
        return result


# Ensemble functionality
def split_matrix(matrix, outdir, n, k):
    """Split the input matrix into n matrices with the original
    number of rows and k columns. Write the resulting matrix to
    the specified output directory"""
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    column_sets = [random.sample(matrix.column_names, k)
                   for _ in range(n)]

    for i, column_names in enumerate(column_sets):
        m = matrix.submatrix_by_name(column_names=column_names)
        path = '%s/ratios-%03d.tsv' % (outdir, i)
        m.write_tsv_file(path)

def prepare_ensemble_matrix(ratiofile, outdir, n, k):
    matrix_factory = DataMatrixFactory([nochange_filter,
                                        center_scale_filter])
    if os.path.exists(ratiofile):
        infile = util.read_dfile(ratiofile, has_header=True, quote='\"')
        matrix = matrix_factory.create_from(infile)
        split_matrix(matrix, outdir, n, k)


__all__ = ['DataMatrix', 'DataMatrixFactory', 'nochange_filter', 'center_scale_filter']
