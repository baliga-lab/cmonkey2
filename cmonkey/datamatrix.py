# vi: sw=4 ts=4 et:
"""datatypes.py - data types for cMonkey

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import scipy
import numpy as np
import operator
import cmonkey.util as util
import logging
import gzip
import os
import random
import pandas

# Python2/Python3 compatibility
try:
    xrange
except NameError:
    xrange = range

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

        if row_names is None:
            self.row_names = ["Row " + str(i) for i in xrange(nrows)]
        else:
            if len(row_names) != nrows:
                raise ValueError("number of row names should be %d" % nrows)
            self.row_names = row_names

        if col_names is None:
            self.column_names = ["Col " + str(i) for i in xrange(ncols)]
        else:
            if len(col_names) != ncols:
                raise ValueError("number of column names should be %d" % ncols)
            self.column_names = col_names

        if values is not None:
            check_values()
            self.values = np.array(values, dtype=np.float64)
        else:
            self.values = np.zeros((nrows, ncols))
            if init_value is not None:
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
        if self.row_indexes is None:
            self.row_indexes = {row: index for index, row in enumerate(self.row_names)}
        return [self.row_indexes[name] if name in self.row_indexes else -1
                for name in row_names]

    def column_indexes_for(self, column_names):
        """returns the column indexes with the matching names"""
        if self.column_indexes is None:
            self.column_indexes = {col: index for index, col in enumerate(self.column_names)}
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
            if row_indexes is None and column_indexes is None:
                return self.values
            elif row_indexes is None:
                return self.values[:, column_indexes]
            elif column_indexes is None:
                return self.values[row_indexes]
            else:
                return self.values[row_indexes][:, column_indexes]

        if row_names is None:
            row_names = self.row_names
            row_indexes = None
        else:
            row_names = [name for name in sorted(row_names)
                         if name in set(self.row_names)]
            row_indexes = self.row_indexes_for(row_names)

        if column_names is None:
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
        self.values -= self.quantile(quantile)

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

    def median(self):
        """returns the mean value"""
        return util.median(self.values)

    def row_variance(self):
        if self.__row_variance is None:
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
        if max_row_variance is not None:
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

        self.values[np.isinf(self.values)] = maxval
        self.values[np.isnan(self.values)] = maxval #Should this actually be 0 or median?

        #01-28-15 reordered to make sure that NAs are removed before this test
        self.values[self.values < min_value] = np.min(masked[masked >= min_value])

    def __repr__(self):
        """returns a string representation of this matrix"""
        return str(self)

    def __str__(self):
        """returns a string representation of this matrix"""
        result = "%10s" % 'Row'
        result += ' '.join([("%10s" % name)
                            for name in self.column_names]) + '\n'
        for row_index, rowname in enumerate(self.row_names):
            result += ("%10s" % rowname) + ' '
            result += ' '.join([("%10f" % value)
                                for value in self.values[row_index]])
            result += '\n'
        return result

    def write_tsv_file(self, path, compressed=True):
        """writes this matrix to tab-separated file"""
        def write_data(outfile):
            titlerow = '\t'.join(self.column_names)
            if compressed:
                outfile.write((titlerow + '\n').encode('utf-8'))
            else:
                outfile.write(titlerow + '\n')
            for row_index, rowname in enumerate(self.row_names):
                row = [rowname]
                row.extend([('%.17f' % value) for value in self.values[row_index]])
                if compressed:
                    outfile.write(('\t'.join(row) + '\n').encode('utf-8'))
                else:
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


FILTER_THRESHOLD = 0.98
ROW_THRESHOLD = 0.17
COLUMN_THRESHOLD = 0.1


def nochange_filter(dataframe):
    """returns a new filtered DataMatrix containing only the columns and
    rows that have large enough measurements"""

    def nochange_filter_rows():
        """subfunction of nochange_filter to filter row-wise"""
        keep = []
        dmvalues = dataframe.values
        for row_index in xrange(dataframe.shape[0]):
            count = 0
            for col_index in xrange(dataframe.shape[1]):
                value = dmvalues[row_index, col_index]
                if np.isnan(value) or abs(value) <= ROW_THRESHOLD:
                    count += 1
            mean = float(count) / dataframe.shape[1]
            if mean < FILTER_THRESHOLD:
                keep.append(row_index)
        return keep

    def nochange_filter_columns():
        """subfunction of nochange_filter to filter column-wise"""
        keep = []
        dmvalues = dataframe.values
        for col_index in xrange(dataframe.shape[1]):
            count = 0
            for row_index in xrange(dataframe.shape[0]):
                value = dmvalues[row_index, col_index]
                if np.isnan(value) or abs(value) <= COLUMN_THRESHOLD:
                    count += 1
            mean = float(count) / dataframe.shape[0]
            if mean < FILTER_THRESHOLD:
                keep.append(col_index)
        return keep

    rows_to_keep = nochange_filter_rows()
    cols_to_keep = nochange_filter_columns()
    colnames = list(map(lambda col: dataframe.columns[col], cols_to_keep))
    rownames = list(map(lambda row: dataframe.index[row], rows_to_keep))
    numrows = len(rows_to_keep)
    numcols = len(cols_to_keep)

    rvalues = np.zeros((numrows, numcols))
    mvalues = dataframe.values
    for row_index in xrange(numrows):
        for col_index in xrange(numcols):
            value = mvalues[rows_to_keep[row_index], cols_to_keep[col_index]]
            rvalues[row_index, col_index] = value
    result = pandas.DataFrame(rvalues, rownames, colnames)
    return result


def row_filter(dataframe, fun):
    """generalize a matrix filter that is applying a function for each row"""
    num_rows = dataframe.shape[0]
    values = np.zeros(dataframe.shape)
    for row_index in xrange(num_rows):
        values[row_index] = fun(dataframe.values[row_index])

    return pandas.DataFrame(values, dataframe.index, dataframe.columns)


def center_scale_filter(dataframe):
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

    return row_filter(dataframe, center_scale)


def create_from_csv(csvpath, filters=[], sep='\t', quotechar='"', case_sensitive=True):
    """creates and returns an initialized, filtered DataMatrix instance"""
    if csvpath.startswith('http://'):
        raise Exception('reading from URL temporarily disabled')

    if os.path.exists(csvpath):
        df = pandas.read_csv(csvpath, index_col=0, sep=sep, quotechar=quotechar)
        df.index = map(str, df.index)
    else:
        raise Exception("File '%s' does not exist" % csvpath)
    for matrix_filter in filters:
        df = matrix_filter(df)
    if not case_sensitive:
        df.index = df.index.str.upper()
    df = df.sort_index()
    # for now, we will make a DataMatrix object from the data frame for compatibility
    return DataMatrix(df.shape[0], df.shape[1], list(df.index), list(df.columns),
                      values=df.values)


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
    if weights is not None:
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
        with util.get_mp_pool() as pool:
            results = pool.map(rank_fun,
                               [(matrix.values, matrix.row_names, matrix.column_names, tmp_mean)
                                for matrix in matrices])
        return results
    else:
        # non-parallelized
        result = []
        for i in xrange(len(matrices)):
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
def split_matrix(matrix, outdir, n, kmin, kmax):
    """Split the input matrix into n matrices with the original
    number of rows and k columns. Write the resulting matrix to
    the specified output directory"""
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    for i in xrange(n):
        k = random.randint(kmin, kmax)
        column_names = random.sample(matrix.column_names, k)
        m = matrix.submatrix_by_name(column_names=column_names)
        path = '%s/ratios-%03d.tsv' % (outdir, i + 1)
        m.write_tsv_file(path)


__all__ = ['DataMatrix', 'nochange_filter', 'center_scale_filter']
