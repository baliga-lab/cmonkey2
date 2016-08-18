# vi: sw=4 ts=4 et:
"""util.py - cMonkey utility module

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import operator
import collections
from collections import defaultdict
import math
import numpy as np
import scipy.stats

# Python2 - Python3 compatibility
try:
    from urllib2 import urlopen
except ImportError:
    from urllib.request import urlopen

try:
    xrange
except NameError:
    xrange = range


import os
import rpy2.robjects as robjects
import gzip
import shelve
import time
import logging
import multiprocessing as mp

# RSAT organism finding is an optional feature, which we can skip in case that
# the user imports all the features through own text files
import bs4

# this tuple structure holds data of a delimited file
DelimitedFile = collections.namedtuple('DelimitedFile', ['lines', 'header'])


def make_delimited_file_from_lines(lines, sep, has_header, comment, quote):
    """Creates a delimited file from a list of lines"""

    def remove_quotes(astring, quote):
        """removes the quote character from the input string if given"""
        if quote:
            return astring.replace(quote, "")
        else:
            return astring

    def next_non_comment_index(lines, comment, line_index):
        """utility method that takes a list of strings and returns the next
        line index of a non-comment line. it also skips empty lines
        """
        if comment:
            num_lines = len(lines)
            line = lines[line_index].lstrip()
            while line_index < num_lines and (len(line) == 0 or
                                              line.startswith(comment)):
                line_index += 1
                if (line_index < num_lines):
                    line = lines[line_index].lstrip()
        return line_index

    file_header = None
    file_lines = []
    line_index = next_non_comment_index(lines, comment, 0)
    if has_header:
        file_header = lines[line_index].rstrip().split(sep)
        file_header = [remove_quotes(elem, quote) for elem in file_header]
        line_index += 1

    num_lines = len(lines)
    while line_index < num_lines:
        line_index = next_non_comment_index(lines, comment, line_index)
        if line_index < num_lines:
            stripped_line = lines[line_index].rstrip('\r\n')
            if len(stripped_line) > 0:  # to catch newline at the end of file
                line = stripped_line.split(sep)
                line = [remove_quotes(elem, quote) for elem in line]
                file_lines.append(line)
            line_index += 1

    return DelimitedFile(file_lines, file_header)


def dfile_from_text(text, sep='\t', has_header=False,
                    comment=None, quote=None):
    """creates a DelimitedFile instance from a text"""
    return make_delimited_file_from_lines(text.split('\n'), sep,
                                          has_header, comment, quote)


def read_dfile(filepath, sep='\t', has_header=False, comment=None,
               quote=None):
    """Creates the reader object"""
    lines = None
    if filepath.endswith('.gz'):
        with gzip.open(filepath) as inputfile:
            lines = [line.decode('utf-8') for line in inputfile.readlines()]
    else:
        with open(filepath, 'r') as inputfile:
            lines = inputfile.readlines()
    return make_delimited_file_from_lines(lines, sep, has_header,
                                          comment, quote)


def make_dfile_map(dfile, key_column, value_column):
    return collections.defaultdict(lambda : None,
                                   [(line[key_column], line[value_column])
                                    for line in dfile.lines])


def levenshtein_distance(str1, str2):
    """computes the Levenshtein distance. This is used in order
    to make approximate string comparisons"""
    strlen1 = len(str1)
    strlen2 = len(str2)

    dist = [[0 for _ in xrange(strlen2 + 1)] for _ in xrange(strlen1 + 1)]
    for row in xrange(strlen1 + 1):
        dist[row][0] = row
    for col in xrange(strlen2 + 1):
        dist[0][col] = col

    for col in xrange(strlen2):
        for row in xrange(strlen1):
            if str1[row] == str2[col]:
                dist[row + 1][col + 1] = dist[row][col]
            else:
                dist[row + 1][col + 1] = min(dist[row][col + 1] + 1,
                                             dist[row + 1][col] + 1,
                                             dist[row][col] + 1)
    return dist[strlen1][strlen2]


# A hyperlink with a Levenshtein score and BeautifulSoup anchor tag
RankedAnchor = collections.namedtuple('RankedAnchor', ['score', 'anchor'])


def best_matching_links(search_string, html):
    """given a search string and an HTML text, extract the best matching
    href"""
    try:
        soup = bs4.BeautifulSoup(html, "lxml")
    except:
        # this is a fallback for ancient sytems like CentOS
        soup = bs4.BeautifulSoup(html)

    links = []
    for anchor in soup.findAll('a'):
        score = levenshtein_distance(search_string, anchor['href'])
        links.append(RankedAnchor(score, anchor))
    result = sorted(links, key=operator.attrgetter('score'))
    num_results = len(result)
    if num_results > 0:
        best_score = result[0].score
        num_best = 1
        while num_best < num_results and result[num_best].score == best_score:
            num_best += 1
    return [entry.anchor['href'] for entry in result[0:num_best]]


def quantile(values, probability):
    """does the same as R's quantile function.
    values a list of numeric values
    probability a value in the range between 0 and 1
    """
    values = np.array(values)
    values = values[np.isfinite(values)]
    if len(values):
        return scipy.stats.scoreatpercentile(values, probability * 100)
    else:
        return np.nan


def r_stddev(values):
    """This is a standard deviation function, adjusted so it will
    return approximately the same value as R's sd() function would"""
    values = np.array(values)
    masked = values[np.isfinite(values)]
    num_values = len(masked)
    return round(np.std(masked) / np.sqrt(float(num_values - 1) /
                                          float(num_values)), 8)


def r_variance_columns(matrix):
    """computes the variance over the columns of a matrix, applying
    a bias of (n/n-1) over the results to match with R"""
    num_rows = len(matrix)
    bias = float(num_rows) / float(num_rows - 1)
    masked = np.ma.masked_array(matrix, np.isnan(matrix))
    result = np.var(masked, 0)
    return [(value * bias) for value in result]


def max_row_var(matrix):
    """computes the maximum row variance of a matrix"""
    masked = np.ma.masked_array(matrix, np.isnan(matrix))
    return np.mean(np.var(masked, 1, ddof=1))


def r_outer(x, y, f):
    """emulates the R "outer" function, calculating the outer product
    with a user-defined function"""
    x = np.array(x)
    y = np.array(y)
    return f(x[:, np.newaxis], y)


def mean(nparray):
    """computes the mean of a numpy array, ignoring NaN values"""
    return np.mean(np.ma.masked_array(nparray, np.isnan(nparray)))


def median(values):
    """computes the mean of a numpy array, ignoring NaN values"""
    values = np.array(values)
    values = values[np.isfinite(values)]
    return np.median(values)


def column_means(matrix):
    """computes the column means of a matrix"""
    return np.ma.filled(np.mean(np.ma.masked_array(matrix, np.isnan(matrix)),
                                axis=0), np.nan)


def row_means(matrix):
    """computes the row means of a matrix"""
    return np.ma.filled(np.mean(np.ma.masked_array(matrix, np.isnan(matrix)),
                                axis=1), np.nan)


class DocumentNotFound(Exception):
    """An exception indicating that the requested document does not exist"""
    pass


def read_url(url):
    """convenience method to read a document from a URL using the
    CMonkeyURLopener"""
    return urlopen(url).read()


def read_url_cached(url, cache_filename):
    """convenience method to read a document from a URL using the
    CMonkeyURLopener, cached version"""
    if not os.path.exists(cache_filename):
        cache_dir = os.path.dirname(cache_filename)
        if not os.path.isdir(cache_dir):
            os.makedirs(cache_dir)
        with open(cache_filename, 'wb') as outfile:
            outfile.write(read_url(url))
    with open(cache_filename, 'rb') as cached_file:
        return cached_file.read()


def get_url_cached(url, cache_filename):
    """convenience method to read a document from a URL using the
    CMonkeyURLopener, cached version, the file is only downloaded"""
    if not os.path.exists(cache_filename):
        with open(cache_filename, 'wb') as outfile:
            outfile.write(read_url(url))


class ThesaurusBasedMap:  # pylint: disable-msg=R0903
    """wrapping a thesaurus and a feature id based map for a flexible
    lookup container that can use any valid gene alias"""

    def __init__(self, synonyms, wrapped_dict):
        """create new instance"""
        self.__thesaurus = synonyms
        self.__wrapped_dict = wrapped_dict

    def __getitem__(self, key):
        """override the __getitem__ method for dictionary-like behaviour"""
        return self.__wrapped_dict[self.__thesaurus[key]]

    def __repr__(self):
        return repr(self.__wrapped_dict)

    def keys(self):
        """Returns the keys of the thesaurus"""
        return self.__thesaurus.keys()


def order2string(order):
    """returns the string representation for an order, e.g. 1st, 2nd etc."""
    if order % 10 == 1 and order != 11:
        return "%dst" % order
    elif order % 10 == 2 and order != 12:
        return "%dnd" % order
    elif order % 10 == 3 and order != 13:
        return "%drd" % order
    else:
        return "%dth" % order


def kcombinations(alist, k):
    """returns all k-combinations of the elements in alist"""
    if k == 0:
        return []
    if k == 1:
        return [[elem] for elem in alist]
    if k == len(alist):
        return [alist]

    ss1 = kcombinations(alist[1:], k - 1)
    ss1 = [[alist[0]] + s for s in ss1]
    ss2 = kcombinations(alist[1:], k)
    return ss1 + ss2


def trim_mean(values, trim):
    """returns the trim mean"""
    if not values or len(values) == 0:
        return 0

    values = sorted(values, reverse=True)
    if trim == 0.5:
        return np.median(values)

    k_floor = int(math.floor(trim * len(values)))
    trim_values_floor = values[k_floor:len(values) - k_floor]
    return np.mean(np.ma.masked_array(trim_values_floor,
                                      np.isnan(trim_values_floor)))


######################################################################
### RPY2 abstraction
######################################################################
def density(kvalues, cluster_values, bandwidth, dmin, dmax):
    """generic function to compute density scores"""
    kwargs = {'bw': bandwidth, 'adjust': 2, 'from': dmin,
              'to': dmax, 'n': 256, 'na.rm': True}
    rdens = robjects.r("""
      rdens <- function(cluster_values, kvalues, ...) {
        d <- density(cluster_values, ...);
        p <- approx(d$x, rev(cumsum(rev(d$y))), kvalues)$y
        p / sum(p, na.rm=T)
      }""")
    return rdens(robjects.FloatVector(cluster_values),
                 robjects.FloatVector(kvalues), **kwargs)


def r_set_seed(value):
    """calls R's set.seed()"""
    set_seed = robjects.r['set.seed']
    set_seed(value)


def r_runif(value):
    """calls R's set.seed()"""
    runif = robjects.r['runif']
    return runif(value)


def rnorm(num_values, std_deviation):
    """returns the result of R's rnorm function"""
    r_rnorm = robjects.r['rnorm']
    kwargs = {'sd': std_deviation}
    return r_rnorm(num_values, **kwargs)


def phyper(q, m, n, k, lower_tail=False):
    """calls the R function phyper"""
    r_phyper = robjects.r['phyper']
    kwargs = {'lower.tail': lower_tail}
    return r_phyper(robjects.FloatVector(q),
                    robjects.FloatVector(m),
                    robjects.FloatVector(n),
                    robjects.FloatVector(k), **kwargs)


def rrank(values):
    """invokes the R function rank"""
    r_rank = robjects.r['rank']
    kwargs = {'ties': 'min', 'na': 'keep'}
    return r_rank(robjects.FloatVector(values), **kwargs)


def mad(values):
    """invokes the R function mad"""
    r_mad = robjects.r['mad']
    kwargs = {'na.rm': False}
    return r_mad(robjects.FloatVector(values), **kwargs)


def sd_rnorm(values, num_rnorm_values, fuzzy_coeff):
    """computes standard deviation on values and then calls rnorm to
    generate the num_rnorm_values. This combines stddev and rnorm
    in one function for reducing rpy2 call overhead"""
    func = robjects.r("""
      sd_rnorm <- function(values, num_out_values, fuzzy_coeff) {
        sdval <- sd(values, na.rm=T) * fuzzy_coeff
        rnorm(num_out_values, sd=sdval)
      }
    """)
    return func(robjects.FloatVector(values), num_rnorm_values,
                fuzzy_coeff)


def rrank_matrix(npmatrix):
    func = robjects.r("""
      rank_mat <- function(values, nrow, ncol) {
        xr <- t(matrix(values, nrow=nrow, ncol=ncol, byrow=T))
        return (rank(xr, ties='min', na='keep') - 1)
      }
    """)
    num_rows, num_cols = npmatrix.shape
    xvec = robjects.FloatVector(npmatrix.ravel())
    res = func(xvec, num_rows, num_cols)
    # Converting the result to an NumPy in array results in a
    # surprisingly nice speedup
    return np.array(res, dtype=np.int32)


def order_fast(values, result_size, reverse=True):
    ranked = zip(values, xrange(1, len(values) + 1))
    ranked.sort(key=operator.itemgetter(0), reverse=reverse)
    return [ranked[i][1] for i in xrange(result_size)]


def rorder(values, result_size):
    """call the R version of order"""
    r_order = robjects.r['order']
    kwargs = {'decreasing': True}
    res = r_order(robjects.FloatVector(values), **kwargs)
    return res[:result_size]


def get_rvec_fun(rvecstr):
    """make scaling function based on an R vector expression string"""
    def scale(iteration):
        rvec = robjects.r(rvecstr)
        if iteration > len(rvec):
            return rvec[-1]
        else:
            return rvec[iteration - 1]
    return scale


def get_iter_fun(params, prefix, num_iterations):
    """returns an iteration function for the given prefix from the configuration parameters"""
    try:
        constval = params[prefix + '_const']
        return lambda i: constval
    except:
        pass
    try:
        rvec = params[prefix + '_rvec']
        return get_rvec_fun(rvec.replace('num_iterations', str(num_iterations)))
    except:
        raise Exception("no rvec found for prefix '%s'" % prefix)

######################################################################
### Misc functionality
######################################################################

class open_shelf:
    """A shelf content manager, so the user does not have to care about
    closing the shelf"""

    def __init__(self, filename):
        """create an instance"""
        self.__filename = filename
        self.__shelf = None

    def __enter__(self):
        """entering the manager"""
        self.__shelf = shelve.open(self.__filename)
        return self.__shelf

    def __exit__(self, type, value, tb):
        """exiting the manager"""
        self.__shelf.sync()
        self.__shelf.close()


def current_millis():
    """returns the current time in milliseconds"""
    return int(math.floor(time.time() * 1000))


def which_multiple(elems):
    result = defaultdict(int)
    for elem in elems:
        result[elem] += 1
    return {elem for elem, count in result.items() if count > 1}


class get_mp_pool:
    """pool manager"""
    def __init__(self, config_params={}):
        """use the configuration to return a pool with user-defined number of cores
        if possible"""
        if 'num_cores' in config_params:
            self.pool = mp.Pool(config_params['num_cores'])
        else:
            self.pool = mp.Pool()

    def __enter__(self):
        return self.pool

    def __exit__(self, type, value, tb):
        self.pool.close()
        self.pool.join()

__all__ = ['DelimitedFile', 'best_matching_links', 'quantile',
           'DocumentNotFound', 'CMonkeyURLopener', 'read_url',
           'read_url_cached', 'ThesaurusBasedMap', 'trim_mean']
