"""util.py - cMonkey utility module

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import BeautifulSoup as bs
import operator
import math
import numpy
import scipy.stats
import urllib
import os
import rpy2.robjects as robjects
import gzip
import shelve


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
            stripped_line = lines[line_index].rstrip()
            if len(stripped_line) > 0:  # to catch newline at the end of file
                line = stripped_line.split(sep)
                line = [remove_quotes(elem, quote) for elem in line]
                file_lines.append(line)
            line_index += 1

    return DelimitedFile(file_lines, file_header)


class DelimitedFile:  # pylint: disable-msg=R0913
    """A file class to read text files that are delimited with certain
    separators. This class offers some flexibility over regular csv
    reading mechanisms by allowing for comments and storing optional
    headers.
    Create a DelimitedFile instance by calling DelimitedFile.read()."""
    def __init__(self, lines, header):
        self.__lines = lines
        self.__header = header

    @classmethod
    def create_from_text(cls, text, sep='\t', has_header=False,
                         comment=None, quote=None):
        """creates a DelimitedFile instance from a text"""
        return make_delimited_file_from_lines(text.split('\n'), sep,
                                              has_header, comment, quote)

    @classmethod
    def read(cls, filepath, sep='\t', has_header=False, comment=None,
             quote=None):
        """Creates the reader object"""
        lines = None
        if filepath.endswith('.gz'):
            with gzip.open(filepath) as inputfile:
                lines = inputfile.readlines()
        else:
            with open(filepath) as inputfile:
                lines = inputfile.readlines()
        return make_delimited_file_from_lines(lines, sep, has_header,
                                              comment, quote)

    def lines(self):
        """returns the lines in the file"""
        return self.__lines

    def header(self):
        """returns the header in the file"""
        return self.__header


class DelimitedFileMapper:
    """A class that linearly searches a key in a DelimitedFile and for
    the first row found, returns the value in the specified column"""

    def __init__(self, delimited_file, key_column, value_column):
        """Creates an instance of the mapper class using a DelimitedFile"""
        self.__values = {}
        for line in delimited_file.lines():
            self.__values[line[key_column]] = line[value_column]

    def __getitem__(self, key):
        """looks for the key in the key column"""
        return self.__values.get(key, None)

    def items(self):
        """returns the key, value pairs"""
        return self.__values.items()

    def keys(self):
        """returns the keys"""
        return self.__values.keys()

    def __str__(self):
        return str(self.__values)


def levenshtein_distance(str1, str2):
    """computes the Levenshtein distance. This is used in order
    to make approximate string comparisons"""
    strlen1 = len(str1)
    strlen2 = len(str2)

    dist = [[0 for _ in range(strlen2 + 1)] for _ in range(strlen1 + 1)]
    for row in range(strlen1 + 1):
        dist[row][0] = row
    for col in range(strlen2 + 1):
        dist[0][col] = col

    for col in range(strlen2):
        for row in range(strlen1):
            if str1[row] == str2[col]:
                dist[row + 1][col + 1] = dist[row][col]
            else:
                dist[row + 1][col + 1] = min(dist[row][col + 1] + 1,
                                             dist[row + 1][col] + 1,
                                             dist[row][col] + 1)
    return dist[strlen1][strlen2]


class RankedAnchor:  # pylint: disable-msg=R0903
    """A hyperlink with a Levenshtein score"""
    def __init__(self, score, anchor):
        """Creates a GenomeListEntry with a given Levenshtein distance
        and a BeautifulSoup anchor tag"""
        self.score = score
        self.anchor = anchor

    def __str__(self):
        return "(score: %d, %s)" % (self.score, str(self.anchor))

    def __repr__(self):
        return str(self)


def best_matching_links(search_string, html):
    """given a search string and an HTML text, extract the best matching
    href"""
    soup = bs.BeautifulSoup(html)
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
    return scipy.stats.scoreatpercentile(values, probability * 100)


def r_stddev(values):
    """This is a standard deviation function, adjusted so it will
    return approximately the same value as R's sd() function would"""
    num_values = len(values)
    return round(numpy.std(values) / numpy.sqrt(float(num_values - 1) /
                                                float(num_values)), 8)


def r_variance_columns(matrix):
    """computes the variance over the columns of a matrix, applying
    a bias of (n/n-1) over the results to match with R"""
    num_rows = len(matrix)
    bias = float(num_rows) / float(num_rows - 1)
    result = numpy.var(matrix, 0)
    return [(value * bias) for value in result]


def column_means(matrix):
    """computes the column means of a matrix"""
    return numpy.mean(matrix, axis=0)


def row_means(matrix):
    """computes the row means of a matrix"""
    return numpy.mean(matrix, axis=1)


def make_matrix(row_names, num_columns, init_value=0):
    """creates a two-dimensional matrix with len(row_names) rows and
    num_cols columns, where all fields are initialized with
    init_value. The rows are accessed by row name"""
    result = {}
    for name in row_names:
        result[name] = [init_value for _ in range(num_columns)]
    return result


class DocumentNotFound(Exception):
    """An exception indicating that the requested document does not exist"""
    pass


class CMonkeyURLopener(urllib.FancyURLopener):
    """An URL opener that can detect 404 errors"""

    def http_error_default(self, url, fp, errcode, errmsg, headers):
        # pylint: disable-msg=R0913
        # pylint: disable-msg=C0103
        """overriding the default error handling method to handle HTTP 404
        errors"""
        if (errcode == 404):
            raise DocumentNotFound(url)

        # call super class handler.
        # note that urllib.FancyURLopener is not a new-style class
        return urllib.FancyURLopener.http_error_default(
            self, url, fp, errcode, errmsg, headers)


def read_url(url):
    """convenience method to read a document from a URL using the
    CMonkeyURLopener"""
    return CMonkeyURLopener().open(url).read()


def read_url_cached(url, cache_filename):
    """convenience method to read a document from a URL using the
    CMonkeyURLopener, cached version"""
    if not os.path.exists(cache_filename):
        CMonkeyURLopener().retrieve(url, cache_filename)
    with open(cache_filename) as cached_file:
        return cached_file.read()


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
        """Returns the keys of the wrapped dictionary"""
        return self.__wrapped_dict.keys()


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
    values = sorted(values, reverse=True)
    if not values or len(values) == 0:
        return 0
    if trim == 0.5:
        return numpy.median(values)

    k_floor = int(math.floor(trim * len(values)))
    trim_values_floor = values[k_floor:len(values) - k_floor]
    return numpy.mean(trim_values_floor)


######################################################################
### RPY2 abstraction
######################################################################
def density(kvalues, cluster_values, bandwidth, dmin, dmax):
    """generic function to compute density scores"""
    r_density = robjects.r['density']
    kwargs = {'bw': bandwidth, 'adjust': 2, 'from': dmin,
              'to': dmax, 'n': 256, 'na.rm': True}
    # x = dvalues[0], y = dvalues[1]
    dvalues = r_density(robjects.FloatVector(cluster_values), **kwargs)
    r_rev = robjects.r['rev']
    r_cumsum = robjects.r['cumsum']
    r_approx = robjects.r['approx']
    r_pvalues = r_approx(dvalues[0], r_rev(r_cumsum(r_rev(dvalues[1]))),
                         kvalues)[1]
    pvalues = [pvalue for pvalue in r_pvalues]
    psum = sum(pvalues)
    return [(pvalue / psum) for pvalue in pvalues]


def rnorm(num_values, std_deviation):
    """returns the result of R's rnorm function"""
    r_rnorm = robjects.r['rnorm']
    kwargs = {'sd': std_deviation}
    return r_rnorm(num_values, **kwargs)


def phyper(q, m, n, k, lower_tail=False):
    r_phyper = robjects.r['phyper']
    kwargs = {'lower.tail': lower_tail}
    return r_phyper(robjects.FloatVector(q),
                    robjects.FloatVector(m),
                    robjects.FloatVector(n),
                    robjects.FloatVector(k), **kwargs)

######################################################################
### Misc functionality
######################################################################


def add_if_unique(sequence, item):
    """add the item to the Python sequence only if it does not exist"""
    if item not in sequence:
        sequence.append(item)


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


__all__ = ['DelimitedFile', 'best_matching_links', 'quantile', 'make_matrix',
           'DocumentNotFound', 'CMonkeyURLopener', 'read_url',
           'read_url_cached', 'ThesaurusBasedMap', 'trim_mean']
