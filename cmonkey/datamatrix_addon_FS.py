'''
Created on Apr 3, 2012

@author: frank
'''

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
import datamatrix as dm_old


class DataMatrix_FS(dm_old.DataMatrix):
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
        dm_old.DataMatrix.__init__(self, nrows, ncols, row_names = row_names,
                                   col_names = col_names, values = values, init_value = init_value)
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
            self.__row_names = np.array(["Row " + str(i)
                                         for i in range(nrows)])
        else:
            if len(row_names) != nrows:
                raise ValueError("number of row names should be %d" % nrows)
            self.__row_names = np.array(row_names)

        
        self.__rownames_original = self.__row_names         # Frank Schmitz
                                                            # becomes important when the rownames are tranlated for
                                                            # internal use


        if col_names == None:
            self.__column_names = np.array(["Col " + str(i)
                                            for i in xrange(ncols)])
        else:
            if len(col_names) != ncols:
                raise ValueError("number of column names should be %d" % ncols)
            self.__column_names = np.array(col_names)

        if values != None:
            check_values()
            self.__values = np.array(values, dtype=np.float64)
        else:
            self.__values = np.zeros((nrows, ncols))
            if init_value != None:
                self.__values.fill(init_value)
        self.__row_variance = None

    
    ### Frank Schmitz
    # want to add addtnl functionality - translate the genes to
    # the common nomenclature, improves performance
     
    def set_row_names(self, new_row_names):
        """
        set the rownames new, leave self.__rownames_original untouched
        
        """
        self.__row_names = new_row_names

    def translateRowNames(self, thesaurus):
        """ Frank Schmitz:
            translate the rownames into common used IDs, for later speed improvement and compatibility
            leave the self.___rownames_originals untouched
        """
        
        logging.info("\x1b[31mmatrix:\t\x1b[0mtranslating matrix")
        nrowL = []
        for rowname in self.__row_names:
            try:
                nrowL.append(thesaurus[rowname])
            except KeyError:
                if rowname in thesaurus.keys():
                    nrowL.append(rowname)
                else:
                    nrowL.append("".join(["unMapped_", rowname]))
        self.set_row_names(np.array(nrowL))


class DataMatrixFactory_FS(dm_old.DataMatrixFactory):
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
        dm_old.DataMatrixFactory.__init__(self, filters)
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
        data_matrix = DataMatrix_FS(nrows, ncols, rownames, colnames,
                                 values=values)

        for matrix_filter in self.filters:
            data_matrix = matrix_filter(data_matrix)
        return data_matrix
