'''
Created on Apr 10, 2012

@author: frank schmitz, sbri, 2012

'''

""" different functions to load the data into cmonkey"""

import util
import datamatrix_addon_FS as dm






def read_matrix(filename, RUG_FILE, RUG_PROPS):
    """reads the data matrix from a file"""
    '''
    '''
    def read_rug(pred):
        """reads the rug file"""
        infile = util.DelimitedFile.read(RUG_FILE, sep='\t', has_header=False)      # Frank Schmitz changed separator to TAB delimited
        return list(set([row[0] for row in infile.lines() if pred(row)]))
    rug = read_rug(lambda row: row[1] in RUG_PROPS)
    columns_to_use = list(set(rug))
    # pass the column filter as the first filter to the DataMatrixFactory,
    # so normalization will be applied to the submatrix
    matrix_factory = dm.DataMatrixFactory_FS([
            lambda matrix: matrix.submatrix_by_name(
                column_names=columns_to_use)])
    infile = util.DelimitedFile.read(filename, sep='\t', has_header=True,   #FS changed to TAB separator
                                     quote="\"")
    matrix = matrix_factory.create_from(infile)

    column_groups = {1: range(matrix.num_columns())}
    return matrix
