#!/usr/bin/python

import cmonkey.datamatrix as dm
import cmonkey.util as util
import sys

def read_matrix(filename):
    dfile = util.DelimitedFile.read(filename, has_header=True, sep='\t', quote='"')
    return dm.DataMatrixFactory([]).create_from(dfile).sorted_by_row_name()

def check_dimensions(matrix1, matrix2):
    if (matrix1.num_rows() != matrix2.num_rows() or
        matrix1.num_columns() != matrix2.num_columns()):
        print ("Dimensions do not match (matrix1: %dx%d, matrix2: %dx%d)" %
               (matrix1.num_rows(), matrix1.num_columns(),
                matrix2.num_rows(), matrix2.num_columns()))
        return False
    return True

def check_names(matrix1, matrix2):
    result = True
    for row in range(matrix1.num_rows()):
        if matrix1.row_name(row) != matrix2.row_name(row):
            print "Row mismatch in %d: '%s' != '%s'" % (row, matrix1.row_name(row),
                                                        matrix2.row_name(row))
            result = False

    return result

#EPS = 1.0e-5
EPS = 0.14

def check_values(matrix1, matrix2):
    result = True
    for row in range(matrix1.num_rows()):
        for col in range(matrix1.num_columns()):
            diff = abs(matrix1[row][col] - matrix2[row][col])
            #print diff
            if diff > EPS:
                print "Value mismatch at (%s, cluster %d): %f != %f (diff = %f)" % (
                    matrix1.row_name(row), col + 1, matrix1[row][col], matrix2[row][col], diff)
                result = False
    return result

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print "Usage: comp_matrices.py <matrix1> <matrix2>"
    else:
        matrix1 = read_matrix(sys.argv[1])
        matrix2 = read_matrix(sys.argv[2])
        result = check_dimensions(matrix1, matrix2)
        if result:
            print "Checking names of (%dx%d) matrices..." % (matrix1.num_rows(), matrix2.num_columns())
            result = check_names(matrix1, matrix2)
        if result:
            print "Checking values..."
            result = check_values(matrix1, matrix2)
    print "Done."
