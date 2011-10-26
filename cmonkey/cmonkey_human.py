"""cmonkey_human.py - cMonkey human module

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import sys
import os
import logging
import datamatrix as dm
import util
import human

LOG_FORMAT = '%(asctime)s %(levelname)-8s %(message)s'
CMONKEY_VERSION = '4.0'
CACHE_DIR = 'humancache'
CONTROLS_FILE = 'human_data/controls.csv'
RUG_FILE = 'human_data/rug.csv'
RUG_PROPS = ['MIXED', 'ASTROCYTOMA', 'GBM', 'OLIGODENDROGLIOMA']

def run_cmonkey():
    """init of the cMonkey system"""
    logging.basicConfig(format=LOG_FORMAT,
                        datefmt='%Y-%m-%d %H:%M:%S',
                        level=logging.DEBUG)
    if not os.path.exists(CACHE_DIR):
        os.mkdir(CACHE_DIR)
    matrix = read_matrix(sys.argv[1])
    logging.info("Filtered input matrix has %d rows and %d columns:",
                 matrix.num_rows(), matrix.num_columns())

def read_controls():
    """reads the controls file"""
    with open(CONTROLS_FILE) as infile:
        return [line.strip() for line in infile.readlines()]


def read_rug(pred):
    """reads the rug file"""
    infile = util.DelimitedFile.read(RUG_FILE, sep=',', has_header=False)
    return list(set([row[0] for row in infile.lines() if pred(row)]))


def read_matrix(filename):
    """reads the data matrix from a file"""
    rug = read_rug(lambda row: row[1] in RUG_PROPS)
    columns_to_use = list(set(rug + read_controls()))

    # pass the column filter as the first filter to the DataMatrixFactory,
    # so normalization will be applied to the submatrix
    matrix_factory = dm.DataMatrixFactory([
            lambda matrix: matrix.submatrix_by_name(column_names=columns_to_use)])
    infile = util.DelimitedFile.read(filename, sep=',', has_header=True,
                                     quote="\"")
    matrix = matrix_factory.create_from(infile)
    return human.select_probes(2000,
                               [1 for _ in range(matrix.num_columns())],
                               matrix)

if __name__ == '__main__':
    print('cMonkey (Python port) (c) 2011, Institute for Systems Biology')
    print('This program is licensed under the General Public License V3.')
    print('See README and LICENSE for details.\n')
    print('(Human version)')
    if len(sys.argv) <= 1:
        print('Usage: ./run_cmonkey_human.sh <ratio-file>')
    else:
        run_cmonkey()
