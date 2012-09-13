# vi: sw=4 ts=4 et:
"""cmonkey.py - cMonkey top-level module

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import os
import os.path
import sys
import logging
import json
import cmonkey_run
import datamatrix as dm
import util

CMONKEY_VERSION = '4.0'
CHECKPOINT_INTERVAL = 100
CHECKPOINT_FILE = None


if __name__ == '__main__':
    print('cMonkey (Python port) (c) 2011-2012, Institute for Systems Biology')
    print('This program is licensed under the General Public License V3.')
    print('See README and LICENSE for details.\n')
    if len(sys.argv) < 3:
        print('Usage: ./run_cmonkey.sh <organism-code> <ratio-file> ' +
              '[string-file] [checkpoint-file]')
    else:
        string_file=None
        if len(sys.argv) > 3:
            string_file = sys.argv[3]
        if len(sys.argv) > 4:
            CHECKPOINT_FILE = sys.argv[4]

        matrix_factory = dm.DataMatrixFactory([dm.nochange_filter, dm.center_scale_filter])
        infile = util.DelimitedFile.read(sys.argv[2], has_header=True, quote='\"')
        matrix = matrix_factory.create_from(infile)
        cmonkey_run = cmonkey_run.CMonkeyRun(sys.argv[1], matrix,
                                             string_file=string_file)  # num_cluster=250 for halo_ref
        cmonkey_run.run()
