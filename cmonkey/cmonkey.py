# vi: sw=4 ts=4 et:
"""cmonkey.py - cMonkey top-level module

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import os.path
import cmonkey_run
import datamatrix as dm
import util
import argparse


if __name__ == '__main__':
    description = """cMonkey (Python port) (c) 2011-2012,
Institute for Systems Biology
This program is licensed under the General Public License V3.
See README and LICENSE for details.\n"""

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--organism', required=True, help='KEGG organism code')
    parser.add_argument('--ratios', required=True,
                        help='tab-separated ratios matrix file')
    parser.add_argument('--out', default='out', help='output directory')
    parser.add_argument('--string',
                        help='tab-separated STRING file for the organism')
    parser.add_argument('--checkpoint', help='checkpoint-file')
    args = parser.parse_args()

    string_file = args.string
    CHECKPOINT_FILE = args.checkpoint

    matrix_factory = dm.DataMatrixFactory([dm.nochange_filter,
                                           dm.center_scale_filter])
    matrix_filename = args.ratios

    if matrix_filename.startswith('http://'):
        indata = util.read_url(matrix_filename)
        infile = util.dfile_from_text(indata, has_header=True, quote='\"')
    else:
        infile = util.read_dfile(matrix_filename, has_header=True, quote='\"')

    matrix = matrix_factory.create_from(infile)
    # num_cluster=250 for halo_ref
    cmonkey_run = cmonkey_run.CMonkeyRun(args.organism, matrix,
                                         string_file=string_file)
    cmonkey_run['output_dir'] = args.out
    cmonkey_run['out_database'] = os.path.join(args.out, 'cmonkey_run.db')
    if args.checkpoint:
        cmonkey_run.run_from_checkpoint(args.checkpoint)
    else:
        cmonkey_run.run()
