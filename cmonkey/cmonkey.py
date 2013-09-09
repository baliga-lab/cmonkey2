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
import logging


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
    parser.add_argument('--no_operons', default="False",
                        help="""don't use operons in retrieving features and networks""")
    parser.add_argument('--checkratios', default="True",
                        help='check gene expression quality')
    parser.add_argument('--remap_network_nodes', default="False",
                        help='network nodes are not named to RSAT primary names')
    parser.add_argument('--rsat_organism', default=None,
                        help="""override the RSAT organism name""")
    parser.add_argument('--logfile', default=None,
                        help="""path to log file""")
    parser.add_argument('--keep_memeout', default="False",
                        help="""keep MEME output files""")
    parser.add_argument('--cachedir', default="cache",
                        help="path to cache directory")
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
    infile = None
    # num_cluster=250 for halo_ref
    cmonkey_run = cmonkey_run.CMonkeyRun(args.organism, matrix,
                                         string_file=string_file,
                                         use_operons=args.no_operons == "False",
                                         rsat_organism=args.rsat_organism,
                                         log_filename=args.logfile,
                                         remap_network_nodes=args.remap_network_nodes == "True")
    cmonkey_run['output_dir'] = args.out
    cmonkey_run['cache_dir'] = args.cachedir
    cmonkey_run['out_database'] = os.path.join(args.out, 'cmonkey_run.db')
    cmonkey_run['keep_memeout'] = args.keep_memeout == "True"

    proceed = True
    checkratios = args.checkratios == "True"
    
    if checkratios:
        thesaurus = cmonkey_run.organism().thesaurus()
        logging.info("Checking the quality of the input matrix names...")
        found = [name for name in matrix.row_names if name in thesaurus]
        num_found = len(found)
        total = len(matrix.row_names)
        percent = (float(num_found) / float(total)) * 100.0
        proceed = percent > 50.0

    if not proceed:
       logging.error("""# genes found: %d, # total: %d, %f %% - please check
 your ratios file""",
                     num_found, total, percent)
    else:
        if args.checkpoint:
            cmonkey_run.run_from_checkpoint(args.checkpoint)
        else:
            cmonkey_run.run()
