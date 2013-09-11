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
import ConfigParser


if __name__ == '__main__':
    description = """cMonkey (Python port) (c) 2011-2012,
Institute for Systems Biology
This program is licensed under the General Public License V3.
See README and LICENSE for details.\n"""

    # read default configuration parameters
    config = ConfigParser.ConfigParser()
    config.read('default.ini')

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--ratios', required=True,
                        help='tab-separated ratios matrix file')

    parser.add_argument('--organism', help='KEGG organism code', default=None)
    parser.add_argument('--out', default=config.get("General", "output_dir"),
                        help='output directory')
    parser.add_argument('--cachedir', default=config.get("General", "cache_dir"),
                        help="path to cache directory")
    parser.add_argument('--string', help='tab-separated STRING file for the organism')
    parser.add_argument('--checkpoint', help='checkpoint-file')
    parser.add_argument('--checkratios', action="store_true",
                        help='check gene expression quality')
    parser.add_argument('--remap_network_nodes', action="store_true",
                        help='network nodes are not named to RSAT primary names')
    parser.add_argument('--rsat_organism', default=None,
                        help="""override the RSAT organism name""")
    parser.add_argument('--logfile', default=None, help="""path to log file""")
    parser.add_argument('--keep_memeout', action="store_true",
                        help="""keep MEME output files""")
    parser.add_argument('--ncbi_code', default=None, help="path to cache directory")

    parser.add_argument('--nomotifs', action="store_true", help="deactivate motif scoring")
    parser.add_argument('--nonetworks', action="store_true", help="deactivate network scoring")
    parser.add_argument('--nostring', action="store_true", help="deactivate STRING network scoring")
    parser.add_argument('--nooperons', action="store_true", help="deactivate operon network scoring")
    parser.add_argument('--config', default=None, help="additional configuration file")
    args = parser.parse_args()

    # no organism provided -> dummy organism
    if args.organism == None:
        print("WARNING - no organism provided - assuming that you want to score ratios only")
        args.nomotifs = True
        args.nonetworks = True

    # user overrides in config files
    if args.config:
        config.read(args.config)

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
                                         string_file=args.string,
                                         rsat_organism=args.rsat_organism,
                                         log_filename=args.logfile,
                                         remap_network_nodes=args.remap_network_nodes,
                                         ncbi_code=args.ncbi_code)
    cmonkey_run['output_dir'] = args.out
    cmonkey_run['cache_dir'] = args.cachedir
    cmonkey_run['num_iterations'] = config.getint("General", "num_iterations")
    cmonkey_run['out_database'] = os.path.join(args.out,
                                               config.get("General", "dbfile_name"))
    cmonkey_run['multiprocessing'] = config.getboolean('General', 'use_multiprocessing')
    cmonkey_run['quantile_normalize'] = config.getboolean('Scoring', 'quantile_normalize')

    cmonkey_run['keep_memeout'] = args.keep_memeout
    cmonkey_run['donetworks'] = not args.nonetworks
    cmonkey_run['domotifs'] = not args.nomotifs
    cmonkey_run['use_string'] = not args.nostring
    cmonkey_run['use_operons'] = not args.nooperons

    proceed = True
    checkratios = args.checkratios
    
    if args.checkratios:
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
