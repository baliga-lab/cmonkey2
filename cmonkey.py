#!/usr/bin/env python
# vi: sw=4 ts=4 et:
"""cmonkey.py - cMonkey top-level module

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import os.path
import cmonkey.cmonkey_run as cmr
import cmonkey.datamatrix as dm
import cmonkey.config as conf
import cmonkey.util as util
import argparse
import logging
from cmonkey.schedule import make_schedule
import ConfigParser
import tempfile
import cmonkey.scoring as scoring
import random


# if we were installed through Debian package management, default.ini is found here
SYSTEM_INI_PATH = '/etc/cmonkey-python/default.ini'
USER_INI_PATH = 'config/default.ini'

if __name__ == '__main__':
    description = """cMonkey (Python port) (c) 2011-2012,
Institute for Systems Biology
This program is licensed under the General Public License V3.
See README and LICENSE for details.\n"""

    # read default configuration parameters
    config = ConfigParser.ConfigParser()
    if os.path.exists(USER_INI_PATH):
        config.read(USER_INI_PATH)
    elif os.path.exists(SYSTEM_INI_PATH):
        config.read(SYSTEM_INI_PATH)
    else:
        raise Exception('could not find default.ini !')

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--ratios', required=True,
                        help='tab-separated ratios matrix file')

    parser.add_argument('--organism', help='KEGG organism code', default=None)
    parser.add_argument('--out', default=config.get("General", "output_dir"),
                        help='output directory')
    parser.add_argument('--cachedir', default=config.get("General", "cache_dir"),
                        help="path to cache directory")
    parser.add_argument('--string', help='tab-separated STRING file for the organism',
                        default=None)
    parser.add_argument('--operons', help='tab-separated STRING file for the organism',
                        default=None)
    parser.add_argument('--checkpoint', help='checkpoint-file')
    parser.add_argument('--checkratios', action="store_true",
                        help='check gene expression quality')
    parser.add_argument('--remap_network_nodes', action="store_true",
                        help='network nodes are not named to RSAT primary names')
    parser.add_argument('--logfile', default=None, help="""path to log file""")
    parser.add_argument('--keep_memeout', action="store_true",
                        help="""keep MEME output files""")
    parser.add_argument('--ncbi_code', default=None, help="NCBI taxonomy id")
    parser.add_argument('--numclusters', type=int,
                        default=None, help="override the number of clusters")

    parser.add_argument('--nomotifs', action="store_true", help="deactivate motif scoring")
    parser.add_argument('--nonetworks', action="store_true", help="deactivate network scoring")
    parser.add_argument('--nostring', action="store_true", help="deactivate STRING network scoring")
    parser.add_argument('--nooperons', action="store_true", help="deactivate operon network scoring")
    parser.add_argument('--config', default=None, help="additional configuration file")
    parser.add_argument('--debug', action="store_true",
                        help="""run in debug mode""")
    parser.add_argument('--random_seed', type=int)

    # RSAT overrides
    parser.add_argument('--rsat_dir', default=None,
                        help="""RSAT override: data directory""")
    parser.add_argument('--rsat_organism', default=None,
                        help="""override the RSAT organism name""")

    args = parser.parse_args()

    # no organism provided -> dummy organism
    if args.organism is None:
        print("WARNING - no organism provided - assuming that you want to score ratios only or don't use automatic download")
        if not args.rsat_dir:
            args.nomotifs = True
        if not args.string and not args.operons:
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

    # override number of clusters either on the command line or through
    # the config file
    try:
        num_clusters = config.getint("General", "num_clusters")
    except:
        num_clusters = args.numclusters

    cmonkey_run = cmr.CMonkeyRun(args.organism, matrix,
                                 string_file=args.string,
                                 rsat_organism=args.rsat_organism,
                                 log_filename=args.logfile,
                                 remap_network_nodes=args.remap_network_nodes,
                                 ncbi_code=args.ncbi_code,
                                 num_clusters=num_clusters,
                                 operon_file=args.operons,
                                 rsat_dir=args.rsat_dir)
    conf.set_config(cmonkey_run, config)

    cmonkey_run['output_dir'] = args.out
    cmonkey_run['cache_dir'] = args.cachedir
    cmonkey_run['debug'] = args.debug
    cmonkey_run['keep_memeout'] = args.keep_memeout or args.debug
    cmonkey_run['nonetworks'] = args.nonetworks
    cmonkey_run['nomotifs'] = args.nomotifs or not cmonkey_run['meme_version']
    cmonkey_run['use_string'] = not args.nostring
    cmonkey_run['use_operons'] = not args.nooperons
    if args.random_seed:
        cmonkey_run['random_seed'] = args.random_seed

    if cmonkey_run['random_seed']:
        random.seed(cmonkey_run['random_seed'])
        util.r_set_seed(cmonkey_run['random_seed'])

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

    # Set update frequency to every iteration, so the full results are written
    if cmonkey_run['debug']:
        cmonkey_run['stats_freq'] = 1
        cmonkey_run['result_freq'] = 1


    if not proceed:
        logging.error("# genes found: %d, # total: %d, %f %% - please check your ratios file",
                      num_found, total, percent)
    else:
        if args.checkpoint:
            cmonkey_run.run_from_checkpoint(args.checkpoint)
        else:
            cmonkey_run.run()
