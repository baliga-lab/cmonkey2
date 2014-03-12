"""config.py - cMonkey configuration module

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import os
import argparse
import ConfigParser
import logging
import tempfile

from cmonkey.schedule import make_schedule
import cmonkey.util as util
import cmonkey.datamatrix as dm
import meme

LOG_FORMAT = '%(asctime)s %(levelname)-8s %(message)s'

DESCRIPTION = """cMonkey (Python port) (c) 2011-2012,
Institute for Systems Biology
This program is licensed under the General Public License V3.
See README and LICENSE for details.\n"""

# if we were installed through Debian package management, default.ini is found here
SYSTEM_INI_PATH = '/etc/cmonkey-python/default.ini'
USER_INI_PATH = 'config/default.ini'


def __set_config(config):
    """Returns a dictionary containing the configuration contained in
    the config parser object. Note that there are only 3 fixed sections:
    General, Membership and Scoring"""
    params = {}

    def set_scaling(section):
        try:
            params[section]['scaling'] = ('scaling_const', config.getfloat(section, 'scaling_const'))
            return
        except:
            pass
        try:
            params[section]['scaling'] = ('scaling_rvec', config.get(section, 'scaling_rvec'))
        except:
            raise Exception("no scaling found for section '%s'" % section)

    # override directories
    tmp_dir = config.get('General', 'tmp_dir')
    if tmp_dir:
        tempfile.tempdir = tmp_dir
    params['output_dir'] = config.get('General', 'output_dir')
    params['cache_dir'] = config.get('General', 'cache_dir')

    params['num_iterations'] = config.getint("General", "num_iterations")
    params['start_iteration'] = config.getint("General", "start_iteration")
    params['multiprocessing'] = config.getboolean('General', 'use_multiprocessing')
    params['postadjust'] = config.getboolean('General', 'postadjust')
    params['log_subresults'] = config.getboolean('General', 'log_subresults')
    params['add_fuzz'] = config.get('General', 'add_fuzz')
    params['checkpoint_interval'] = config.getint('General', 'checkpoint_interval')
    try:
        params['random_seed'] = config.getint('General', 'random_seed')
    except:
        params['random_seed'] = None

    # Quantile normalization is false by default in cMonkey-R
    params['quantile_normalize'] = config.getboolean('Scoring', 'quantile_normalize')

    # membership default parameters
    params['memb.min_cluster_rows_allowed'] = config.getint('Membership', 'min_cluster_rows_allowed')
    params['memb.max_cluster_rows_allowed'] = config.getint('Membership', 'max_cluster_rows_allowed')
    params['memb.prob_row_change'] = config.getfloat('Membership', 'probability_row_change')
    params['memb.prob_col_change'] = config.getfloat('Membership', 'probability_column_change')
    params['memb.max_changes_per_row'] = config.getint('Membership', 'max_changes_per_row')
    params['memb.max_changes_per_col'] = config.getint('Membership', 'max_changes_per_column')

    ids = [section for section in config.sections()
           if section not in {'General', 'Scoring', 'Membership'}]

    # processing scoring function specific stuff
    for id in ids:
        params[id] = {}
        for option, value in config.items(id):
            if option == 'schedule':
                params[id]['schedule'] = make_schedule(value)
            elif option.startswith('scaling_'):
                set_scaling(id)
            else:
                params[id][option] = value

    params['sequence_types'] = config.get('Motifs', 'sequence_types').split(',')
    params['search_distances'] = {}
    params['scan_distances'] = {}
    for seqtype in params['sequence_types']:
        cat = "SequenceType-%s" % seqtype
        params['search_distances'][seqtype] = tuple(
            map(int, config.get(cat, 'search_distance').split(',')))
        params['scan_distances'][seqtype] = tuple(
            map(int, config.get(cat, 'scan_distance').split(',')))

    params['stats_freq'] = config.getint('General', 'stats_frequency')
    params['result_freq'] = config.getint('General', 'result_frequency')
    
    return params

def __get_config_parser():
    # read default configuration parameters
    config = ConfigParser.ConfigParser()
    if os.path.exists(USER_INI_PATH):
        config.read(USER_INI_PATH)
    elif os.path.exists(SYSTEM_INI_PATH):
        config.read(SYSTEM_INI_PATH)
    else:
        raise Exception('could not find default.ini !')
    return config


def __get_arg_parser():
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    parser.add_argument('--ratios', required=True,
                        help='tab-separated ratios matrix file')

    parser.add_argument('--organism', help='KEGG organism code', default=None)
    parser.add_argument('--out', help='output directory')
    parser.add_argument('--cachedir', help="path to cache directory")
    parser.add_argument('--string', help='tab-separated STRING file for the organism',
                        default=None)
    parser.add_argument('--operons', help='tab-separated operons file for the organism',
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
    parser.add_argument('--pipeline', default=None,
                        help="""override the scoring pipeline""")
    return parser


def __num_clusters(config, args, ratios):
    """override number of clusters either on the command line or through
    the config file"""
    try:
        num_clusters = config.getint("General", "num_clusters")
    except:
        num_clusters = args.numclusters
    if num_clusters is None:
        return int(round(ratios.num_rows * args.clusters_per_row / 20.0))
    else:
        return num_clusters
    

def setup():
    """main configuration function - does everything. It reads and configures
    everything that can be derived from the configuration files and input
    parameters and returns the configuration in a dictionary as well as
    the normalized ratios matrix"""

    config_parser = __get_config_parser()
    arg_parser = __get_arg_parser()
    args = arg_parser.parse_args()
    logging.basicConfig(format=LOG_FORMAT, datefmt='%Y-%m-%d %H:%M:%S',
                        level=logging.DEBUG, filename=args.logfile)

    # no organism provided -> dummy organism
    if args.organism is None:
        print("WARNING - no organism provided - assuming that you want to score ratios only or don't use automatic download")
        if not args.rsat_dir:
            args.nomotifs = True
        if not args.string and not args.operons:
            args.nonetworks = True

    # user overrides in config files
    if args.config:
        config_parser.read(args.config)

    # Initial configuration from default + user config
    params = __set_config(config_parser)
    matrix_factory = dm.DataMatrixFactory([dm.nochange_filter,
                                           dm.center_scale_filter])
    matrix_filename = args.ratios

    if matrix_filename.startswith('http://'):
        indata = util.read_url(matrix_filename)
        infile = util.dfile_from_text(indata, has_header=True, quote='\"')
    else:
        infile = util.read_dfile(matrix_filename, has_header=True, quote='\"')

    ratios = matrix_factory.create_from(infile)
    infile = None

    args.clusters_per_row = 2
    
    """The overrides dictionary holds all the values that will overwrite or add
    to the settings defined in the default and user-defined ini files
    """
    overrides = {'organism_code': args.organism, 'string_file': args.string,
                 'logfile': args.logfile, 'rsat_organism': args.rsat_organism,
                 'num_clusters': __num_clusters(config_parser, args, ratios),
                 'memb.clusters_per_row': args.clusters_per_row,
                 'remap_network_nodes': args.remap_network_nodes,
                 'ncbi_code': args.ncbi_code,
                 'operon_file': args.operons,
                 'rsat_dir': args.rsat_dir,
                 'use_operons': True, 'use_string': True, 'global_background': True,
                 'meme_version': meme.check_meme_version(),
                 'debug': args.debug,
                 'keep_memeout': args.debug or args.keep_memeout,
                 'nomotifs': False,
                 'nonetworks': args.nonetworks,
                 'checkratios': args.checkratios,
                 'checkpoint': args.checkpoint, 'random_seed': args.random_seed,
                 'pipeline_file': args.pipeline}

    # membership update default parameters
    # these come first, since a lot depends on clustering numbers
    num_clusters = overrides['num_clusters']
    if ratios.num_columns >= 60:
        overrides['memb.clusters_per_col'] = int(round(num_clusters / 2.0))
    else:
        overrides['memb.clusters_per_col'] = int(round(num_clusters * 2.0 / 3.0))

    overrides['nomotifs'] = args.nomotifs or not overrides['meme_version']
    overrides['use_string'] = not args.nostring
    overrides['use_operons'] = not args.nooperons
    if args.out:
        overrides['output_dir'] = args.out
    if args.cachedir:
        overrides['cache_dir'] = args.cachedir

    if overrides['random_seed']:
        random.seed(overrides['random_seed'])
        util.r_set_seed(overrides['random_seed'])

    # Set update frequency to every iteration, so the full results are written
    if overrides['debug']:
        overrides['stats_freq'] = 1
        overrides['result_freq'] = 1
    for key, value in overrides.items():
        params[key] = value

    params['out_database'] = os.path.join(params['output_dir'],
                                          config_parser.get("General", "dbfile_name"))

    return params, ratios
