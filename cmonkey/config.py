"""config.py - cMonkey configuration module

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import os
import sys
import argparse
import ConfigParser
import logging
import tempfile
import json
import random

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

"""
supported debug options:

  - all: all options turned on
  - keep_memeout: keeps meme output files
  - dump_results: dump results into cmresults files, implies keep_memeout
  - dump_scores: dumps score matrices as received from the individual scoring functions
  - profile_mem: write memory profile
  - random_seed: fixed random seed
"""
ALL_DEBUG_OPTIONS = {'keep_memeout', 'dump_results', 'dump_scores', 'profile_mem',
                     'random_seed', 'keep_mastout'}


def get_config_boolean(config, section, option, default_value):
    """tries to retrieve a boolean value from config or returns the default"""
    try:
        return config.getboolean(section, option)
    except:
        return default_value


def get_config_int(config, section, option, default_value=None):
    """tries to retrieve an int value from config or returns the default"""
    try:
        return config.getint(section, option)
    except:
        return default_value

def get_config_str(config, section, option, default_value=None):
    """tries to retrieve a str value from config or returns the default"""
    try:
        return config.get(section, option)
    except:
        return default_value


def set_config(config):
    """Returns a dictionary containing the configuration contained in
    the config parser object. Note that there are only 3 fixed sections:
    General, Membership and Scoring"""
    params = {}
    set_config_general(config, params)
    params['quantile_normalize'] = config.getboolean('Scoring', 'quantile_normalize')
    set_config_membership(config, params)
    set_config_scoring_functions(config, params)
    set_config_motifs(config, params)
    return params


def set_config_general(config, params):
    """Process General section"""
    # override directories
    tmp_dir = config.get('General', 'tmp_dir')
    if tmp_dir:
        tempfile.tempdir = tmp_dir
        
    params['output_dir'] = config.get('General', 'output_dir')
    params['cache_dir'] = config.get('General', 'cache_dir')
    params['tmp_dir'] = tmp_dir
    params['dbfile_name'] = config.get('General', 'dbfile_name')
    params['normalize_ratios'] = config.getboolean('General', 'normalize_ratios')
    params['num_iterations'] = config.getint("General", "num_iterations")
    params['start_iteration'] = config.getint("General", "start_iteration")
    params['multiprocessing'] = config.getboolean('General', 'use_multiprocessing')
    params['case_sensitive'] = config.getboolean('General', 'case_sensitive')
    params['num_cores'] = get_config_int('General', 'num_cores', None)
    params['postadjust'] = config.getboolean('General', 'postadjust')
    params['log_subresults'] = config.getboolean('General', 'log_subresults')
    params['add_fuzz'] = config.get('General', 'add_fuzz')

    # python can have large seeds, R, however has a 32 bit limit it seems
    params['random_seed'] = get_config_int(config, 'General', 'random_seed',
                                           util.current_millis() % 2147483647)
    params['stats_freq'] = config.getint('General', 'stats_frequency')
    params['result_freq'] = config.getint('General', 'result_frequency')
    params['debug_freq'] = config.getint('General', 'debug_frequency')

    # implicit parameters for compatibility
    params['use_operons'] = get_config_boolean(config, 'General', 'use_operons', True)
    params['use_string'] = get_config_boolean(config, 'General', 'use_string', True)
    params['checkratios'] = get_config_boolean(config, 'General', 'checkratios', False)
    params['organism_code'] = get_config_str(config, 'General', 'organism_code', None)


def set_config_membership(config, params):
    """membership default parameters"""
    params['memb.min_cluster_rows_allowed'] = config.getint('Membership',
                                                            'min_cluster_rows_allowed')
    params['memb.max_cluster_rows_allowed'] = config.getint('Membership',
                                                            'max_cluster_rows_allowed')
    params['memb.prob_row_change'] = config.getfloat('Membership', 'probability_row_change')
    params['memb.prob_col_change'] = config.getfloat('Membership', 'probability_column_change')
    params['memb.max_changes_per_row'] = config.getint('Membership', 'max_changes_per_row')
    params['memb.max_changes_per_col'] = config.getint('Membership', 'max_changes_per_column')
    params['memb.clusters_per_row'] = get_config_int(config, 'Membership',
                                                     'clusters_per_row')
    params['memb.clusters_per_col'] = get_config_int(config, 'Membership',
                                                     'clusters_per_column')


def set_config_scoring_functions(config, params):
    """processing scoring function specific stuff"""
    def set_scaling(section):
        try:
            params[section]['scaling'] = ('scaling_const',
                                          config.getfloat(section, 'scaling_const'))
            return
        except:
            pass
        try:
            params[section]['scaling'] = ('scaling_rvec', config.get(section, 'scaling_rvec'))
        except:
            raise Exception("no scaling found for section '%s'" % section)

    ids = [section for section in config.sections()
           if section not in {'General', 'Scoring', 'Membership'}]

    for id in ids:
        params[id] = {}
        for option, value in config.items(id):
            if option == 'schedule':
                params[id]['schedule'] = make_schedule(value)
            elif option.startswith('scaling_'):
                set_scaling(id)
            else:
                params[id][option] = value


def set_config_motifs(config, params):
    """process the Motifs section"""
    params['sequence_types'] = config.get('Motifs', 'sequence_types').split(',')
    params['search_distances'] = {}
    params['scan_distances'] = {}
    for seqtype in params['sequence_types']:
        cat = "SequenceType-%s" % seqtype
        params['search_distances'][seqtype] = tuple(
            map(int, config.get(cat, 'search_distance').split(',')))
        params['scan_distances'][seqtype] = tuple(
            map(int, config.get(cat, 'scan_distance').split(',')))


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


def __get_arg_parser(arg_ext):
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    parser.add_argument('--ratios',
                        help='tab-separated ratios matrix file')

    parser.add_argument('--verbose', action="store_true")
    parser.add_argument('--organism', help='KEGG organism code', default=None)
    parser.add_argument('--out', help='output directory')
    parser.add_argument('--cachedir', help="path to cache directory")
    parser.add_argument('--string', help='tab-separated STRING file for the organism',
                        default=None)
    parser.add_argument('--operons', help='tab-separated operons file for the organism',
                        default=None)
    parser.add_argument('--checkratios', action="store_true",
                        help='check gene expression quality')
    parser.add_argument('--remap_network_nodes', action="store_true",
                        help='network nodes are not named to RSAT primary names')
    parser.add_argument('--logfile', default=None, help="""path to log file""")
    parser.add_argument('--ncbi_code', default=None, help="NCBI taxonomy id")
    parser.add_argument('--numclusters', type=int,
                        default=None, help="override the number of clusters")

    parser.add_argument('--nomotifs', action="store_true", help="deactivate motif scoring")
    parser.add_argument('--nonetworks', action="store_true", help="deactivate network scoring")
    parser.add_argument('--nostring', action="store_true", help="deactivate STRING network scoring")
    parser.add_argument('--nooperons', action="store_true", help="deactivate operon network scoring")
    parser.add_argument('--config', default=None, nargs='*', help="additional configuration file(s)")
    parser.add_argument('--debug', default=None,  help="""run in debug mode, can be keep_memeout,
dump_results, dump_scores, profile_mem, random_seed, keep_mastout, all or a combination""")
    parser.add_argument('--random_seed', type=int)
    parser.add_argument('--num_cores', type=int, default=None)
    parser.add_argument('--minimize_io', action="store_true",
                        help='minimal io setting')

    # RSAT overrides
    parser.add_argument('--rsat_dir', default=None,
                        help="""RSAT override: data directory""")
    parser.add_argument('--rsat_organism', default=None,
                        help="""override the RSAT organism name""")
    parser.add_argument('--rsat_features', default='feature',
                        help="""Gene look up table.  Aternative 'cds', 'protein_coding' or 'gene' """)
    parser.add_argument('--rsat_URL', default='http://rsat.bigre.ulb.ac.be/rsat/',
                        help="""RSAT mirror. Alternative 'http://embnet.ccg.unam.mx/rsa-tools'""")

    # Synonym override
    parser.add_argument('--synonym_file', default=None, help="synonyms file")
    parser.add_argument('--fasta_file', default=None, help="take sequences from here instead")

    parser.add_argument('--pipeline', default=None,
                        help="""override the scoring pipeline""")

    parser.add_argument('--case_sensitive', action='store_true',
                        help="""override the case sensitive default""")
    parser.add_argument('--interactive', action='store_true',
                        help="""initialize, but don't run""")
    parser.add_argument('--resume', action='store_true',
                        help="""initialize from out directory""")
    parser.add_argument('--num_iterations', type=int, default=None,
                        help="""specify number of iterations""")


    if arg_ext is not None:
        arg_ext(parser)
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
    

def setup(arg_ext=None):
    """main configuration function - does everything. It reads and configures
    everything that can be derived from the configuration files and input
    parameters and returns the configuration in a dictionary as well as
    the normalized ratios matrix.
    The arg_ext parameter can be used to add additional argements to the
    argparser that can be processed outside of the core application."""
    config_parser = __get_config_parser()
    arg_parser = __get_arg_parser(arg_ext)
    args = arg_parser.parse_args()
    args.command_line = ' '.join(sys.argv)
    if args.verbose:
        loglevel = logging.DEBUG
    else:
        loglevel = logging.INFO

    logging.basicConfig(format=LOG_FORMAT, datefmt='%Y-%m-%d %H:%M:%S',
                        level=loglevel, filename=args.logfile)

    if args.resume:
        return setup_resume(args, config_parser)
    else:
        return setup_default(args, config_parser)

def setup_resume(args, config_parser):
    """setup from out directory"""
    outdir = config_parser.get('General', 'output_dir')
    if args.out is not None:
        outdir = args.out
    logging.info("Reading configuration from '%s'...", outdir)
    config_parser.read(os.path.join(outdir, 'final.ini'))
    logging.info('done.')
    params = set_config(config_parser)
    args.ratios = os.path.join(outdir, 'ratios.tsv.gz')
    ratios = read_ratios(params, args)

    params['resume'] = True
    params['out_database'] = os.path.join(params['output_dir'], params['dbfile_name'])
    params['num_clusters'] = config_parser.getint('General', 'num_clusters')
    num_clusters = params['num_clusters']

    # TODO: these need to be restored or it will crash, need to move stuff around
    # needs rework
    params['rsat_dir'] = None
    params['rsat_organism'] = None
    params['operon_file'] = None
    params['string_file'] = None
    params['ncbi_code'] = None
    params['nonetworks'] = False
    params['nomotifs'] = False
    params['synonym_file'] = None
    params['fasta_file'] = None
    params['pipeline_file'] = None
    params['debug'] = []
    params['interactive'] = True
    # TODO END

    return args, params, ratios

def setup_default(args, config_parser):
    """default configuration method"""
    # no organism provided -> dummy organism
    if args.organism is None:
        logging.warn("no organism provided - assuming that you want to score ratios only or don't use automatic download")
        if not args.rsat_dir:
            args.nomotifs = True
        if not args.string and not args.operons:
            args.nonetworks = True

    # user overrides in config files
    if args.config is not None:
        config_parser.read(args.config)

    # Initial configuration from default + user config
    params = set_config(config_parser)
    ratios = read_ratios(params, args)
    args.clusters_per_row = 2

    # debug options
    debug_options = set(args.debug.split(',')) if args.debug is not None else set()
    if 'dump_results' in debug_options:
        debug_options.add('keep_memeout')
    if debug_options == {'all'}:
        debug_options = ALL_DEBUG_OPTIONS
    
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
                 'rsat_URL': args.rsat_URL,
                 'rsat_features': args.rsat_features,
                 'use_operons': True, 'use_string': True,
                 'debug': debug_options,
                 'nomotifs': False,
                 'minimize_io': args.minimize_io,
                 'nonetworks': args.nonetworks,
                 'checkratios': args.checkratios,
                 'random_seed': args.random_seed,
                 'pipeline_file': args.pipeline,
                 'synonym_file': args.synonym_file,
                 'fasta_file': args.fasta_file,
                 'interactive': args.interactive,
                 'resume': args.resume,
                 'case_sensitive': args.case_sensitive,
                 'command_line': args.command_line}

    if overrides['random_seed'] is None:
        del overrides['random_seed']
    if overrides['case_sensitive'] is None:
        del overrides['case_sensitive']

    # membership update default parameters
    # these come first, since a lot depends on clustering numbers
    num_clusters = overrides['num_clusters']
    if ratios.num_columns >= 60:
        overrides['memb.clusters_per_col'] = int(round(num_clusters / 2.0))
    else:
        overrides['memb.clusters_per_col'] = int(round(num_clusters * 2.0 / 3.0))

    params['MEME']['version'] = meme.check_meme_version()
    overrides['nomotifs'] = args.nomotifs or not params['MEME']['version']
    overrides['use_string'] = not args.nostring
    overrides['use_operons'] = not args.nooperons

    if args.num_cores is not None:
        overrides['num_cores'] = args.num_cores
    if args.out:
        overrides['output_dir'] = args.out
    if args.cachedir:
        overrides['cache_dir'] = args.cachedir

    if args.num_iterations is not None:
        overrides['num_iterations'] = args.num_iterations

    for key, value in overrides.iteritems():
        params[key] = value

    if params['random_seed'] is not None:
        random.seed(params['random_seed'])
        util.r_set_seed(params['random_seed'])

    params['out_database'] = os.path.join(params['output_dir'], params['dbfile_name'])

    # we return the args here, too, in case we have some parameters not
    # processed
    return args, params, ratios


def read_ratios(params, args):
    """reading ratios matrix"""
    if params['normalize_ratios']:
        ratio_filters = [dm.nochange_filter, dm.center_scale_filter]
    else:
        ratio_filters = []

    matrix_factory = dm.DataMatrixFactory(ratio_filters)
    matrix_filename = args.ratios

    if matrix_filename.startswith('http://'):
        indata = util.read_url(matrix_filename)
        infile = util.dfile_from_text(indata, has_header=True, quote='\"')
    else:
        infile = util.read_dfile(matrix_filename, has_header=True, quote='\"')

    if params['case_sensitive'] or args.case_sensitive:
        ratios = matrix_factory.create_from(infile, True)
    else:
        ratios = matrix_factory.create_from(infile, False)
    return ratios


def write_setup(config_params):
    """Write all active configuration settings to <outputdir>/final.ini"""
    with open(os.path.join(config_params['output_dir'], 'final.ini'), 'w') as outfile:
        write_general_settings(outfile, config_params)
        write_membership_settings(outfile, config_params)
        outfile.write('\n[Scoring]\n')
        outfile.write('quantile_normalize = %s\n' % str(config_params['quantile_normalize']))
        for key, value in config_params.iteritems():
            if key != 'pipeline' and type(value) is dict:
                write_section(outfile, key, value)

    with open(os.path.join(config_params['output_dir'], 'pipeline.json'), 'w') as outfile:
        outfile.write(json.dumps(config_params['pipeline']))


def write_general_settings(outfile, config_params):
    """Writes the General section of the final.ini file"""
    def strparam(value):
        return '' if not value else value

    outfile.write('[General]\n')
    outfile.write('command_line = %s\n' % config_params['command_line'])
    outfile.write('num_iterations = %d\n' % config_params['num_iterations'])
    outfile.write('start_iteration = %d\n' % config_params['start_iteration'])
    outfile.write('output_dir = %s\n' % config_params['output_dir'])
    outfile.write('cache_dir = %s\n' % config_params['cache_dir'])
    outfile.write('tmp_dir = %s\n' % config_params['tmp_dir'])
    outfile.write('dbfile_name = %s\n' % config_params['dbfile_name'])
    outfile.write('use_multiprocessing = %s\n' % str(config_params['multiprocessing']))
    outfile.write('case_sensitive = %s\n' % str(config_params['case_sensitive']))
    if config_params['num_cores'] is None:
        outfile.write('num_cores =\n')
    else:
        outfile.write('num_cores = %d\n' % config_params['num_cores'])

    outfile.write('stats_frequency = %d\n' % config_params['stats_freq'])
    outfile.write('result_frequency = %d\n' % config_params['result_freq'])
    outfile.write('debug_frequency = %d\n' % config_params['debug_freq'])
    outfile.write('postadjust = %s\n' % str(config_params['postadjust']))
    outfile.write('add_fuzz = %s\n' % str(config_params['add_fuzz']))
    outfile.write('num_clusters = %d\n' % config_params['num_clusters'])
    outfile.write('random_seed = %s\n' % strparam(config_params['random_seed']))
    outfile.write('log_subresults = %s\n' % str(config_params['log_subresults']))

    # compatibility
    outfile.write('organism_code = %s\n' % str(config_params['organism_code']))
    outfile.write('use_operons = %s\n' % str(config_params['use_operons']))
    outfile.write('use_string = %s\n' % str(config_params['use_string']))
    outfile.write('checkratios = %s\n' % str(config_params['checkratios']))


def write_membership_settings(outfile, config_params):
    """Writes the Membership section of the final.ini file"""
    outfile.write('\n[Membership]\n')
    outfile.write('probability_row_change = %f\n' % config_params['memb.prob_row_change'])
    outfile.write('probability_column_change = %f\n' % config_params['memb.prob_col_change'])
    outfile.write('max_changes_per_row = %d\n' % config_params['memb.max_changes_per_row'])
    outfile.write('max_changes_per_column = %d\n' % config_params['memb.max_changes_per_col'])
    outfile.write('min_cluster_rows_allowed = %d\n' % config_params['memb.min_cluster_rows_allowed'])
    outfile.write('max_cluster_rows_allowed = %d\n' % config_params['memb.max_cluster_rows_allowed'])
    outfile.write('clusters_per_row = %d\n' % config_params['memb.clusters_per_row'])
    outfile.write('clusters_per_column = %d\n' % config_params['memb.clusters_per_col'])


def write_section(outfile, section, settings):
    """Writes a scoring function specific section to the final.ini file"""
    outfile.write('\n[%s]\n' % section)
    for setting, param in settings.iteritems():
        if setting == 'scaling':
            if param[0] == 'scaling_const':
                outfile.write('scaling_const = %f\n' % param[1])
            elif param[0] == 'scaling_rvec':
                outfile.write('scaling_rvec = %s\n' % param[1])
        else:
            outfile.write('%s = %s\n' % (setting, str(param)))
