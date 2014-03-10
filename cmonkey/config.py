"""config.py - cMonkey configuration module

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import os
from cmonkey.schedule import make_schedule

def set_config(params, config):
    def set_scaling(section):
        try:
            params['scaling'][section] = ('scaling_const', config.getfloat(section, 'scaling_const'))
            return
        except:
            pass
        try:
            params['scaling'][section] = ('scaling_rvec', config.get(section, 'scaling_rvec'))
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
    params['out_database'] = os.path.join(params['output_dir'],
                                               config.get("General", "dbfile_name"))
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

    params['sequence_types'] = config.get('Motifs', 'sequence_types').split(',')
    params['search_distances'] = {}
    params['scan_distances'] = {}
    for seqtype in params['sequence_types']:
        cat = "SequenceType-%s" % seqtype
        params['search_distances'][seqtype] = tuple(
            map(int, config.get(cat, 'search_distance').split(',')))
        params['scan_distances'][seqtype] = tuple(
            map(int, config.get(cat, 'scan_distance').split(',')))

    params['schedule'] = {}
    params['schedule']['Rows'] = make_schedule(config.get("Rows", "schedule"))
    params['schedule']['Columns'] = make_schedule(config.get("Columns", "schedule"))
    params['schedule']['MEME'] = make_schedule(config.get("MEME", "schedule"))
    params['schedule']['Motifs'] = make_schedule(config.get("Motifs", "schedule"))
    params['schedule']['Networks'] = make_schedule(config.get("Networks", "schedule"))

    params['stats_freq'] = config.getint('General', 'stats_frequency')
    params['result_freq'] = config.getint('General', 'result_frequency')

    # parse the scalings
    params['scaling'] = {}
    set_scaling('Motifs')
    set_scaling('Rows')
    set_scaling('Networks')

    try:
        params['nmotifs_rvec'] = config.get('MEME', 'nmotifs_rvec')
    except:
        raise Exception("no setting found to retrieve the MEME nmotifs function")


