#!/usr/bin/env python
# vi: sw=4 ts=4 et:
"""cmonkey_ensemble.py - cMonkey frontend for ensemble runs

This was created to keep ensemble specific aspects out of the core
codebase and tailor behavior easily

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import os.path
import cmonkey.cmonkey_run as cmr
import cmonkey.config as conf
import cmonkey.util as util
import logging
import cmonkey.scoring as scoring
import random

def arg_ext(argparser):
    """additional parameters added for ensembles"""
    argparser.add_argument('--ensemble_run_id', type=int, default=None)

BGORDER = [None, 0, 1, 2, 3, 4, 5]

if __name__ == '__main__':
    """process configuration"""
    args, params, ratios = conf.setup(arg_ext)
    if args.ensemble_run_id is not None:
        # setup for ensemble run
        params['random_seed'] = args.ensemble_run_id
        random.seed(params['random_seed'])
        util.r_set_seed(params['random_seed'])

        # these currently are experimental and very eco-centric
        params['num_clusters'] = random.randint(150, 550)

        # MEME parameters
        maxw = random.randint(12, 30)  # meme parameter
        mmotifs = random.randint(1, 3)
        motif_upstream_scan = (ramdom.randint(-50, 0), random.randint(150, 250))
        motif_upstream_search = (ramdom.randint(-20, 0), random.randint(100, 200))
        string_weight = random.uniform(0.0, 1.0) * 0.5 + 0.2
        operon_weight = random.uniform(0.0, 1.0) * 0.5 + 0.2

        # 50% probability for each network to be included in scoring
        params['use_string'] = random.randint(0, 1) == 1
        params['use_operons'] = random.randint(0, 1) == 1

        # background order None -> no background order
        bgorder = BGORDER[random.randint(0, 6)

    cmonkey_run = cmr.CMonkeyRun(ratios, params)
    cmonkey_run.run()
    cmonkey_run.cleanup()
