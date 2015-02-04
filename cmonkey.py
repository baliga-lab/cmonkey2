#!/usr/bin/env python
# vi: sw=4 ts=4 et:
"""cmonkey.py - cMonkey top-level module

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

if __name__ == '__main__':
    """process configuration"""
    args, params, ratios = conf.setup()
    
    cmonkey_run = cmr.CMonkeyRun(ratios, params)

    proceed = True

    if params['checkratios']:
        thesaurus = cmonkey_run.organism().thesaurus()
        logging.info("Checking the quality of the input matrix names...")
        found = [name for name in ratios.row_names if name in thesaurus]
        num_found = len(found)
        total = len(ratios.row_names)
        percent = (float(num_found) / float(total)) * 100.0
        proceed = percent > 50.0

    if not proceed:
        logging.error("# genes found: %d, # total: %d, %f %% - please check your ratios file",
                      num_found, total, percent)
    else:
        cmonkey_run.run()
    if not args.interactive:
        cmonkey_run.cleanup()
