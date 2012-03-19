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


def run_cmonkey(config):
    """init of the cMonkey system"""
    gene_scoring = config.row_scoring()
    cond_scoring = config.column_scoring()
    output_dir   = config.output_directory()
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    for iteration in range(config.start_iteration(),
                           config.num_iterations()):
        logging.info("Iteration # %d", iteration)
        iteration_result = {'iteration': iteration}
        config.membership().update(config.matrix(),
                                   gene_scoring.compute(iteration_result),
                                   cond_scoring.compute(iteration_result),
                                   iteration, config.num_iterations())
        #if iteration > 0 and  iteration % CHECKPOINT_INTERVAL == 0:
        #    config.save_checkpoint_data(iteration)

        # Write a snapshot
        iteration_result['columns'] = {}
        iteration_result['rows'] = {}
        for cluster in range(1, config.membership().num_clusters() + 1):
            iteration_result['columns'][cluster] = config.membership().columns_for_cluster(cluster)
            iteration_result['rows'][cluster] = config.membership().rows_for_cluster(cluster)

        with open('%s/%d-results.json' % (output_dir, iteration), 'w') as outfile:
            outfile.write(json.dumps(iteration_result))
    print "Done !!!!"
    print "cluster\t# rows"
    for cluster in range(1, config.membership().num_clusters() + 1):
        print "%d\t%d" % (cluster,
                          config.membership().num_row_members(cluster))
    #print membership


if __name__ == '__main__':
    print('cMonkey (Python port) (c) 2011-2012, Institute for Systems Biology')
    print('This program is licensed under the General Public License V3.')
    print('See README and LICENSE for details.\n')
    if len(sys.argv) <= 3:
        print('Usage: ./run_cmonkey.sh <organism-code> <ratio-file> ' +
              '<string-file> [checkpoint-file]')
    else:
        if len(sys.argv) > 4:
            CHECKPOINT_FILE = sys.argv[4]

        matrix_factory = dm.DataMatrixFactory([dm.nochange_filter, dm.center_scale_filter])
        infile = util.DelimitedFile.read(sys.argv[2], has_header=True, quote='\"')
        matrix = matrix_factory.create_from(infile)
        cmonkey_run = cmonkey_run.CMonkeyRun(sys.argv[1], matrix, 43)
        cmonkey_run['string_file'] = sys.argv[3]
        cmonkey_run.run()
