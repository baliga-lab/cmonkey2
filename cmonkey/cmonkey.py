# vi: sw=4 ts=4 et:
"""cmonkey.py - cMonkey top-level module

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import sys
import logging
import microbial_config as microbe
import tps_config as tps
import human
import json

CMONKEY_VERSION = '4.0'
#CHECKPOINT_INTERVAL = 3
CHECKPOINT_INTERVAL = 100
CHECKPOINT_FILE = None


def run_cmonkey(config):
    """init of the cMonkey system"""
    gene_scoring = config.row_scoring()
    cond_scoring = config.column_scoring()

    for iteration in range(config.start_iteration(),
                           config.num_iterations()):
        logging.info("Iteration # %d", iteration)
        iteration_result = {'iteration': iteration}
        config.membership().update(config.matrix(),
                                   gene_scoring.compute(iteration_result),
                                   cond_scoring.compute(iteration_result),
                                   iteration, config.num_iterations())
        if iteration > 0 and  iteration % CHECKPOINT_INTERVAL == 0:
            config.save_checkpoint_data(iteration)

        # Write a snapshot
        iteration_result['columns'] = {}
        iteration_result['rows'] = {}
        for cluster in range(1, config.membership().num_clusters() + 1):
            iteration_result['columns'][cluster] = config.membership().columns_for_cluster(cluster)
            iteration_result['rows'][cluster] = config.membership().rows_for_cluster(cluster)

        with open('out/%d-results.json' % iteration, 'w') as outfile:
            outfile.write(json.dumps(iteration_result))
    print "Done !!!!"
    print "cluster\t# rows"
    for cluster in range(1, config.membership().num_clusters() + 1):
        print "%d\t%d" % (cluster,
                          config.membership().num_row_members(cluster))
    #print membership


if __name__ == '__main__':
    print('cMonkey (Python port) (c) 2011, Institute for Systems Biology')
    print('This program is licensed under the General Public License V3.')
    print('See README and LICENSE for details.\n')
    if len(sys.argv) <= 2:
        print('Usage: ./run_cmonkey.sh <organism-code> <ratio-file> ' +
              '[checkpoint-file]')
    else:
        if len(sys.argv) > 3:
            CHECKPOINT_FILE = sys.argv[3]

        if sys.argv[1] == 'hsa':
            run_cmonkey(human.CMonkeyConfiguration.create(
                    sys.argv[2], checkpoint_file=CHECKPOINT_FILE))
        elif sys.argv[1] == 'tps':
            run_cmonkey(tps.CMonkeyConfiguration.create(
                    sys.argv[1], sys.argv[2],
                    checkpoint_file=CHECKPOINT_FILE))
        else:
            run_cmonkey(microbe.CMonkeyConfiguration.create(
                    sys.argv[1], sys.argv[2],
                    checkpoint_file=CHECKPOINT_FILE))
