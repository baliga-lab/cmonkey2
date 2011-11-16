"""cmonkey.py - cMonkey top-level module

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import sys
import logging
import microbial_config as microbe
import human

CMONKEY_VERSION = '4.0'


def run_cmonkey(config):
    """init of the cMonkey system"""
    gene_scoring = config.row_scoring()
    cond_scoring = config.column_scoring()

    for iteration in range(config.num_iterations()):
        logging.info("Iteration # %d", iteration)
        config.membership().update(config.matrix(),
                                   gene_scoring.compute(iteration),
                                   cond_scoring.compute(iteration),
                                   iteration, config.num_iterations())
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
        print('Usage: ./run_cmonkey.sh <organism-code> <ratio-file>')
    else:
        if sys.argv[1] == 'hsa':
            run_cmonkey(human.CMonkeyConfiguration(sys.argv[2]))
        else:
            run_cmonkey(microbe.CMonkeyConfiguration(sys.argv[1], sys.argv[2]))
