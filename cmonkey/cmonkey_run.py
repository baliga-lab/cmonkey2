# vi: sw=4 ts=4 et:
import logging
import microarray
import membership as memb
import meme
import motif
import util
import rsat
import microbes_online
import organism as org
import scoring
import network as nw
import stringdb
import os
from datetime import date
import json
import numpy as np

KEGG_FILE_PATH = 'testdata/KEGG_taxonomy'
GO_FILE_PATH = 'testdata/proteome2taxid'
RSAT_BASE_URL = 'http://rsat.ccb.sickkids.ca'
COG_WHOG_URL = 'ftp://ftp.ncbi.nih.gov/pub/COG/COG/whog'
CACHE_DIR = 'cache'

LOG_FORMAT = '%(asctime)s %(levelname)-8s %(message)s'
STATS_FREQ = 10
RESULT_FREQ = 10

class CMonkeyRun:
    def __init__(self, organism_code, ratio_matrix, num_clusters):
        logging.basicConfig(format=LOG_FORMAT,
                            datefmt='%Y-%m-%d %H:%M:%S',
                            level=logging.DEBUG)
        self.__membership = None
        self.__organism = None
        self.config_params = {}
        self.ratio_matrix = ratio_matrix.sorted_by_row_name()

        self['organism_code'] = organism_code
        self['num_clusters'] = num_clusters

        # defaults
        self.row_seeder = memb.make_kmeans_row_seeder(num_clusters)
        self.column_seeder = microarray.seed_column_members
        self['row_scaling'] = 6.0
        self['string_file'] = None
        self['cache_dir'] = CACHE_DIR
        self['output_dir'] = 'out'
        self['start_iteration'] = 1
        self['num_iterations'] = 2000
        self['multiprocessing'] = True

        # used to select sequences and MEME
        self['sequence_types'] = ['upstream']
        self['search_distances'] = {'upstream': (-20, 150)}
        # used for background distribution and MAST
        self['scan_distances'] = {'upstream': (-30, 250)}

        # membership update default parameters
        self['memb.clusters_per_row'] = 2
        self['memb.clusters_per_col'] = int(round(num_clusters * 2.0 / 3.0))
        self['memb.prob_row_change'] = 0.5
        self['memb.prob_col_change'] = 1.0
        self['memb.max_changes_per_row'] = 1
        self['memb.max_changes_per_col'] = 5

        # motifing default parameters
        self['motif.min_cluster_rows_allowed'] = 3
        self['motif.max_cluster_rows_allowed'] = 70

        today = date.today()
        self.CHECKPOINT_INTERVAL = None
        self.__checkpoint_basename = "cmonkey-checkpoint-%s-%d%d%d" % (
            organism_code, today.year, today.month, today.day)


    def __getitem__(self, key):
        return self.config_params[key]

    def __setitem__(self, key, value):
        self.config_params[key] = value

    def __make_membership(self):
        """returns the seeded membership on demand"""
        return memb.ClusterMembership.create(
            self.ratio_matrix.sorted_by_row_name(),
            self.row_seeder, self.column_seeder,
            self.config_params)

    def make_column_scoring(self):
        """returns the column scoring function"""
        return scoring.ColumnScoringFunction(
            self.membership(), self.ratio_matrix,
                config_params=self.config_params)

    def make_row_scoring(self):
        """makes a row scoring function on demand"""
        # Default row scoring functions
        row_scoring = microarray.RowScoringFunction(
            self.membership(), self.ratio_matrix,
            scaling_func=lambda iteration: self['row_scaling'],
            config_params=self.config_params)

        meme_suite = meme.MemeSuite430()
        sequence_filters = [
            motif.unique_filter,
            motif.get_remove_low_complexity_filter(meme_suite),
            motif.get_remove_atgs_filter(self['search_distances']['upstream'])]

        motif_scaling_fun = scoring.get_default_motif_scaling(self['num_iterations'])
        motif_scoring = motif.MemeScoringFunction(
            self.organism(),
            self.membership(),
            self.ratio_matrix,
            meme_suite,
            sequence_filters=sequence_filters,
            pvalue_filter=motif.MinPValueFilter(-20.0),
            scaling_func=motif_scaling_fun,
            run_in_iteration=scoring.default_motif_iterations,
            config_params=self.config_params)

        network_scaling_fun = scoring.get_default_network_scaling(self['num_iterations'])
        network_scoring = nw.ScoringFunction(self.organism(),
                                             self.membership(),
                                             self.ratio_matrix,
                                             scaling_func=network_scaling_fun,
                                             run_in_iteration=scoring.default_network_iterations,
                                             config_params=self.config_params)

        row_scoring_functions = [row_scoring, motif_scoring, network_scoring]
        return scoring.ScoringFunctionCombiner(self.membership(),
                                               row_scoring_functions,
                                               log_subresults=True)

    def membership(self):
        if self.__membership == None:
            logging.info("creating and seeding memberships")
            self.__membership = self.__make_membership()
        return self.__membership

    def organism(self):
        """returns the organism object to work on"""
        if self.__organism == None:
            self.__organism = self.make_microbe()
        return self.__organism

    def make_microbe(self):
        """returns the organism object to work on"""
        keggfile = util.DelimitedFile.read(KEGG_FILE_PATH, comment='#')
        gofile = util.DelimitedFile.read(GO_FILE_PATH)
        rsatdb = rsat.RsatDatabase(RSAT_BASE_URL, self['cache_dir'])
        mo_db = microbes_online.MicrobesOnline()
        stringfile = self.config_params['string_file']

        nw_factories = []
        if stringfile != None:
            nw_factories.append(stringdb.get_network_factory2(stringfile, 0.5))
        else:
            logging.warn("no STRING file specified !")

        nw_factories.append(microbes_online.get_network_factory(
                mo_db, max_operon_size=self.ratio_matrix.num_rows() / 20, weight=0.5))

        org_factory = org.MicrobeFactory(org.make_kegg_code_mapper(keggfile),
                                         org.make_rsat_organism_mapper(rsatdb),
                                         org.make_go_taxonomy_mapper(gofile),
                                         mo_db,
                                         nw_factories)
        return org_factory.create(self['organism_code'], self['search_distances'],
                                  self['scan_distances'])

    def __make_dirs_if_needed(self):
        output_dir = self['output_dir']
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)

        cache_dir = self['cache_dir']
        if not os.path.exists(cache_dir):
            os.mkdir(cache_dir)

        # write the normalized ratio matrix for stats and visualization
        if not os.path.exists(output_dir + '/ratios.tsv'):
            self.ratio_matrix.write_tsv_file(output_dir + '/ratios.tsv')

    def run(self):
        row_scoring = self.make_row_scoring()
        col_scoring = self.make_column_scoring()
        self.__make_dirs_if_needed()
        self.run_iterations(row_scoring, col_scoring)

    def run_from_checkpoint(self,checkpoint_filename):
        row_scoring = self.make_row_scoring()
        col_scoring = self.make_column_scoring()
        self.__make_dirs_if_needed()
        self.init_from_checkpoint(checkpoint_filename, row_scoring, col_scoring)
        self.run_iterations(row_scoring, col_scoring)

    def run_iterations(self, row_scoring, col_scoring):
        output_dir = self['output_dir']

        for iteration in range(self['start_iteration'],
                               self['num_iterations'] + 1):
            logging.info("Iteration # %d", iteration)
            iteration_result = {'iteration': iteration}
            self.membership().update(self.ratio_matrix,
                                     row_scoring.compute(iteration_result),
                                     col_scoring.compute(iteration_result),
                                     iteration, self['num_iterations'])

            if iteration > 0 and self.CHECKPOINT_INTERVAL and iteration % self.CHECKPOINT_INTERVAL == 0:
                self.save_checkpoint_data(iteration, row_scoring, col_scoring)

            if iteration == 1 or (iteration % RESULT_FREQ == 0):
                # Write a snapshot
                iteration_result['columns'] = {}
                iteration_result['rows'] = {}
                for cluster in range(1, self['num_clusters'] + 1):
                    iteration_result['columns'][cluster] = self.membership().columns_for_cluster(cluster)
                    iteration_result['rows'][cluster] = self.membership().rows_for_cluster(cluster)

                # write results
                with open('%s/%d-results.json' % (output_dir, iteration), 'w') as outfile:
                    outfile.write(json.dumps(iteration_result))

            if iteration == 1 or (iteration % STATS_FREQ == 0):
                # write stats for this iteration
                residuals = []
                cluster_stats = {}
                for cluster in range(1, self['num_clusters'] + 1):
                    row_names = iteration_result['rows'][cluster]
                    column_names = iteration_result['columns'][cluster]
                    if len(column_names) <= 1 or len(row_names) <= 1:
                        residual = 1.0
                    else:
                        matrix = self.ratio_matrix.submatrix_by_name(row_names, column_names)
                        residual = matrix.residual()
                    residuals.append(residual)
                    cluster_stats[cluster] = {'num_rows': len(row_names),
                                              'num_columns': len(column_names),
                                              'residual': residual
                                              }
                stats = {'cluster': cluster_stats, 'median_residual': np.median(residuals) }
                with open('%s/%d-stats.json' % (output_dir, iteration), 'w') as outfile:
                    outfile.write(json.dumps(stats))

        print "Done !!!!"

    ############################################################
    ###### CHECKPOINTING
    ##############################

    def save_checkpoint_data(self, iteration, row_scoring, col_scoring):
        """save checkpoint data for the specified iteration"""
        with util.open_shelf("%s.%d" % (self.__checkpoint_basename,
                                        iteration)) as shelf:
            shelf['config'] = self.config_params
            shelf['iteration'] = iteration
            self.membership().store_checkpoint_data(shelf)
            row_scoring.store_checkpoint_data(shelf)
            col_scoring.store_checkpoint_data(shelf)

    def init_from_checkpoint(self, checkpoint_filename, row_scoring, col_scoring):
        """initialize this object from a checkpoint file"""
        logging.info("Continue run using checkpoint file '%s'",
                     checkpoint_filename)
        with util.open_shelf(checkpoint_filename) as shelf:
            self.config_params = shelf['config']
            self['start_iteration'] = shelf['iteration'] + 1
            self.__membership = memb.ClusterMembership.restore_from_checkpoint(
                self.config_params, shelf)
            row_scoring.restore_checkpoint_data(shelf)
            col_scoring.restore_checkpoint_data(shelf)
            #return row_scoring, col_scoring necessary??
