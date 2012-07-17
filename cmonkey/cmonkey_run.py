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
from datetime import date, datetime
import json
import numpy as np
import gc
import sizes

KEGG_FILE_PATH = 'testdata/KEGG_taxonomy'
GO_FILE_PATH = 'testdata/proteome2taxid'
RSAT_BASE_URL = 'http://rsat.ccb.sickkids.ca'
COG_WHOG_URL = 'ftp://ftp.ncbi.nih.gov/pub/COG/COG/whog'
CACHE_DIR = 'cache'

LOG_FORMAT = '%(asctime)s %(levelname)-8s %(message)s'
STATS_FREQ = 10
RESULT_FREQ = 10

class CMonkeyRun:
    def __init__(self, organism_code, ratio_matrix, num_clusters=None):
        logging.basicConfig(format=LOG_FORMAT,
                            datefmt='%Y-%m-%d %H:%M:%S',
                            level=logging.DEBUG)
        self.__membership = None
        self.__organism = None
        self.config_params = {}
        self.ratio_matrix = ratio_matrix.sorted_by_row_name()

        # membership update default parameters
        # these come first, since a lot depends on clustering numbers
        self['memb.clusters_per_row'] = 2
        if num_clusters == None:
            num_clusters = int(round(self.ratio_matrix.num_rows() *
                                     self['memb.clusters_per_row'] / 20.0))
        self['memb.clusters_per_col'] = int(round(num_clusters * 2.0 / 3.0))
        self['memb.prob_row_change'] = 0.5
        self['memb.prob_col_change'] = 1.0
        self['memb.max_changes_per_row'] = 1
        self['memb.max_changes_per_col'] = 5

        self['organism_code'] = organism_code
        self['num_clusters'] = num_clusters
        logging.info("# CLUSTERS: %d", self['num_clusters'])

        # defaults
        self.row_seeder = memb.make_kmeans_row_seeder(num_clusters)
        self.column_seeder = microarray.seed_column_members
        self['row_scaling'] =  6.0
        self['string_file'] = None
        self['cache_dir'] = CACHE_DIR
        self['output_dir'] = 'out'
        self['start_iteration'] = 1
        self['num_iterations'] = 2000
        self['multiprocessing'] = True
        # Quantile normalization is false by default in cMonkey-R
        self['quantile_normalize'] = True

        # used to select sequences and MEME
        self['sequence_types'] = ['upstream']
        self['search_distances'] = {'upstream': (-20, 150)}
        # used for background distribution and MAST
        self['scan_distances'] = {'upstream': (-30, 250)}

        # membership default parameters
        self['memb.min_cluster_rows_allowed'] = 3
        self['memb.max_cluster_rows_allowed'] = 70

        today = date.today()
        self.CHECKPOINT_INTERVAL = None
        self.__checkpoint_basename = "cmonkey-checkpoint-%s-%d%d%d" % (
            organism_code, today.year, today.month, today.day)

    def report_params(self):
        logging.info('cmonkey_run config_params:')
        for param,value in self.config_params.items():
            logging.info('%s=%s' %(param,str(value)))

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
        self.row_scoring = row_scoring

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
            scaling_func=motif_scaling_fun,
            num_motif_func=motif.default_nmotif_fun,
            #update_in_iteration=scoring.schedule(601, 3),
            #motif_in_iteration=scoring.schedule(600, 100),
            update_in_iteration=scoring.schedule(100, 10),
            motif_in_iteration=scoring.schedule(100, 100),
            config_params=self.config_params)
        self.motif_scoring = motif_scoring

        network_scaling_fun = scoring.get_default_network_scaling(self['num_iterations'])
        network_scoring = nw.ScoringFunction(self.organism(),
                                             self.membership(),
                                             self.ratio_matrix,
                                             scaling_func=network_scaling_fun,
                                             run_in_iteration=scoring.schedule(1, 7),
                                             config_params=self.config_params)
        self.network_scoring = network_scoring

        row_scoring_functions = [row_scoring, motif_scoring, network_scoring]
        return scoring.ScoringFunctionCombiner(self.membership(),
                                               row_scoring_functions,
                                               config_params=self.config_params,
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
        self.__make_dirs_if_needed()
        row_scoring = self.make_row_scoring()
        col_scoring = self.make_column_scoring()
        self.run_iterations(row_scoring, col_scoring)

    def run_from_checkpoint(self,checkpoint_filename):
        row_scoring = self.make_row_scoring()
        col_scoring = self.make_column_scoring()
        self.__make_dirs_if_needed()
        self.init_from_checkpoint(checkpoint_filename, row_scoring, col_scoring)
        self.run_iterations(row_scoring, col_scoring)

    def residual_for(self, row_names, column_names):
        if len(column_names) <= 1 or len(row_names) <= 1:
            return 1.0
        else:
            matrix = self.ratio_matrix.submatrix_by_name(row_names, column_names)
            return matrix.residual()


    def write_results(self, iteration_result):
        # Write a snapshot
        iteration = iteration_result['iteration']
        iteration_result['columns'] = {}
        iteration_result['rows'] = {}
        iteration_result['residuals'] = {}
        for cluster in range(1, self['num_clusters'] + 1):
            column_names = self.membership().columns_for_cluster(cluster)
            row_names = self.membership().rows_for_cluster(cluster)
            iteration_result['columns'][cluster] = column_names
            iteration_result['rows'][cluster] = row_names
            residual = self.residual_for(row_names, column_names)
            iteration_result['residuals'][cluster] = residual

        # write results
        with open('%s/%d-results.json' % (self['output_dir'], iteration), 'w') as outfile:
            outfile.write(json.dumps(iteration_result))

    def write_stats(self, iteration_result):
        # write stats for this iteration
        iteration = iteration_result['iteration']
        residuals = []
        cluster_stats = {}
        network_scores = iteration_result['networks']
        if 'motif-pvalue' in iteration_result:
            motif_pvalue = iteration_result['motif-pvalue']
        else:
            motif_pvalue = 0.0

        if 'fuzzy-coeff' in iteration_result:
            fuzzy_coeff = iteration_result['fuzzy-coeff']
        else:
            fuzzy_coeff = 0.0

        for cluster in range(1, self['num_clusters'] + 1):
            row_names = iteration_result['rows'][cluster]
            column_names = iteration_result['columns'][cluster]
            residual = self.residual_for(row_names, column_names)
            residuals.append(residual)
            cluster_stats[cluster] = {'num_rows': len(row_names),
                                      'num_columns': len(column_names),
                                      'residual': residual }
        stats = {'cluster': cluster_stats, 'median_residual': np.median(residuals),
                 'motif-pvalue': motif_pvalue, 'network-scores': network_scores,
                 'fuzzy-coeff': fuzzy_coeff}
        with open('%s/%d-stats.json' % (self['output_dir'], iteration), 'w') as outfile:
            try:
                outfile.write(json.dumps(stats))
            except:
                logging.error("Could not write stats - probably non-serializable values found")
                # print stats object, likely there is something that is not serializable
                print stats

    def write_runlog(self, row_scoring, iteration):
        logging.info("Writing run map for this iteration")
        run_infos = [run_log.to_json() for run_log in row_scoring.run_logs()]
        with open('%s/%d-runlog.json' % (self['output_dir'], iteration), 'w') as outfile:
            try:
                outfile.write(json.dumps(run_infos))
            except:
                logging.error("Could not run map - probably non-serializable values found")
                # print run_infos object, likely there is something that is not serializable
                print run_infos

    def write_start_info(self):
        start_info = { 'start_time': str(datetime.now()),
                       'num_iterations': self['num_iterations'],
                       'organism-code': self.organism().code,
                       'species': self.organism().species(),
                       'num_rows': self.ratio_matrix.num_rows(),
                       'num_columns': self.ratio_matrix.num_columns()
                       }
        with open('%s/start.json' % self['output_dir'], 'w') as outfile:
            outfile.write(json.dumps(start_info))

    def write_finish_info(self):
        finish_info = { 'finish_time': str(datetime.now()) }
        with open('%s/finish.json' % self['output_dir'], 'w') as outfile:
            outfile.write(json.dumps(finish_info))


    def run_iterations(self, row_scoring, col_scoring):
        self.report_params()
        self.write_start_info()

        for iteration in range(self['start_iteration'],
                               self['num_iterations'] + 1):
            logging.info("Iteration # %d", iteration)
            iteration_result = {'iteration': iteration}

            rscores = row_scoring.compute(iteration_result)
            start_time = util.current_millis()
            cscores = col_scoring.compute(iteration_result)
            elapsed = util.current_millis() - start_time
            logging.info("computed column_scores in %f s.", elapsed / 1000.0)

            self.membership().update(self.ratio_matrix, rscores, cscores,
                                     self['num_iterations'], iteration_result)

            if iteration > 0 and self.CHECKPOINT_INTERVAL and iteration % self.CHECKPOINT_INTERVAL == 0:
                self.save_checkpoint_data(iteration, row_scoring, col_scoring)
            mean_net_score = 0.0
            mean_mot_pvalue = 0.0
            if 'networks' in iteration_result.keys():
                mean_net_score = iteration_result['networks']
            mean_mot_pvalue = "NA"
            if 'motif-pvalue' in iteration_result.keys():
                mean_mot_pvalue = ""
                mean_mot_pvalues = iteration_result['motif-pvalue']
                mean_mot_pvalue = ""
                for seqtype in mean_mot_pvalues.keys():
                    mean_mot_pvalue = mean_mot_pvalue + (" '%s' = %f" % (seqtype, mean_mot_pvalues[seqtype]))
                    
            logging.info('mean net = %s | mean mot = %s', str(mean_net_score), mean_mot_pvalue)

            if iteration == 1 or (iteration % RESULT_FREQ == 0):
                self.write_results(iteration_result)

            if iteration == 1 or (iteration % STATS_FREQ == 0):
                self.write_stats(iteration_result)
                # run infos should be written with the same frequency as stats
                self.write_runlog(row_scoring, iteration)

            gc.collect()
            print "# ROW SCORING: ", sizes.asizeof(self.row_scoring)
            print "# MOT SCORING: ", sizes.asizeof(self.motif_scoring)
            print "# NET SCORING: ", sizes.asizeof(self.network_scoring)
            print "# COL SCORING: ", sizes.asizeof(col_scoring)
            print "# MEMBERSHIP: ", sizes.asizeof(self.membership())

        logging.info("Postprocessing: Adjusting the clusters....")
        self.membership().postadjust()
        iteration = self['num_iterations'] + 1
        iteration_result = {'iteration': iteration }
        logging.info("Adjusted. Now re-run scoring (iteration: %d)", iteration_result['iteration'])
        row_scoring.compute_force(iteration_result)
        self.write_results(iteration_result)
        self.write_stats(iteration_result)
        self.write_finish_info()
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
