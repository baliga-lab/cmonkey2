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
import gzip
import sqlite3
from decimal import Decimal

KEGG_FILE_PATH = 'testdata/KEGG_taxonomy'
GO_FILE_PATH = 'testdata/proteome2taxid'
RSAT_BASE_URL = 'http://rsat.ccb.sickkids.ca'
COG_WHOG_URL = 'ftp://ftp.ncbi.nih.gov/pub/COG/COG/whog'
STRING_URL_PATTERN = "http://networks.systemsbiology.net/string9/%s.gz"
CACHE_DIR = 'cache'

LOG_FORMAT = '%(asctime)s %(levelname)-8s %(message)s'

# these parameters have a strong impact on the size of the result file
# Please set it to something >= 10 (which can be in the 10s of gigabytes)
STATS_FREQ = 50
RESULT_FREQ = 50


class CMonkeyRun:
    def __init__(self, organism_code, ratio_matrix,
                 string_file=None,
                 num_clusters=None,
                 use_operons=True):
        logging.basicConfig(format=LOG_FORMAT,
                            datefmt='%Y-%m-%d %H:%M:%S',
                            level=logging.DEBUG,
                            filename='cmonkey.log')
        self.__membership = None
        self.__organism = None
        self.config_params = {}
        self.ratio_matrix = ratio_matrix

        # membership update default parameters
        # these come first, since a lot depends on clustering numbers
        self['memb.clusters_per_row'] = 2
        if num_clusters == None:
            num_clusters = int(round(self.ratio_matrix.num_rows *
                                     self['memb.clusters_per_row'] / 20.0))
        self['memb.clusters_per_col'] = int(round(num_clusters * 2.0 / 3.0))
        #self['memb.clusters_per_col'] = int(round(num_clusters / 2.0))
        self['memb.prob_row_change'] = 0.5
        self['memb.prob_col_change'] = 1.0
        self['memb.max_changes_per_row'] = 1
        self['memb.max_changes_per_col'] = 5

        self['organism_code'] = organism_code
        self['num_clusters'] = num_clusters
        self['use_operons'] = use_operons
        logging.info("# CLUSTERS: %d", self['num_clusters'])
        logging.info("use operons: %d", self['use_operons'])

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
        self['string_file'] = string_file
        self['out_database'] = self['output_dir'] + '/cmonkey_run.db'

        self['keep_memeout'] = False

        today = date.today()
        self['checkpoint_interval'] = 100
        self.__checkpoint_basename = "cmonkey-checkpoint-%s-%d%d%d" % (
            organism_code, today.year, today.month, today.day)

    def __dbconn(self, isolation_level='DEFERRED'):
        """returns an autocommit database connection"""
        conn = sqlite3.connect(self['out_database'])
        conn.isolation_level = isolation_level
        return conn

    def __create_output_database(self):
        conn = self.__dbconn()
        c = conn.cursor()
        # these are the tables for storing cmonkey run information.
        # run information
        c.execute('''create table run_infos (start_time timestamp,
                     finish_time timestamp,
                     num_iterations int, last_iteration int,
                     organism text, species text, num_rows int,
                     num_columns int, num_clusters int)''')

        # stats tables
        # Note: there is some redundancy with the result tables here.
        # ----- I measured the cost for creating those on the fly and
        #       it is more expensive
        #       than I expected, so I left the tables in-place
        c.execute('''create table iteration_stats (iteration int,
                     median_residual decimal,
                     fuzzy_coeff decimal)''')
        c.execute('''create table cluster_stats (iteration int, cluster int,
                     num_rows int, num_cols int, residual decimal)''')
        c.execute('''create table network_stats (iteration int, network text,
                     score decimal)''')
        c.execute('''create table motif_stats (iteration int, seqtype text,
                     pval decimal)''')
        c.execute('''create table row_names (order_num int, name text)''')
        c.execute('''create table column_names (order_num int, name text)''')

        # result tables
        c.execute('''create table row_members (iteration int, cluster int,
                     order_num int)''')
        c.execute('''create table column_members (iteration int, cluster int,
                     order_num int)''')
        c.execute('''create table cluster_residuals (iteration int,
                     cluster int, residual decimal)''')

        # motif results: TODO: we might want to have the motif scoring function
        # ------------- write it
        # in case you are wondering about the redundant iteration field here -
        # it allows for much faster database access when selecting by iteration
        c.execute('''create table motif_infos (iteration int, cluster int,
                     seqtype text, motif_num int, evalue decimal)''')
        c.execute('''create table motif_pssm_rows (motif_info_id int,
                     iteration int, row int, a decimal, c decimal, g decimal,
                     t decimal)''')

        c.execute('''create table motif_annotations (motif_info_id int,
                     iteration int, gene_num int,
                     position int, reverse boolean, pvalue decimal)''')
        c.execute('''create table motif_pvalues (iteration int, cluster int,
                     gene_num int, pvalue decimal)''')
        c.execute('''create index if not exists colmemb_iter_index
                     on column_members (iteration)''')
        c.execute('''create index if not exists rowmemb_iter_index
                     on row_members (iteration)''')
        c.execute('''create index if not exists clustresid_iter_index
                     on cluster_residuals (iteration)''')
        logging.info("created output database schema")

        # all cluster members are stored relative to the base ratio matrix
        for index in xrange(len(self.ratio_matrix.row_names)):
            c.execute('''insert into row_names (order_num, name) values
                         (?,?)''',
                      (index, self.ratio_matrix.row_names[index]))
        for index in xrange(len(self.ratio_matrix.column_names)):
            c.execute('''insert into column_names (order_num, name) values
                         (?,?)''',
                      (index, self.ratio_matrix.column_names[index]))
        logging.info("added row and column names to output database")
        conn.commit()
        c.close()
        conn.close()

    def report_params(self):
        logging.info('cmonkey_run config_params:')
        for param, value in self.config_params.items():
            logging.info('%s=%s' % (param, str(value)))

    def __getitem__(self, key):
        return self.config_params[key]

    def __setitem__(self, key, value):
        self.config_params[key] = value

    def __make_membership(self):
        """returns the seeded membership on demand"""
        return memb.ClusterMembership.create(
            self.ratio_matrix,
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

        motif_scaling_fun = scoring.get_default_motif_scaling(
            self['num_iterations'])
        motif_scoring = motif.MemeScoringFunction(
            self.organism(),
            self.membership(),
            self.ratio_matrix,
            meme_suite,
            sequence_filters=sequence_filters,
            scaling_func=motif_scaling_fun,
            num_motif_func=motif.default_nmotif_fun,
            update_in_iteration=scoring.schedule(601, 3),
            motif_in_iteration=scoring.schedule(600, 100),
            #update_in_iteration=scoring.schedule(100, 10),
            #motif_in_iteration=scoring.schedule(100, 100),
            config_params=self.config_params)
        self.motif_scoring = motif_scoring

        network_scaling_fun = scoring.get_default_network_scaling(
            self['num_iterations'])
        network_scoring = nw.ScoringFunction(
            self.organism(),
            self.membership(),
            self.ratio_matrix,
            scaling_func=network_scaling_fun,
            run_in_iteration=scoring.schedule(1, 7),
            config_params=self.config_params)
        self.network_scoring = network_scoring

        row_scoring_functions = [row_scoring, motif_scoring, network_scoring]
        return scoring.ScoringFunctionCombiner(
            self.membership(),
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
        keggfile = util.read_dfile(KEGG_FILE_PATH, comment='#')
        gofile = util.read_dfile(GO_FILE_PATH)
        rsatdb = rsat.RsatDatabase(RSAT_BASE_URL, self['cache_dir'])
        mo_db = microbes_online.MicrobesOnline()
        stringfile = self.config_params['string_file']
        kegg_mapper = org.make_kegg_code_mapper(keggfile)
        rsat_mapper = org.make_rsat_organism_mapper(rsatdb)

        # automatically download STRING file
        if stringfile == None:
            rsat_info = rsat_mapper(kegg_mapper(self['organism_code']))
            ncbi_code = rsat_info.taxonomy_id
            logging.info("NCBI CODE IS: %s", ncbi_code)
            url = STRING_URL_PATTERN % ncbi_code
            stringfile = "%s/%s.gz" % (self['cache_dir'], ncbi_code)
            self['string_file'] = stringfile
            logging.info("Automatically using STRING file in '%s'", stringfile)
            util.get_url_cached(url, stringfile)

        nw_factories = []
        if stringfile != None:
            nw_factories.append(stringdb.get_network_factory2(
                    self['organism_code'], stringfile, 0.5))
        else:
            logging.warn("no STRING file specified !")

        if self['use_operons']:
            logging.info('adding operon network factory')
            nw_factories.append(microbes_online.get_network_factory(
                    mo_db, max_operon_size=self.ratio_matrix.num_rows / 20,
                    weight=0.5))
            

        org_factory = org.MicrobeFactory(kegg_mapper,
                                         rsat_mapper,
                                         org.make_go_taxonomy_mapper(gofile),
                                         mo_db,
                                         nw_factories)
        return org_factory.create(self['organism_code'],
                                  self['search_distances'],
                                  self['scan_distances'])

    def __make_dirs_if_needed(self):
        output_dir = self['output_dir']
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)

        cache_dir = self['cache_dir']
        if not os.path.exists(cache_dir):
            os.mkdir(cache_dir)

    def __clear_output_dir(self):
        output_dir = self['output_dir']
        if os.path.exists(output_dir):
            outfiles = os.listdir(output_dir)
            for filename in outfiles:
                os.remove('/'.join([output_dir, filename]))

    def run(self):
        self.__make_dirs_if_needed()
        self.__clear_output_dir()
        self.__create_output_database()
        # write the normalized ratio matrix for stats and visualization
        output_dir = self['output_dir']
        if not os.path.exists(output_dir + '/ratios.tsv'):
            self.ratio_matrix.write_tsv_file(output_dir + '/ratios.tsv')

        row_scoring = self.make_row_scoring()
        col_scoring = self.make_column_scoring()
        self.run_iterations(row_scoring, col_scoring)

    def run_from_checkpoint(self, checkpoint_filename):
        row_scoring = self.make_row_scoring()
        col_scoring = self.make_column_scoring()
        self.__make_dirs_if_needed()
        self.init_from_checkpoint(checkpoint_filename, row_scoring,
                                  col_scoring)
        self.run_iterations(row_scoring, col_scoring)

    def residual_for(self, row_names, column_names):
        if len(column_names) <= 1 or len(row_names) <= 1:
            return 1.0
        else:
            matrix = self.ratio_matrix.submatrix_by_name(row_names,
                                                         column_names)
            return matrix.residual()

    def write_results(self, iteration_result, compressed=True):
        """write iteration results to database"""
        iteration = iteration_result['iteration']
        conn = self.__dbconn()
        with conn:
            for cluster in range(1, self['num_clusters'] + 1):
                column_names = self.membership().columns_for_cluster(cluster)
                for order_num in self.ratio_matrix.column_indexes(column_names):
                    conn.execute('''insert into column_members (iteration,cluster,order_num)
                                    values (?,?,?)''', (iteration, cluster, order_num))

                row_names = self.membership().rows_for_cluster(cluster)
                for order_num in self.ratio_matrix.row_indexes(row_names):
                    conn.execute('''insert into row_members (iteration,cluster,order_num)
                                    values (?,?,?)''', (iteration, cluster, order_num))
                try:
                    residual = self.residual_for(row_names, column_names)
                    conn.execute('''insert into cluster_residuals (iteration,cluster,residual)
                               values (?,?,?)''', (iteration, cluster, residual))
                except:
                    # apparently computing the mean residual led to a numpy masked
                    # value. We set it to 1.0 to avoid crashing out
                    conn.execute('''insert into cluster_residuals (iteration,cluster,residual)
                               values (?,?,?)''', (iteration, cluster, 1.0))

        # write motif infos: TODO: we might want the motif scoring function writing
        # this part
        if 'motifs' in iteration_result:
            motifs = iteration_result['motifs']
            with conn:
                c = conn.cursor()
                for seqtype in motifs:
                    for cluster in motifs[seqtype]:
                        motif_infos = motifs[seqtype][cluster]['motif-info']
                        for motif_info in motif_infos:
                            c.execute('''insert into motif_infos (iteration,cluster,seqtype,motif_num,evalue)
                                        values (?,?,?,?,?)''',
                                        (iteration, cluster, seqtype, motif_info['motif_num'],
                                        motif_info['evalue']))
                            motif_info_id = c.lastrowid
                            pssm_rows = motif_info['pssm']
                            for row in xrange(len(pssm_rows)):
                                pssm_row = pssm_rows[row]
                                conn.execute('''insert into motif_pssm_rows (motif_info_id,iteration,row,a,c,g,t)
                                                values (?,?,?,?,?,?,?)''',
                                            (motif_info_id, iteration, row, pssm_row[0], pssm_row[1],
                                            pssm_row[2], pssm_row[3]))
                            annotations = motif_info['annotations']
                            for annotation in annotations:
                                gene_num = self.gene_indexes[annotation['gene']]
                                conn.execute('''insert into motif_annotations (motif_info_id,
                                                iteration,gene_num,
                                                position,reverse,pvalue) values (?,?,?,?,?,?)''',
                                        (motif_info_id, iteration, gene_num,
                                        annotation['position'],
                                        annotation['reverse'], annotation['pvalue']))

                        pvalues = motifs[seqtype][cluster]['pvalues']
                        for gene in pvalues:
                            gene_num = self.gene_indexes[gene]
                            conn.execute('''insert into motif_pvalues (iteration,cluster,gene_num,pvalue)
                                            values (?,?,?,?)''',
                                        (iteration, cluster, 4711, pvalues[gene]))
                c.close()
            conn.close()

    def write_stats(self, iteration_result):
        # write stats for this iteration
        iteration = iteration_result['iteration']

        network_scores = iteration_result['networks']
        motif_pvalues = iteration_result['motif-pvalue'] if 'motif-pvalue' in iteration_result else {}
        fuzzy_coeff = iteration_result['fuzzy-coeff'] if 'fuzzy-coeff' in iteration_result else 0.0

        residuals = []
        conn = self.__dbconn()
        with conn:
            for cluster in range(1, self['num_clusters'] + 1):
                row_names = self.membership().rows_for_cluster(cluster)
                column_names = self.membership().columns_for_cluster(cluster)
                residual = self.residual_for(row_names, column_names)
                residuals.append(residual)
                try:
                    conn.execute('''insert into cluster_stats (iteration, cluster, num_rows,
                                    num_cols, residual) values (?,?,?,?,?)''',
                                 (iteration, cluster, len(row_names), len(column_names),
                                  residual))
                except:
                    # residual is messed up, insert with 1.0
                    logging.warn('STATS: residual was messed up, insert with 1.0')
                    conn.execute('''insert into cluster_stats (iteration, cluster, num_rows,
                                    num_cols, residual) values (?,?,?,?,?)''',
                                 (iteration, cluster, len(row_names), len(column_names),
                                  1.0))

            median_residual = np.median(residuals)
            try:
                conn.execute('''insert into iteration_stats (iteration, median_residual,
                                fuzzy_coeff) values (?,?,?)''',
                             (iteration, median_residual, fuzzy_coeff))
            except:
                logging.warn('STATS: median was messed up, insert with 1.0')
                conn.execute('''insert into iteration_stats (iteration, median_residual,
                                fuzzy_coeff) values (?,?,?)''',
                             (iteration, 1.0, fuzzy_coeff))

        with conn:
            for network, score in network_scores.items():
                conn.execute('''insert into network_stats (iteration, network, score)
                                values (?,?,?)''', (iteration, network, score))

        with conn:
            for seqtype, pval in motif_pvalues.items():
                conn.execute('''insert into motif_stats (iteration, seqtype, pval)
                                values (?,?,?)''', (iteration, seqtype, pval))
        conn.close()

    def write_start_info(self):
        conn = self.__dbconn()
        c = conn.cursor()
        with conn:
            conn.execute('''insert into run_infos (start_time, num_iterations, organism,
                            species, num_rows, num_columns, num_clusters) values (?,?,?,?,?,?,?)''',
                         (datetime.now(), self['num_iterations'], self.organism().code,
                          self.organism().species(), self.ratio_matrix.num_rows,
                          self.ratio_matrix.num_columns, self['num_clusters']))
        conn.close()

    def update_iteration(self, iteration):
        conn = self.__dbconn()
        with conn:
            conn.execute('''update run_infos set last_iteration = ?''', (iteration,))
        conn.close()

    def write_finish_info(self):
        conn = self.__dbconn()
        with conn:
            conn.execute('''update run_infos set finish_time = ?''', (datetime.now(),))
        conn.close()

    def run_iterations(self, row_scoring, col_scoring):
        self.report_params()
        self.write_start_info()
        thesaurus = self.organism().thesaurus()

        genes = [thesaurus[row_name] if row_name in thesaurus else row_name
                 for row_name in self.ratio_matrix.row_names]
        self.gene_indexes = {genes[index]: index
                             for index in xrange(len(genes))}
        # mark available for GC
        genes = None
        thesaurus = None
        membership = self.membership()

        for iteration in range(self['start_iteration'],
                               self['num_iterations'] + 1):
            logging.info("Iteration # %d", iteration)
            iteration_result = {'iteration': iteration}

            rscores = row_scoring.compute(iteration_result)
            start_time = util.current_millis()
            cscores = col_scoring.compute(iteration_result)
            elapsed = util.current_millis() - start_time
            logging.info("computed column_scores in %f s.", elapsed / 1000.0)

            membership.update(self.ratio_matrix, rscores, cscores,
                                     self['num_iterations'], iteration_result)
            membership.reseed_small_row_clusters(self.ratio_matrix.row_names)
            membership.reseed_small_column_clusters(self.ratio_matrix.column_names)

            if (iteration > 0 and self['checkpoint_interval']
                and iteration % self['checkpoint_interval'] == 0):
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
                self.update_iteration(iteration)

            gc.collect()
            #print "# ROW SCORING: ", sizes.asizeof(self.row_scoring)
            #print "# MOT SCORING: ", sizes.asizeof(self.motif_scoring)
            #print "# NET SCORING: ", sizes.asizeof(self.network_scoring)
            #print "# COL SCORING: ", sizes.asizeof(col_scoring)
            #print "# MEMBERSHIP: ", sizes.asizeof(membership)

        logging.info("Postprocessing: Adjusting the clusters....")
        membership.postadjust()
        iteration = self['num_iterations'] + 1
        iteration_result = {'iteration': iteration}
        logging.info("Adjusted. Now re-run scoring (iteration: %d)", iteration_result['iteration'])
        row_scoring.compute_force(iteration_result)
        self.write_results(iteration_result)
        self.write_stats(iteration_result)
        self.update_iteration(self['num_iterations'] + 1)
        self.write_finish_info()
        logging.info("Done !!!!")

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
