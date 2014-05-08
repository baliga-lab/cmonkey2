# vi: sw=4 ts=4 et:
import re
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
import debug
import os
from datetime import date, datetime
import json
import numpy as np
import gc
import sizes
import gzip
import sqlite3
from decimal import Decimal
import cPickle
import bz2
import config

USER_KEGG_FILE_PATH = 'config/KEGG_taxonomy'
USER_GO_FILE_PATH = 'config/proteome2taxid'
SYSTEM_KEGG_FILE_PATH = '/etc/cmonkey-python/KEGG_taxonomy'
SYSTEM_GO_FILE_PATH = '/etc/cmonkey-python/proteome2taxid'

# pipeline paths
PIPELINE_USER_PATHS = {
    'default': 'config/default_pipeline.json',
    'rows': 'config/rows_pipeline.json',
    'rowsandmotifs': 'config/rows_and_motifs_pipeline.json',
    'rowsandnetworks': 'config/rows_and_networks_pipeline.json'
}
PIPELINE_SYSTEM_PATHS = {
    'default': '/etc/cmonkey-python/default_pipeline.json',
    'rows': '/etc/cmonkey-python/rows_pipeline.json',
    'rowsandmotifs': '/etc/cmonkey-python/rows_and_motifs_pipeline.json',
    'rowsandnetworks': '/etc/cmonkey-python/rows_and_networks_pipeline.json'
}

COG_WHOG_URL = 'ftp://ftp.ncbi.nih.gov/pub/COG/COG/whog'
STRING_URL_PATTERN = "http://networks.systemsbiology.net/string9/%s.gz"

# We support non-microbes the easy way for now, until we have a working
# database scheme
VERTEBRATES = {'hsa', 'mmu', 'rno'}

class CMonkeyRun:
    def __init__(self, ratios, args):
        self.__membership = None
        self.__organism = None
        self.config_params = args
        self.ratios = ratios
        self.row_seeder = memb.make_kmeans_row_seeder(args['num_clusters'])
        self.column_seeder = microarray.seed_column_members
        self.__conn = None

        today = date.today()
        logging.info("# clusters/row: %d", args['memb.clusters_per_row'])
        logging.info("# clusters/column: %d", args['memb.clusters_per_col'])
        logging.info("# CLUSTERS: %d", args['num_clusters'])
        logging.info("use operons: %d", args['use_operons'])

        if args['meme_version']:
            logging.info('using MEME version %s', args['meme_version'])
        else:
            logging.error('MEME not detected - please check')

    def cleanup(self):
        """cleanup this run object"""
        if self.__conn is not None:
            self.__conn.close()
            self.__conn = None

    def __dbconn(self):
        """Returns an autocommit database connection. We maintain a single database
        connection throughout the life of this run objec"""
        if self.__conn is None:
            self.__conn = sqlite3.connect(self['out_database'], 15, isolation_level='DEFERRED')
        return self.__conn

    def __create_output_database(self):
        conn = self.__dbconn()
        # these are the tables for storing cmonkey run information.
        # run information
        conn.execute('''create table run_infos (start_time timestamp,
                        finish_time timestamp,
                        num_iterations int, last_iteration int,
                        organism text, species text, num_rows int,
                        num_columns int, num_clusters int)''')

        # stats tables
        # Note: there is some redundancy with the result tables here.
        # ----- I measured the cost for creating those on the fly and
        #       it is more expensive
        #       than I expected, so I left the tables in-place
        conn.execute('''create table iteration_stats (iteration int,
                        median_residual decimal,
                        fuzzy_coeff decimal)''')
        conn.execute('''create table cluster_stats (iteration int, cluster int,
                        num_rows int, num_cols int, residual decimal)''')
        conn.execute('''create table network_stats (iteration int, network text,
                        score decimal)''')
        conn.execute('''create table motif_stats (iteration int, seqtype text,
                        pval decimal)''')
        conn.execute('''create table row_names (order_num int, name text)''')
        conn.execute('''create table column_names (order_num int, name text)''')

        # result tables
        conn.execute('''create table row_members (iteration int, cluster int,
                        order_num int)''')
        conn.execute('''create table column_members (iteration int, cluster int,
                        order_num int)''')
        conn.execute('''create table cluster_residuals (iteration int,
                        cluster int, residual decimal)''')

        # in case you are wondering about the redundant iteration field here -
        # it allows for much faster database access when selecting by iteration
        conn.execute('''create table motif_infos (iteration int, cluster int,
                        seqtype text, motif_num int, evalue decimal)''')
        conn.execute('''create table motif_pssm_rows (motif_info_id int,
                        iteration int, row int, a decimal, c decimal, g decimal,
                        t decimal)''')

        # Additional info: MEME generated top matching sites
        conn.execute('''create table meme_motif_sites (motif_info_id int,
                        seq_name text,
                        reverse boolean, start int, pvalue decimal,
                        flank_left text, seq text, flank_right text)''')

        conn.execute('''create table motif_annotations (motif_info_id int,
                        iteration int, gene_num int,
                        position int, reverse boolean, pvalue decimal)''')
        conn.execute('''create index if not exists colmemb_iter_index
                        on column_members (iteration)''')
        conn.execute('''create index if not exists rowmemb_iter_index
                        on row_members (iteration)''')
        conn.execute('''create index if not exists clustresid_iter_index
                        on cluster_residuals (iteration)''')
        logging.info("created output database schema")

        # all cluster members are stored relative to the base ratio matrix
        with conn:
            for index in xrange(len(self.ratios.row_names)):
                conn.execute('''insert into row_names (order_num, name) values
                                (?,?)''',
                             (index, self.ratios.row_names[index]))
            for index in xrange(len(self.ratios.column_names)):
                conn.execute('''insert into column_names (order_num, name) values
                                (?,?)''',
                             (index, self.ratios.column_names[index]))
        logging.info("added row and column names to output database")

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
        if 'random_seed' in self['debug']:
            util.r_set_seed(10)

        return memb.create_membership(self.ratios,
                                      self.row_seeder, self.column_seeder,
                                      self.config_params)

    def membership(self):
        if self.__membership is None:
            logging.info("creating and seeding memberships")
            self.__membership = self.__make_membership()

            # debug: write seed into an analytical file for iteration 0
            if 'random_seed' in self['debug']:
                conn = self.__dbconn()
                with conn:
                    self.write_memberships(conn, 0)
                # write complete result into a cmresults.tsv
                path =  os.path.join(self['output_dir'], 'cmresults-0000.tsv.bz2')
                with bz2.BZ2File(path, 'w') as outfile:
                    debug.write_iteration(conn, outfile, 0,
                                          self['num_clusters'], self['output_dir'])

        return self.__membership

    def organism(self):
        """returns the organism object to work on"""
        if self['dummy_organism']:
            self.__organism = org.DummyOrganism()
        elif self.__organism is None:
            self.__organism = self.make_organism()
        return self.__organism

    def make_organism(self):
        """returns the organism object to work on"""
        self.__make_dirs_if_needed()

        if os.path.exists(USER_KEGG_FILE_PATH):
            keggfile = util.read_dfile(USER_KEGG_FILE_PATH, comment='#')
        elif os.path.exists(SYSTEM_KEGG_FILE_PATH):
            keggfile = util.read_dfile(SYSTEM_KEGG_FILE_PATH, comment='#')
        else:
            raise Exception('KEGG file not found !!')

        if os.path.exists(USER_GO_FILE_PATH):
            gofile = util.read_dfile(USER_GO_FILE_PATH)
        elif os.path.exists(SYSTEM_GO_FILE_PATH):
            gofile = util.read_dfile(SYSTEM_GO_FILE_PATH)
        else:
            raise Exception('GO file not found !!')

        if self['rsat_dir']:
            if not self['rsat_organism']:
                raise Exception('override RSAT loading: please specify --rsat_organism')
            logging.info("using RSAT files for '%s'", self['rsat_organism'])
            rsatdb = rsat.RsatFiles(self['rsat_dir'], self['rsat_organism'], self['ncbi_code'])
        else:
            rsatdb = rsat.RsatDatabase(rsat.RSAT_BASE_URL, self['cache_dir'])

        if self['operon_file']:
            logging.info("using operon file at '%s'", self['operon_file'])
            mo_db = microbes_online.MicrobesOnlineOperonFile(self['operon_file'])
        else:
            logging.info("attempting automatic download of operons from Microbes Online")
            mo_db = microbes_online.MicrobesOnline(self['cache_dir'])

        stringfile = self['string_file']
        kegg_mapper = org.make_kegg_code_mapper(keggfile)
        ncbi_code = self['ncbi_code']
        nw_factories = []
        is_microbe = self['organism_code'] not in VERTEBRATES

        # do we use STRING ?
        if not self['nonetworks'] and self['use_string']:
            # download if not provided
            if stringfile is None:
                if ncbi_code is None:
                    rsat_info = org.RsatSpeciesInfo(rsatdb,
                                                    kegg_mapper(self['organism_code']),
                                                    self['rsat_organism'], None)
                    ncbi_code = rsat_info.taxonomy_id

                logging.info("NCBI CODE IS: %s", ncbi_code)
                url = STRING_URL_PATTERN % ncbi_code
                stringfile = "%s/%s.gz" % (self['cache_dir'], ncbi_code)
                self['string_file'] = stringfile
                logging.info("Automatically using STRING file in '%s'", stringfile)
                util.get_url_cached(url, stringfile)
            else:
                logging.info("Loading STRING file at '%s'", stringfile)

            # create and add network
            nw_factories.append(stringdb.get_network_factory(
                self['organism_code'], stringfile, 0.5))

        # do we use operons ?
        if not self['nonetworks'] and self['use_operons']:
            logging.info('adding operon network factory')
            nw_factories.append(microbes_online.get_network_factory(
                mo_db, max_operon_size=self.ratios.num_rows / 20,
                weight=0.5))

        if is_microbe:
            orgcode = self['organism_code']
            logging.info("Creating Microbe object for '%s'", orgcode)
            keggorg = kegg_mapper(orgcode)
            rsat_info = org.RsatSpeciesInfo(rsatdb, keggorg, self['rsat_organism'],
                                            self['ncbi_code'])
            gotax = org.make_go_taxonomy_mapper(gofile)(rsat_info.go_species())
            return org.Microbe(orgcode, keggorg, rsat_info, gotax, mo_db, nw_factories,
                               self['search_distances'], self['scan_distances'],
                               self['use_operons'], self.ratios)

    def __make_dirs_if_needed(self):
        logging.info('creating aux directories')
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

    def __check_parameters(self):
        """ensure that we all required parameters before we start running"""
        PARAM_NAMES = ['num_iterations', 'start_iteration', 'multiprocessing',
                       'quantile_normalize',
                       'memb.min_cluster_rows_allowed', 'memb.max_cluster_rows_allowed',
                       'memb.prob_row_change', 'memb.prob_col_change',
                       'memb.max_changes_per_row', 'memb.max_changes_per_col',
                       'sequence_types', 'search_distances', 'scan_distances']

        for param in PARAM_NAMES:
            if param not in self.config_params:
                raise Exception("required parameter not found in config: '%s'" % param)

    
    def __setup_pipeline(self):
        """Reading pipeline setup
        By default, this uses the default pipelines defined in config
        The default pipeline can be modified by
        1. nomotifs switch
        2. nonetworks switch
        
        User-defined pipelines can be provided using a JSON file, which is
        specified using the --pipeline switch on the command line
        """
        pipeline_id = 'default'
        if self['nonetworks'] and self['nomotifs']:
            pipeline_id = 'rows'
        elif self['nonetworks']:
            pipeline_id = 'rowsandmotifs'
        elif self['nomotifs']:
            pipeline_id = 'rowsandnetworks'

        if self['pipeline_file']:
            pipeline_file = self['pipeline_file']
            if os.path.exists(pipeline_file):
                with open(pipeline_file) as infile:
                    self['pipeline'] = json.load(infile)
            else:
                raise Exception("Pipeline file '%s' does not exist" % pipeline_file)
        else:
            if os.path.exists(PIPELINE_USER_PATHS[pipeline_id]):
                with open(PIPELINE_USER_PATHS[pipeline_id]) as infile:
                    self['pipeline'] = json.load(infile)

        # TODO: for now, we always assume the top level of row scoring is a combiner
        class_ = get_function_class(self['pipeline']['row-scoring']['function'])
        if class_.__name__ == 'ScoringFunctionCombiner':
            funs = [get_function_class(fun['function'])(self.organism(),
                                                       self.membership(),
                                                       self.ratios,
                                                       self.config_params)
                    for fun in self['pipeline']['row-scoring']['args']['functions']]
            row_scoring = class_(self.organism(), self.membership(), funs, self.config_params)
        else:
            raise Exception('Row scoring top level must be ScoringFunctionCombiner')

        # column scoring
        class_ = get_function_class(self['pipeline']['column-scoring']['function'])
        col_scoring = class_(self.organism(), self.membership(), self.ratios,
                             config_params=self.config_params)
        return row_scoring, col_scoring

    def prepare_run(self, check_params=True):
        """Setup output directories and scoring functions for the scoring.
        Separating setup and actual run facilitates testing"""
        self['dummy_organism'] = (self['organism_code'] is None and
                                  self['nonetworks'] and self['nomotifs'])
        if check_params:
            self.__check_parameters()
        self.__make_dirs_if_needed()
        self.__clear_output_dir()
        self.__create_output_database()
        # write the normalized ratio matrix for stats and visualization
        output_dir = self['output_dir']
        if not os.path.exists(output_dir + '/ratios.tsv'):
            self.ratios.write_tsv_file(output_dir + '/ratios.tsv')

        # gene index map is used for writing statistics
        thesaurus = self.organism().thesaurus()
        genes = [thesaurus[row_name] if row_name in thesaurus else row_name
                 for row_name in self.ratios.row_names]
        self.gene_indexes = {genes[index]: index
                             for index in xrange(len(genes))}
        row_scoring, col_scoring = self.__setup_pipeline()
        row_scoring.check_requirements()
        col_scoring.check_requirements()

        config.write_setup(self.config_params)
        return row_scoring, col_scoring

    def run(self):
        row_scoring, col_scoring = self.prepare_run()
        self.run_iterations(row_scoring, col_scoring)

    def residual_for(self, row_names, column_names):
        if len(column_names) <= 1 or len(row_names) <= 1:
            return 1.0
        else:
            matrix = self.ratios.submatrix_by_name(row_names, column_names)
            return matrix.residual()

    def write_memberships(self, conn, iteration):
        for cluster in range(1, self['num_clusters'] + 1):
            column_names = self.membership().columns_for_cluster(cluster)
            for order_num in self.ratios.column_indexes_for(column_names):
                conn.execute('''insert into column_members (iteration,cluster,order_num)
                                values (?,?,?)''', (iteration, cluster, order_num))

            row_names = self.membership().rows_for_cluster(cluster)
            for order_num in self.ratios.row_indexes_for(row_names):
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

    def write_results(self, iteration_result):
        """write iteration results to database"""
        iteration = iteration_result['iteration']
        conn = self.__dbconn()
        with conn:
            self.write_memberships(conn, iteration)

        if 'motifs' in iteration_result:
            motifs = iteration_result['motifs']
            with conn:
                for seqtype in motifs:
                    for cluster in motifs[seqtype]:
                        motif_infos = motifs[seqtype][cluster]['motif-info']
                        for motif_info in motif_infos:
                            c = conn.cursor()
                            c.execute('''insert into motif_infos (iteration,cluster,seqtype,motif_num,evalue)
                                        values (?,?,?,?,?)''',
                                      (iteration, cluster, seqtype, motif_info['motif_num'],
                                       motif_info['evalue']))
                            motif_info_id = c.lastrowid
                            c.close()
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

                            sites = motif_info['sites']
                            if len(sites) > 0 and isinstance(sites[0], tuple):
                                for seqname, strand, start, pval, flank_left, seq, flank_right in sites:
                                    conn.execute('''insert into meme_motif_sites (motif_info_id, seq_name, reverse, start, pvalue, flank_left, seq, flank_right)
                                                    values (?,?,?,?,?,?,?,?)''',
                                                 (motif_info_id, seqname, strand == '-',
                                                  start, pval, flank_left, seq,
                                                  flank_right))

    def write_stats(self, iteration_result):
        # write stats for this iteration
        iteration = iteration_result['iteration']

        network_scores = iteration_result['networks'] if 'networks' in iteration_result else {}
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

    def write_start_info(self):
        conn = self.__dbconn()
        with conn:
            conn.execute('''insert into run_infos (start_time, num_iterations, organism,
                            species, num_rows, num_columns, num_clusters) values (?,?,?,?,?,?,?)''',
                         (datetime.now(), self['num_iterations'], self.organism().code,
                          self.organism().species(), self.ratios.num_rows,
                          self.ratios.num_columns, self['num_clusters']))

    def update_iteration(self, iteration):
        conn = self.__dbconn()
        with conn:
            conn.execute('''update run_infos set last_iteration = ?''', (iteration,))

    def write_finish_info(self):
        conn = self.__dbconn()
        with conn:
            conn.execute('''update run_infos set finish_time = ?''', (datetime.now(),))

    def combined_rscores_pickle_path(self):
        return "%s/combined_rscores_last.pkl" % self.config_params['output_dir']

    def run_iteration(self, row_scoring, col_scoring, iteration):
        logging.info("Iteration # %d", iteration)
        iteration_result = {'iteration': iteration}
        rscores = row_scoring.compute(iteration_result)
        start_time = util.current_millis()
        cscores = col_scoring.compute(iteration_result)
        elapsed = util.current_millis() - start_time
        if elapsed > 0.0001:
            logging.info("computed column_scores in %f s.", elapsed / 1000.0)

        self.membership().update(self.ratios, rscores, cscores,
                                 self['num_iterations'], iteration_result)

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

        if iteration == 1 or (iteration % self['result_freq'] == 0):
            self.write_results(iteration_result)

        if iteration == 1 or (iteration % self['stats_freq'] == 0):
            self.write_stats(iteration_result)
            self.update_iteration(iteration)

        if 'dump_results' in self['debug'] and (iteration == 1 or
                                                (iteration % self['debug_freq'] == 0)):
            # write complete result into a cmresults.tsv
            conn = self.__dbconn()
            path =  os.path.join(self['output_dir'], 'cmresults-%04d.tsv.bz2' % iteration)
            with bz2.BZ2File(path, 'w') as outfile:
                debug.write_iteration(conn, outfile, iteration,
                                      self['num_clusters'], self['output_dir'])

    def write_mem_profile(self, outfile, row_scoring, col_scoring, iteration):
        membsize = sizes.asizeof(self.membership()) / 1000000.0
        orgsize = sizes.asizeof(self.organism()) / 1000000.0
        colsize = sizes.asizeof(col_scoring) / 1000000.0
        funs = row_scoring.scoring_functions
        rowsize = sizes.asizeof(funs[0]) / 1000000.0
        netsize = sizes.asizeof(funs[1]) / 1000000.0
        motsize = sizes.asizeof(funs[2]) / 1000000.0
        outfile.write('%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n' % (iteration, membsize, orgsize, colsize, rowsize, netsize, motsize))

    def run_iterations(self, row_scoring, col_scoring):
        self.report_params()
        self.write_start_info()
        if 'profile_mem' in self['debug']:
            with open(os.path.join(self['output_dir'], 'memprofile.tsv'), 'w') as outfile:
                outfile.write('Iteration\tMembership\tOrganism\tCol\tRow\tNetwork\tMotif\n')

        for iteration in range(self['start_iteration'],
                               self['num_iterations'] + 1):
            start_time = util.current_millis()
            self.run_iteration(row_scoring, col_scoring, iteration)
            # garbage collection after everything in iteration went out of scope
            gc.collect()
            elapsed = util.current_millis() - start_time
            logging.info("performed iteration %d in %f s.", iteration, elapsed / 1000.0)
            
            if 'profile_mem' in self['debug'] and (iteration == 1 or iteration % 100 == 0):
                with open(os.path.join(self['output_dir'], 'memprofile.tsv'), 'a') as outfile:
                    self.write_mem_profile(outfile, row_scoring, col_scoring, iteration)


        """run post processing after the last iteration. We store the results in
        num_iterations + 1 to have a clean separation"""
        if self['postadjust']:
            logging.info("Postprocessing: Adjusting the clusters....")
            # run combiner using the weights of the last iteration
            rscores = row_scoring.combine_cached(self['num_iterations'])
            rd_scores = memb.get_row_density_scores(self.membership(), rscores)
            logging.info("Recomputed combined + density scores.")
            memb.postadjust(self.membership(), rd_scores)
            logging.info("Adjusted. Now re-run scoring (iteration: %d)",
                         self['num_iterations'])

            iteration_result = {'iteration': self['num_iterations'] + 1}
            combined_scores = row_scoring.compute_force(iteration_result)

            # write the combined scores for benchmarking/diagnostics
            with open(self.combined_rscores_pickle_path(), 'w') as outfile:
                cPickle.dump(combined_scores, outfile)

            self.write_results(iteration_result)
            self.write_stats(iteration_result)
            self.update_iteration(iteration)

            if 'dump_results' in self['debug']:
                # write complete result into a cmresults.tsv
                conn = self.__dbconn()
                path =  os.path.join(self['output_dir'], 'cmresults-postproc.tsv.bz2')
                with bz2.BZ2File(path, 'w') as outfile:
                    debug.write_iteration(conn, outfile,
                                          self['num_iterations'] + 1,
                                          self['num_clusters'], self['output_dir'])


        self.write_finish_info()
        logging.info("Done !!!!")


def get_function_class(scorefun):
    modulepath = scorefun['module'].split('.')
    if len(modulepath) > 1:
        module = __import__(scorefun['module'], fromlist=[modulepath[1]])
    else:
        module = __import__(modulepath[0])
    return getattr(module, scorefun['class'])        

