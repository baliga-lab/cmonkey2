'''
wrapper to run cMonkey Python with MY DATA
Created on Mar 22, 2012
modified from the original human.py wrapper
@author / modifier: frank schmitz, sbri, 2012
'''

from Environment.Wrappers.FileSystem import AutoInitFileStructure       # init my FileStructure

LOG_FORMAT = '%(asctime)s %(levelname)-8s %(message)s'
import logging
logging.basicConfig(format=LOG_FORMAT,level=logging.DEBUG)
logging.debug("starting cMonkey Python")


import util
import datamatrix_addon_FS as dm
import cmonkey_run
import copy
from datetime import date

import Kapsel as Kap
import FS_dataload as dataload
import stringDB_addon_FS
import network_addon_FS as nw
import motif_FS as motif_FS
import meme_addon_FS as meme_FS
import scoring
import membership as memb
import microarray
import pprint as pp
import json
import os
import numpy as np

from BioBabelFish.Servers.MonkeyPythonServers.MPS_HumanThesaurus import GetThesaurus_HSA_MP, BabelThesaurus
from BioBabelFish.Servers.MonkeyPythonServers.MPS_HumanNetwork import GetNetwork_HSA_MP, BabelNetwork
import Environment.Variables.BioInfo as EVB

RUG_PROPS = ['Pam', 'LPS', 'TNFa', 'Poly']
NUM_CLUSTERS = 50  # 133
ROW_WEIGHT = 6.0
NUM_ITERATIONS = 2000
MULTIPROCESSING = True
NETWORK_WEIGHT = 0.5
MAX_MOTIF_WIDTH = 12
STATS_FREQ = 10
RESULT_FREQ = 10

alphabet = EVB.Nalphabet        # the sequence alphabet to use, for Markov Model
alphabet_replacement = EVB.Nalphabet_replacement

class FScMonkeyRun(cmonkey_run.CMonkeyRun):

    def __init__(self, organism_code, ratio_matrix, num_clusters=None):
        cmonkey_run.CMonkeyRun.__init__(self, organism_code, ratio_matrix, num_clusters)

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
        logging.info("\x1b[31mMain:\t\x1b[0m# CLUSTERS: %d", self['num_clusters'])

        # defaults
        self.row_seeder = memb.make_kmeans_row_seeder(num_clusters)
        self.column_seeder = microarray.seed_column_members
        self['row_scaling'] = 6.0
        self['string_file'] = None
        self['cache_dir'] = 'cache'
        self['output_dir'] = 'out'
        self['start_iteration'] = 1
        self['num_iterations'] = 2000
        self['multiprocessing'] = True

        # membership default parameters
        self['memb.min_cluster_rows_allowed'] = 10
        self['memb.max_cluster_rows_allowed'] = 70

        self['sequence_types'] = ['Promoter', '3pUTR']
        self['search_distances'] = {'Promoter': (-1000, 200), '3pUTR': (0, 500)}
        # used for background distribution and MAST
        self['scan_distances'] = {'Promoter': (-2000, 750), '3pUTR': (0, 750)}
        logging.info("\x1b[31mMain:\t\x1b[0mcM object initialized")


        today = date.today()
        self.CHECKPOINT_INTERVAL = None
        self.__checkpoint_basename = "cmonkey-checkpoint-%s-%d%d%d" % (
            organism_code, today.year, today.month, today.day)

    def initKapsel(self):
        self.__organism.networks()
        logging.info("\x1b[31mKapsel:\t\x1b[0minitialized Kapsel %s" % (self.__organism.Kapsel_name))
    def make_Kapsel(self):
        """returns a human organism object"""
        OrgKapsel = Kap.Kapsel("hsa", "human", "Mensch in Kapsel")
        # add Network Factory to List
        OrgKapsel.put_thesaurus(GetThesaurus_HSA_MP(self['ThesaurusClass']).getBrachiosaurus())
        OrgKapsel.put_ReSaurus(GetThesaurus_HSA_MP(self['ThesaurusClass']).getOutDict())
        OrgKapsel.addNWfactory(stringDB_addon_FS.Babel_network_factory1_FS(self['NetworkClass'], NETWORK_WEIGHT), NETWORK_WEIGHT, nwfactname = "human STRING Db")
        # add a sequence environment, based on RSAT info
        OrgKapsel.addRSATenvir(self['RSATDir'], self['RSATfeaturefile'], self['RSATUTRfile'])
        OrgKapsel.addRSAT_feature("CDS", self['search_distances']["Promoter"], self['scan_distances']["Promoter"])            # add the feature from CDS, for promoter
        OrgKapsel.addRSAT_UTR("3'UTR", self['search_distances']["3pUTR"], self['scan_distances']["3pUTR"])              # add the feature 3'UTR
        OrgKapsel.CollectSequences("repeat_masked")         # grabs the seqs from the DNA contig files
                                                            # either normal or repeat_masked
                                                            # and stores it within Einheit and Sequence
        self.__organism = OrgKapsel
        return



    def organism(self):
        if self.__organism == None:
            self.__organism = self.make_hsa()
        return self.__organism
    def get_organism(self):
        return self.__organism



    def make_row_scoring(self):
        """returns the row scoring function"""
        row_scoring = microarray.RowScoringFunction(
            self.membership(), self.ratio_matrix,
            scaling_func=lambda iteration: self['row_scaling'],
            config_params=self.config_params)

        background_file_prom = meme_FS.global_background_file_FS(
            self.organism(), self.ratio_matrix.row_names(), 'Promoter',
            alphabet, alphabet_replacement,
            use_revcomp=True, resultsDir = "/".join([self['meme_dir'], 'backgroundFiles']))
        meme_suite = meme_FS.MemeSuite481(background_file=background_file_prom,
                                          remove_tempfiles=True,
                                          resultsDir = self['meme_dir'])
        """
        background_file_p3utr = meme_FS.global_background_file_FS(
            self.organism(), self.ratio_matrix.row_names(), '3pUTR',
            alphabet, alphabet_replacement,
            use_revcomp=True)
        meme_suite_p3utr = meme_FS.MemeSuite430(background_file=background_file_p3utr)
        """
        
        sequence_filters = [
            motif_FS.unique_filter,
            motif_FS.get_remove_low_complexity_filter(meme_suite),
            motif_FS.get_remove_atgs_filter(self['search_distances']['Promoter'])]


        sequence_filters = [motif_FS.get_remove_atgs_filter(self['search_distances']['Promoter'])]


        motif_scaling_fun = scoring.get_default_motif_scaling(self['num_iterations'])
        motif_scoring = motif_FS.MemeScoringFunction(
            self.organism(),
            self.membership(),
            self.ratio_matrix,
            meme_suite,
            seqtype='Promoter',
            sequence_filters=sequence_filters,
            pvalue_filter=motif_FS.MinPValueFilter(-9),
            scaling_func=motif_scaling_fun,
            update_in_iteration=scoring.schedule(11, 3),
            motif_in_iteration=scoring.schedule(10, 100),
            config_params=self.config_params)
        """
        motif_scoring_3pUTR = motif_FS.MemeScoringFunction(
            self.organism(),
            self.membership(),
            self.ratio_matrix,
            meme_suite,
            seqtype='3pUTR',
            sequence_filters=sequence_filters,
            pvalue_filter=motif_FS.MinPValueFilter(-9),
            scaling_func=motif_scaling_fun,
            run_in_iteration=scoring.schedule(100, 50),
            config_params=self.config_params)
        motif_combiner = scoring.ScoringFunctionCombiner(
            self.membership(),
            [motif_scoring, motif_scoring_3pUTR],
            scaling_func=lambda iteration: 0.5)
        """

        '''
        weeder_scoring = motif.WeederScoringFunction(
            self.organism(), self.membership(), self.ratio_matrix,
            meme_suite_p3utr, '3pUTR',
            pvalue_filter=motif.MinPValueFilter(-20.0),
            scaling_func=lambda iteration: 0.0,
            run_in_iteration=scoring.default_motif_iterations,
            config_params=self.config_params)

        motif_combiner = scoring.ScoringFunctionCombiner(
            self.membership(),
            [motif_scoring, weeder_scoring],
            scaling_func=lambda iteration: 0.5)
        '''

        network_scaling_fun = scoring.get_default_network_scaling(self['num_iterations'])
        network_scoring = nw.ScoringFunction(self.organism(),
                                             self.membership(),
                                             self.ratio_matrix,
                                             scaling_func=network_scaling_fun,
                                             run_in_iteration=scoring.schedule(1, 7),
                                             config_params=self.config_params)

        row_scoring_functions = [row_scoring, motif_scoring, network_scoring]
        #row_scoring_functions = [row_scoring]
        return scoring.ScoringFunctionCombiner(self.membership(),
                                               row_scoring_functions,
                                               log_subresults=True)
    def __make_dirs_if_needed(self):
        output_dir = self['output_dir']
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)

        cache_dir = self['cache_dir']
        if not os.path.exists(cache_dir):
            os.mkdir(cache_dir)

        meme_dir = self['meme_dir']
        if not os.path.exists(meme_dir):
            os.mkdir(meme_dir)


        # write the normalized ratio matrix for stats and visualization
        #if not os.path.exists(output_dir + '/ratios.tsv'):
        self.ratio_matrix.write_tsv_file_FS(output_dir + '/ratios.tsv', self.__organism.get_ReSaurus()['EID2GS'])



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

    def run_iterations(self, row_scoring, col_scoring):
        def residual_for(row_names, column_names):
            if len(column_names) <= 1 or len(row_names) <= 1:
                return 1.0
            else:
                matrix = self.ratio_matrix.submatrix_by_name(row_names, column_names)
                return matrix.residual()
        self.report_params()
        output_dir = self['output_dir']

        for iteration in range(self['start_iteration'],
                               self['num_iterations'] + 1):
            logging.info("\x1b[31mIteration:\t\x1b[0m # %d", iteration)
            iteration_result = {'iteration': iteration}
            exp_iteration_result = {'iteration': iteration}
            """
            """
            self.membership().update(self.ratio_matrix,
                                     row_scoring.compute(iteration_result),
                                     col_scoring.compute(iteration_result),
                                     iteration, self['num_iterations'])
            #pp.pprint (self.membership().get_rowassignment())

            if iteration > 0 and self.CHECKPOINT_INTERVAL and iteration % self.CHECKPOINT_INTERVAL == 0:
                self.save_checkpoint_data(iteration, row_scoring, col_scoring)

            if iteration == 1 or (iteration % RESULT_FREQ == 0):
                # Write a snapshot
                iteration_result['columns'] = {}
                iteration_result['rows'] = {}
                iteration_result['residuals'] = {}
                temprows = {}
                tempcols = {}
                for cluster in range(1, self['num_clusters'] + 1):
                    column_names = self.membership().columns_for_cluster(cluster)
                    row_names = self.membership().rows_for_cluster(cluster)
                    iteration_result['columns'][cluster] = column_names
                    iteration_result['rows'][cluster] = self.TranslateMembership(
                                                        row_names)
                    residual = residual_for(row_names, column_names)
                    iteration_result['residuals'][cluster] = residual
                    temprows[cluster] = row_names 


                # write results
                with open('%s/%d-results.json' % (output_dir, iteration), 'w') as outfile:
                    outfile.write(json.dumps(iteration_result))

                iteration_result['rows'] = copy.deepcopy(temprows)
                temprows[cluster] = None
            if iteration == 1 or (iteration % STATS_FREQ == 0):
                # write stats for this iteration
                residuals = []
                cluster_stats = {}
                for cluster in range(1, self['num_clusters'] + 1):
                    row_names = iteration_result['rows'][cluster]
                    column_names = iteration_result['columns'][cluster]
                    residual = residual_for(row_names, column_names)
                    residuals.append(residual)
                    
                    cluster_stats[cluster] = {'num_rows': len(row_names),
                                              'num_columns': len(column_names),
                                              'residual': residual
                                              }
                stats = {'cluster': cluster_stats, 'median_residual': np.median(residuals) }
                with open('%s/%d-stats.json' % (output_dir, iteration), 'w') as outfile:
                    outfile.write(json.dumps(stats))
            run_infos = row_scoring.run_logs()
            logging.info("Writing run map for this iteration")
            #print run_infos

        print "Done !!!!"

    def TranslateMembership(self, memberlist):
        
        newmemberL = []
        for member in memberlist:
            try:
                newmemberL.append(self.__organism.get_ReSaurus()['EID2GS'][member])
            except KeyError:
                newmemberL.append(member)
        return newmemberL
                
        
        
        

if __name__ == '__main__':
    print('cMonkey (Python port) (c) 2011, Institute for Systems Biology')
    print('running a modified version by Frank Schmitz, SBRI, 2012 - ')
    print('This program is licensed under the General Public License V3.')
    print('See README and LICENSE for details.\n')

    
    
    Dirs = AutoInitFileStructure()
    DataBase = "".join([Dirs.Resources, "experimental/cMonkey/MonkeyPython/"])
    CHECKPOINT_FILE = "".join([DataBase, "temp/Checkpoint.txt"])
    MatrixFile = "".join([DataBase, "ratiofile.txt"])
    CONTROLS_FILE =  "".join([DataBase, "controlsFile.txt"])
    RUG_FILE =  "".join([DataBase, "rugfile.txt"])

    cmonkey_run = FScMonkeyRun('hsa', dataload.read_matrix(MatrixFile, RUG_FILE, RUG_PROPS))
    cmonkey_run['thesaurus_filetype'] = "RSAT"
    cmonkey_run['RSATDir'] = Dirs.RSATHumanPath
    cmonkey_run['RSATfeaturefile'] = "".join([Dirs.RSATHumanPath, "cds.tab"])
    cmonkey_run['RSATUTRfile'] = "".join([Dirs.RSATHumanPath, "utr.tab"])
    cmonkey_run['output_dir'] = "".join([DataBase, "out"])
    cmonkey_run['cache_dir'] = "".join([DataBase, "CacheDir/"])
    cmonkey_run['meme_dir'] = "".join([DataBase, "meme_results"])
    cmonkey_run['ThesaurusClass'] = BabelThesaurus(Dirs, "EID", "hGS")
    cmonkey_run['NetworkClass'] = BabelNetwork(Dirs, [9606])

    
    cmonkey_run.make_Kapsel()                                      # generate Organism - human

    
    cmonkey_run.ratio_matrix.translateRowNames(cmonkey_run.get_organism().thesaurus())
                                                                # translate network node names 
    for Network in cmonkey_run.organism().networks(): Network.ShrinkNetwork(cmonkey_run.ratio_matrix.row_names())
                                                                # reduce network according to present Genes
    cmonkey_run.initKapsel()                                    # initialize Kapsel pre-maturely

    cmonkey_run.run()
    
    