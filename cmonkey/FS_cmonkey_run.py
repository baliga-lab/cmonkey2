'''
wrapper to run cMonkey Python with MY DATA
Created on Mar 22, 2012
modified from the original human.py wrapper
@author / modifier: frank schmitz, sbri, 2012
'''

from Environment.Wrappers.FileSystem import AutoInitFileStructure       # init my FileStructure

import util
from datetime import date
import Kapsel as Kap
import datamatrix as dm
import scoring
import logging
import microarray
import stringdb
import network as nw
import motif
import meme
import cmonkey_run
import pprint as pp


RUG_PROPS = ['Pam', 'LPS']
NUM_CLUSTERS = 50  # 133
ROW_WEIGHT = 6.0
NUM_ITERATIONS = 2000
MULTIPROCESSING = True
MAX_MOTIF_WIDTH = 12
LOG_FORMAT = '%(asctime)s %(levelname)-8s %(message)s'
STATS_FREQ = 10
RESULT_FREQ = 10



def read_rug(pred):
    """reads the rug file"""
    infile = util.DelimitedFile.read(RUG_FILE, sep='\t', has_header=False)      # Frank Schmitz changed separator to TAB delimited
    return list(set([row[0] for row in infile.lines() if pred(row)]))


def read_matrix(filename):
    """reads the data matrix from a file"""
    rug = read_rug(lambda row: row[1] in RUG_PROPS)
    columns_to_use = list(set(rug))

    # pass the column filter as the first filter to the DataMatrixFactory,
    # so normalization will be applied to the submatrix
    matrix_factory = dm.DataMatrixFactory([
            lambda matrix: matrix.submatrix_by_name(
                column_names=columns_to_use)])
    infile = util.DelimitedFile.read(filename, sep='\t', has_header=True,   #FS changed to TAB separator
                                     quote="\"")
    matrix = matrix_factory.create_from(infile)

    column_groups = {1: range(matrix.num_columns())}
    return matrix


class FScMonkeyRun(cmonkey_run.CMonkeyRun):

    def __init__(self, organism_code, ratio_matrix, num_clusters):
        cmonkey_run.CMonkeyRun.__init__(self, organism_code, ratio_matrix, num_clusters)
        self.__organism = None
        self.config_params = {'num_clusters' : NUM_CLUSTERS}
        #self.ratio_matrix = ratio_matrix.sorted_by_row_name()            #Frank Schmitz - I am not so sure...

        self['row_scaling'] = 6.0
        self['start_iteration'] = 0
        self['num_iterations'] = NUM_ITERATIONS
        self['multiprocessing'] = MULTIPROCESSING

        # used to select sequences and MEME
        self['sequence_types'] = ['Promoter', '3pUTR']
        self['search_distances'] = {'Promoter': (-1500, 500), '3pUTR': (0, 500)}
        # used for background distribution and MAST
        self['scan_distances'] = {'Promoter': (-2000, 750), '3pUTR': (0, 750)}



        # membership update default parameters
        self['memb.clusters_per_row'] = 2
        self['memb.clusters_per_col'] = int(round(num_clusters * 2.0 / 3.0))
        self['memb.prob_row_change'] = 0.5
        self['memb.prob_col_change'] = 1.0
        self['memb.max_changes_per_row'] = 1
        self['memb.max_changes_per_col'] = 5
        self["memb.min_cluster_rows_allowed"] = 5
        # motifing default parameters
        self['motif.min_cluster_rows_allowed'] = 3
        self['motif.max_cluster_rows_allowed'] = 70

        
        today = date.today()
        self.CHECKPOINT_INTERVAL = 100
        self.__checkpoint_basename = "cmonkey-checkpoint-%s-%d%d%d" % (
            organism_code, today.year, today.month, today.day)


        logging.info("\x1b[31mMain:\t\x1b[0mcM object initialized")

    def organism(self):
        if self.__organism == None:
            self.__organism = self.make_hsa()
        return self.__organism

    def get_organism(self):
        #if self.__organism == None:
        #    self.__organism = self.make_hsa()
        # here to test whether a Kapsel can be pickled...
        return self.__organism

    def initKapsel(self):
        self.__organism.networks()
        logging.info("\x1b[31mKapsel:\t\x1b[0minitialized Kapsel %s" % (self.__organism.Kapsel_name))

    def make_hsa(self):
        """returns a human organism object"""
        OrgKapsel = Kap.Kapsel("hsa", "human", "Mensch in Kapsel")
        # add Network Factory to List
        OrgKapsel.addThesaurusFile(self['thesaurus_file'], self['thesaurus_filetype'])
        OrgKapsel.addNWfactory(stringdb.get_network_factory2_FS(self['string_file'], 0.5, sep=" "), "human STRING Db")
        # add thesaurus file(s) to Kapsel

        # add a sequence environment, based on RSAT info
        OrgKapsel.addRSATenvir(self['RSATDir'], self['RSATfeaturefile'], self['RSATUTRfile'])
        OrgKapsel.addRSAT_feature("CDS", self['search_distances']["Promoter"], self['scan_distances']["Promoter"])            # add the feature from CDS, for promoter
        OrgKapsel.addRSAT_UTR("3'UTR", self['search_distances']["3pUTR"], self['scan_distances']["3pUTR"])              # add the feature 3'UTR
        OrgKapsel.CollectSequences("repeat_masked")        # grabs the seqs from the DNA contig files
                                                    # either normal or repeat_masked
                                                    # and stores it within Einheit and Sequence
        #OrgKapsel.dumpPickle("/Users/Frank/Desktop/KapselPickle.cMonkey")
        #done, return
        self.__organism = OrgKapsel
        return



    def make_row_scoring(self):
        """returns the row scoring function"""

        row_scoring = microarray.RowScoringFunction(
            self.membership(), self.ratio_matrix,
            scaling_func=lambda iteration: self['row_scaling'],
            config_params=self.config_params)




        logging.info("starting meme collections")

        background_file_prom = meme.global_background_file_FS(
            self.organism(), self.ratio_matrix.row_names(), 'Promoter',
            use_revcomp=True)


        meme_suite = meme.MemeSuite481(
            max_width=MAX_MOTIF_WIDTH,
            background_file=background_file_prom)
        logging.info("intermediate meme collections")
        
        sequence_filters = [
            motif.unique_filter,
            motif.get_remove_low_complexity_filter(meme_suite),
            motif.get_remove_atgs_filter(self['search_distances']['Promoter'])]

        sequence_filters = [motif.get_remove_atgs_filter(self['search_distances']['Promoter'])]


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
        return scoring.ScoringFunctionCombiner(
            self.membership(),
            [row_scoring, motif_scoring, network_scoring])


if __name__ == '__main__':
    print('cMonkey (Python port) (c) 2011, Institute for Systems Biology')
    print('running a modified version by Frank Schmitz, SBRI, 2012 - ')
    print('This program is licensed under the General Public License V3.')
    print('See README and LICENSE for details.\n')
    
    Dirs = AutoInitFileStructure()
    STRINGFILE = "".join([Dirs.StringPath, "HumanProteinLinks.txt"])
    THESAURUSFILE = "".join([Dirs.RSATHumanPath, "cds_names.tab"])
    DataBase = "".join([Dirs.expRespath, "cMonkey/MonkeyPython/"])
    CHECKPOINT_FILE = "".join([DataBase, "temp/Checkpoint.txt"])
    MatrixFile = "".join([DataBase, "ratiofile.txt"])
    CONTROLS_FILE =  "".join([DataBase, "controlsFile.txt"])
    RUG_FILE =  "".join([DataBase, "rugfile.txt"])

    cmonkey_run = FScMonkeyRun('hsa', read_matrix(MatrixFile), NUM_CLUSTERS)
    cmonkey_run['string_file'] = STRINGFILE
    cmonkey_run['thesaurus_file'] = THESAURUSFILE
    cmonkey_run['thesaurus_filetype'] = "RSAT"
    cmonkey_run['RSATDir'] = Dirs.RSATHumanPath
    cmonkey_run['RSATfeaturefile'] = "".join([Dirs.RSATHumanPath, "cds.tab"])
    cmonkey_run['RSATUTRfile'] = "".join([Dirs.RSATHumanPath, "utr.tab"])
    cmonkey_run['output_dir'] = "".join([DataBase, "out"])
    cmonkey_run['cache_dir'] = "".join([DataBase, "CacheDir/"])
    cmonkey_run.make_hsa()                                      # generate Organism - human
    cmonkey_run.ratio_matrix.translateRowNames(cmonkey_run.get_organism().thesaurus())
                                                                # translate network node names 
    for Network in cmonkey_run.organism().networks(): Network.ShrinkNetwork(cmonkey_run.ratio_matrix.row_names())
                                                                # reduce network according to present Genes
    cmonkey_run.initKapsel()                                    # initialize Kapsel pre-maturely

    cmonkey_run.run()
    
    