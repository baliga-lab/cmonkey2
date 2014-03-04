# vi: sw=4 ts=4 et:
"""motif.py - cMonkey motif related processing
This module captures the motif-specific scoring component
of cMonkey.

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import logging
import multiprocessing as mp
import numpy as np
import scoring
import datamatrix as dm
import weeder
import meme
import tempfile
import seqtools as st
import util
import os
import cPickle
import collections

ComputeScoreParams = collections.namedtuple('ComputeScoreParams',
                                            ['iteration',
                                             'cluster', 'feature_ids', 'seqs', 'used_seqs',
                                             'meme_runner', 'min_cluster_rows',
                                             'max_cluster_rows', 'num_motifs',
                                             'previous_motif_infos',
                                             'keep_memeout', 'outdir', 'num_iterations',
                                             'debug'])


def num_meme_motif_fun(params):
    """returns a function that returns the number of motifs to be returned by MEME"""
    return util.get_iter_fun(params, "nmotifs", params['num_iterations'])


# Applicable sequence filters
def unique_filter(seqs, feature_ids):
    """returns a map that contains only the keys that are in
    feature_ids and only contains unique sequences"""
    unique_seqs = {}
    for feature_id in feature_ids:
        if feature_id in seqs and seqs[feature_id] not in unique_seqs.values():
            unique_seqs[feature_id] = seqs[feature_id]
    return unique_seqs


def get_remove_low_complexity_filter(meme_suite):
    """Factory method that returns a low complexity filter"""
    def remove_low_complexity(seqs, feature_ids):
        """low-complexity filter that depends on meme"""
        return meme_suite.remove_low_complexity(seqs)
    return remove_low_complexity


def get_remove_atgs_filter(distance):
    """returns a remove ATG filter"""
    def remove_atgs_filter(seqs, feature_ids):
        """a filter removes the ATG's from the sequence, this
        just masks a window of 4 letters with N's"""
        for feature_id in seqs:
            chars = [c for c in seqs[feature_id]]
            chars[distance[1]:distance[1] + 4] = "NNNN"
            seqs[feature_id] = "".join(chars)
        return seqs
    return remove_atgs_filter


def compute_mean_score(pvalue_matrix, membership, organism):
    """cluster-specific mean scores"""
    if pvalue_matrix is None:
        return 0.0
    values = []
    pvalues = pvalue_matrix.values
    for cluster in xrange(1, membership.num_clusters() + 1):
        cluster_rows = membership.rows_for_cluster(cluster)
        row_indexes = pvalue_matrix.row_indexes_for(cluster_rows)
        for row in row_indexes:
            values.append(pvalues[row][cluster - 1])
    return np.median(values)

# Readonly structure to avoid passing it to the forked child processes for efficiency.
# non-serializable parameters go here, too
MOTIF_PARAMS = None
SEQUENCE_FILTERS = None
ORGANISM = None
MEMBERSIP = None


def pvalues2matrix(all_pvalues, num_clusters, gene_names, reverse_map):
    """converts a map from {cluster: {feature: pvalue}} to a scoring matrix
    """
    row_map = {gene: index for index, gene in enumerate(gene_names)}

    # convert remapped to an actual scoring matrix
    matrix = dm.DataMatrix(len(gene_names), num_clusters,
                           gene_names)
    mvalues = matrix.values
    for cluster, feature_pvals in all_pvalues.items():
        for feature_id, pval in feature_pvals.items():
            ridx = row_map[reverse_map[feature_id]]
            mvalues[ridx][cluster - 1] = pval

    matrix.apply_log()
    return matrix


class MotifScoringFunctionBase(scoring.ScoringFunctionBase):
    """Base class for motif scoring functions that use MEME"""

    def __init__(self, organism, membership, ratios,
                 meme_suite, seqtype,
                 sequence_filters=[],
                 scaling_func=None,
                 num_motif_func=None,
                 update_in_iteration=lambda iteration: True,
                 motif_in_iteration=lambda iteration: True,
                 config_params=None):
        """creates a ScoringFunction"""
        # run_in_iteration does not apply here, since we actually have
        # two schedules, motif_in_iteration and update_in_iteration here
        scoring.ScoringFunctionBase.__init__(self, organism, membership,
                                             ratios, scaling_func,
                                             schedule=None,
                                             config_params=config_params)
        # attributes accessible by subclasses
        self.meme_suite = meme_suite
        self.seqtype = seqtype
        self.update_in_iteration = update_in_iteration
        self.motif_in_iteration = motif_in_iteration
        self.num_motif_func = num_motif_func

        self.__sequence_filters = sequence_filters
        self.__last_motif_infos = None
        self.__last_iteration_result = {}
        self.all_pvalues = None
        self.last_result = None

        self.update_log = scoring.RunLog("motif-score-" + seqtype, config_params)
        self.motif_log = scoring.RunLog("motif-motif-" + seqtype, config_params)

        used_genes = sorted(ratios.row_names)
        self.used_seqs = organism.sequences_for_genes_scan(
            used_genes, seqtype=self.seqtype)

        logging.info("building reverse map...")
        start_time = util.current_millis()
        self.reverse_map = self.__build_reverse_map(ratios)
        logging.info("reverse map built in %d ms.",
                     util.current_millis() - start_time)

    def run_logs(self):
        return [self.update_log, self.motif_log]

    def __build_reverse_map(self, ratios):
        """build a map that reconstructs the original row name from
        a feature id"""
        def feature_id_for(gene):
            """convenience method to return the feature id for a gene"""
            feature_ids = self.organism.feature_ids_for([gene])
            if len(feature_ids) > 0:
                return feature_ids[0]
            else:
                return None

        result = {}
        num_not_found = 0
        for row_name in ratios.row_names:
            feature_id = feature_id_for(row_name)
            if feature_id is not None:
                result[feature_id] = row_name
            else:
                num_not_found += 1
        if num_not_found > 0:
            logging.warn("%d genes not found in synonyms.", num_not_found)
        return result

    def name(self):
        """returns the name of this scoring function"""
        return "Motif"""

    def compute(self, iteration_result, ref_matrix=None):
        """override base class compute() method, behavior is more complicated,
        since it nests Motif and MEME runs"""
        return self.__compute(iteration_result, False, ref_matrix)

    def compute_force(self, iteration_result, ref_matrix=None):
        """override base class compute() method, behavior is more complicated,
        since it nests Motif and MEME runs"""
        result = self.__compute(iteration_result, True, ref_matrix)
        # and pickle the last result for diagnostics
        with open(self.pickle_path(), 'w') as outfile:
            cPickle.dump(result, outfile)
        return result

    def last_cached(self):
        return self.last_result

    def matrix_pickle_path(self):
        return "%s/%s_matrix_last.pkl" % (self.config_params['output_dir'],
                                          self.name())

    def __compute(self, iteration_result, force, ref_matrix=None):
        """compute method for the specified iteration
        Note: will return None if not computed yet and the result of a previous
        scoring if the function is not supposed to actually run in this iteration
        """
        iteration = iteration_result['iteration']
        if force or self.motif_in_iteration(iteration):  # meme.iter in R
            logging.info("Running Motifing for sequence type '%s'...", self.seqtype)
            # running MEME and store the result for the non-motifing iterations
            # to reuse
            # Note: currently, iteration results are only computed here
            num_motifs = int(self.num_motif_func(iteration))
            self.__last_iteration_result = {'iteration': iteration}
            self.all_pvalues = self.compute_pvalues(self.__last_iteration_result,
                                                    num_motifs)

        if self.all_pvalues is not None and (force or self.update_in_iteration(iteration)):  # mot.iter in R
            logging.info("UPDATING MOTIF SCORES in iteration %d with scaling: %f",
                         iteration, self.scaling(iteration))
            self.last_result = pvalues2matrix(self.all_pvalues, self.num_clusters(),
                                              self.gene_names(), self.reverse_map)

        self.update_log.log(iteration, self.update_in_iteration(iteration),
                            self.scaling(iteration))
        self.motif_log.log(iteration, self.motif_in_iteration(iteration),
                           self.scaling(iteration))

        if 'motifs' not in iteration_result:
            iteration_result['motifs'] = {}

        # remove the 'iteration' key that we inserted
        if "iteration" in self.__last_iteration_result:
            del self.__last_iteration_result['iteration']

        iteration_result['motifs'][self.seqtype] = self.__last_iteration_result

        # for stats, support multiple sequence types
        if 'motif-pvalue' not in iteration_result:
            iteration_result['motif-pvalue'] = {}

        iteration_result['motif-pvalue'][self.seqtype] = compute_mean_score(
            self.last_result, self.membership, self.organism)
        return self.last_result

    def compute_pvalues(self, iteration_result, num_motifs):
        """Compute motif scores.
        The result is a dictionary from cluster -> (feature_id, pvalue)
        containing a sparse gene-to-pvalue mapping for each cluster

        In order to influence the sequences
        that go into meme, the user can specify a list of sequence filter
        functions that have the signature
        (seqs, feature_ids, distance) -> seqs
        These filters are applied in the order they appear in the list.
        """
        global MOTIF_PARAMS, SEQUENCE_FILTERS, ORGANISM, MEMBERSHIP

        cluster_pvalues = {}
        min_cluster_rows_allowed = self.config_params['memb.min_cluster_rows_allowed']
        max_cluster_rows_allowed = self.config_params['memb.max_cluster_rows_allowed']
        use_multiprocessing = self.config_params[
            scoring.KEY_MULTIPROCESSING]

        # extract the sequences for each cluster, slow
        start_time = util.current_millis()
        SEQUENCE_FILTERS = self.__sequence_filters
        ORGANISM = self.organism
        MEMBERSHIP = self.membership

        pool = mp.Pool()
        cluster_seqs_params = [(cluster, self.seqtype)
                               for cluster in xrange(1, self.num_clusters() + 1)]
        seqs_list = pool.map(cluster_seqs, cluster_seqs_params)
        SEQUENCE_FILTERS = None
        ORGANISM = None
        MEMBERSHIP = None
        logging.info("prepared sequences in %d ms.", util.current_millis() - start_time)

        # Make the parameters, this is fast enough
        start_time = util.current_millis()
        params = []
        for cluster in xrange(1, self.num_clusters() + 1):
            # Pass the previous run's seed if possible
            if self.__last_motif_infos is not None:
                previous_motif_infos = self.__last_motif_infos.get(cluster, None)
            else:
                previous_motif_infos = None

            seqs, feature_ids = seqs_list[cluster - 1]
            params.append(ComputeScoreParams(iteration_result['iteration'], cluster,
                                             feature_ids,
                                             seqs,
                                             self.used_seqs,
                                             self.meme_runner(),
                                             min_cluster_rows_allowed,
                                             max_cluster_rows_allowed,
                                             num_motifs,
                                             previous_motif_infos,
                                             self.config_params.get('keep_memeout', False),
                                             self.config_params['output_dir'],
                                             self.config_params['num_iterations'],
                                             self.config_params['debug']))

        logging.info("prepared MEME parameters in %d ms.",
                     util.current_millis() - start_time)

        # create motif result map if necessary
        for cluster in xrange(1, self.num_clusters() + 1):
            if not cluster in iteration_result:
                iteration_result[cluster] = {}

        # compute and store motif results
        MOTIF_PARAMS = params
        self.__last_motif_infos = {}
        if use_multiprocessing:
            pool = mp.Pool()
            results = pool.map(compute_cluster_score, xrange(1, self.num_clusters() + 1))

            for cluster in xrange(1, self.num_clusters() + 1):
                pvalues, run_result = results[cluster - 1]
                cluster_pvalues[cluster] = pvalues
                if run_result:
                    self.__last_motif_infos[cluster] = run_result.motif_infos
                iteration_result[cluster]['motif-info'] = meme_json(run_result)
                iteration_result[cluster]['pvalues'] = pvalues
            pool.close()
            pool.join()
        else:
            for cluster in xrange(1, self.num_clusters() + 1):
                pvalues, run_result = compute_cluster_score(cluster)
                cluster_pvalues[cluster] = pvalues
                if run_result:
                    self.__last_motif_infos[cluster] = run_result.motif_infos
                iteration_result[cluster]['motif-info'] = meme_json(run_result)
                iteration_result[cluster]['pvalues'] = pvalues

        # cleanup
        MOTIF_PARAMS = None
        return cluster_pvalues


def cluster_seqs(params):
    """Retrieves the sequences for a cluster. Designed to run in in pool.map()"""
    global SEQUENCE_FILTERS, ORGANISM, MEMBERSHIP
    cluster, seqtype = params
    genes = sorted(MEMBERSHIP.rows_for_cluster(cluster))
    feature_ids = ORGANISM.feature_ids_for(genes)
    seqs = ORGANISM.sequences_for_genes_search(
        feature_ids, seqtype=seqtype)
    for sequence_filter in SEQUENCE_FILTERS:
        seqs = sequence_filter(seqs, feature_ids)
    if len(seqs) == 0:
        logging.warn('Cluster %i with %i genes: no sequences!',
                     cluster, len(seqs))
    return (seqs, feature_ids)


def meme_json(run_result):
    result = []
    if run_result is not None:
        motif_annotations = {}  # map motif_num -> [annotations]
        for gene in run_result.annotations:
            for annotation in run_result.annotations[gene]:
                # motif_num is either positive or negative, indicating forward/reverse
                motif_num = annotation[2]
                key = abs(motif_num)
                reverse = motif_num < 0
                if key not in motif_annotations:
                    motif_annotations[key] = []

                motif_annotations[key].append(
                    {'gene': gene,
                     'position': annotation[1],
                     'pvalue': annotation[0],
                     'reverse': reverse})
        for motif_info in run_result.motif_infos:
            motif_num = motif_info.motif_num
            motif_annot = []
            if motif_num in motif_annotations:
                motif_annot = motif_annotations[motif_num]
            result.append({'motif_num': motif_num,
                           'pssm': motif_info.pssm,
                           'evalue': motif_info.evalue,
                           'annotations': motif_annot,
                           'sites': motif_info.sites})
    return result


def compute_cluster_score(cluster):
    """This function computes the MEME score for a cluster"""
    global MOTIF_PARAMS
    params = MOTIF_PARAMS[cluster - 1]
    pvalues = {}
    run_result = None
    nseqs = len(params.seqs)
    logging.info('Cluster %d, # sequences: %d', params.cluster, nseqs)
    if (nseqs >= params.min_cluster_rows and nseqs <= params.max_cluster_rows):
        run_result = params.meme_runner(params)
        pe_values = run_result.pe_values
        for feature_id, pvalue, evalue in pe_values:
            pvalues[feature_id] = pvalue
    else:
        logging.info("# seqs (= %d) outside of defined limits, "
                     "skipping cluster %d", len(params.seqs), params.cluster)
    return pvalues, run_result


class MemeScoringFunction(MotifScoringFunctionBase):
    """Scoring function for motifs"""

    def __init__(self, organism, membership, ratios,
                 meme_suite,
                 seqtype='upstream',
                 sequence_filters=[],
                 scaling_func=None,
                 num_motif_func=None,
                 update_in_iteration=None,
                 motif_in_iteration=None,
                 config_params=None):
        """creates a ScoringFunction"""
        MotifScoringFunctionBase.__init__(self, organism, membership,
                                          ratios, meme_suite, seqtype,
                                          sequence_filters,
                                          scaling_func,
                                          num_motif_func,
                                          update_in_iteration,
                                          motif_in_iteration,
                                          config_params)

    def name(self):
        """returns the name of this scoring function"""
        return "MEME"

    def meme_runner(self):
        """returns the MEME runner object"""
        return self.meme_suite


class WeederScoringFunction(MotifScoringFunctionBase):
    """Motif scoring function that runs Weeder instead of MEME"""

    def __init__(self, organism, membership, ratios,
                 meme_suite, seqtype,
                 sequence_filters=[],
                 scaling_func=None,
                 num_motif_func=None,
                 update_in_iteration=None,
                 motif_in_iteration=None,
                 config_params=None):
        """creates a scoring function"""
        MotifScoringFunctionBase.__init__(self, organism, membership, ratios,
                                          meme_suite, seqtype,
                                          sequence_filters,
                                          scaling_func,
                                          num_motif_func,
                                          update_in_iteration,
                                          motif_in_iteration,
                                          config_params)

    def name(self):
        """returns the name of this scoring function"""
        return "Weeder"

    def meme_runner(self):
        """returns the MEME runner object"""
        return WeederRunner(self.meme_suite)


class WeederRunner:
    """Wrapper around Weeder so we can use the multiprocessing module.
    The function basically runs Weeder ont the specified set of sequences,
    converts its output to a MEME output file and runs MAST on the MEME output
    to generate a MEME run result.
    """

    def __init__(self, meme_suite, background_file=None, remove_tempfiles=True):
        """create a runner object"""
        self.meme_suite = meme_suite
        self.__background_file = background_file
        self.__remove_tempfiles = remove_tempfiles

    def __call__(self, params):
        """call the runner like a function"""
        with tempfile.NamedTemporaryFile(prefix='weeder.fasta',
                                         delete=False) as outfile:
            filename = outfile.name
            logging.info("Run Weeder on FASTA file: '%s'", filename)
            st.write_sequences_to_fasta_file(outfile, params.seqs.items())

        pssms = weeder.run_weeder(filename)
        meme_outfile = '%s.meme' % filename
        dbfile = self.meme_suite.make_sequence_file(
            [(feature_id, locseq[1])
             for feature_id, locseq in params.used_seqs.items()])
        logging.info("# PSSMS created: %d %s", len(pssms), str([i.consensus_motif() for i in pssms]))
        logging.info("run MAST on '%s'", meme_outfile)

        motif_infos = []
        for i in xrange(len(pssms)):
            pssm = pssms[i]
            motif_infos.append(meme.MemeMotifInfo(pssm.values, i + 1,
                                                  pssm.sequence_length(),
                                                  len(pssm.sites),
                                                  None, pssm.e_value,
                                                  pssm.sites))

        try:
            mast_out = self.meme_suite.mast(
                meme_outfile, dbfile,
                self.meme_suite.global_background_file())
            pe_values, annotations = meme.read_mast_output(mast_out,
                                                           params.seqs.keys())
            return meme.MemeRunResult(pe_values, annotations, motif_infos)
        except:
            return meme.MemeRunResult([], {}, [])
        finally:
            if self.__remove_tempfiles:
                for fileExtension in ['', '.wee', '.mix', '.html', '.meme', '.1.f1', '.1.f2', '.2.f1', '.2.f2']:
                    tmpName = filename+fileExtension
                    if os.path.exists(tmpName):
                        try:
                            os.remove(tmpName)
                        except:
                            logging.warn("could not remove tmp file:'%s'", tmpName)
            try:
                os.remove(dbfile)
            except:
                logging.warn("could not remove tmp file:'%s'", dbfile)
            #if self.__background_file==None:
            #    try:
            #        os.remove(bgFile)
            #    except:
            #        logging.warn("could not remove tmp file: '%s'", bgfile)
