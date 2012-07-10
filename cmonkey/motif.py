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


def default_nmotif_fun(iteration, num_iterations):
    """default function to compute the nmotif parameter for MEME dependent on
    the iteration"""
    if iteration <= (num_iterations * 2 / 3):
        return 1
    else:
        return 2


# Applicable sequence filters
def unique_filter(seqs, feature_ids):
    """returns a map that contains only the keys that are in
    feature_ids and only contains unique sequences"""
    unique_seqs = {}
    for feature_id in feature_ids:
        if (feature_id in seqs
            and seqs[feature_id] not in unique_seqs.values()):
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
    if pvalue_matrix == None:
        return 0.0
    values = []
    pvalues = pvalue_matrix.values
    for cluster in xrange(1, membership.num_clusters() + 1):
        cluster_rows = membership.rows_for_cluster(cluster)
        row_indexes = pvalue_matrix.row_indexes(cluster_rows)
        for row in row_indexes:
            values.append(pvalues[row][cluster - 1])
    return np.median(values)

# Readonly structure to avoid passing it to the forked child processes
MOTIF_PARAMS = None

class MotifScoringFunctionBase(scoring.ScoringFunctionBase):
    """Base class for motif scoring functions that use MEME"""

    def __init__(self, organism, membership, matrix,
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
        scoring.ScoringFunctionBase.__init__(self, membership,
                                             matrix, scaling_func,
                                             run_in_iteration=None,
                                             config_params=config_params)
        # attributes accessible by subclasses
        self.organism = organism
        self.meme_suite = meme_suite
        self.seqtype = seqtype
        self.update_in_iteration = update_in_iteration
        self.motif_in_iteration = motif_in_iteration
        self.num_motif_func = num_motif_func

        self.__sequence_filters = sequence_filters
        self.__last_computed_result = None
        self.__last_run_results = None
        self.__last_pvalues = None
        self.__last_iteration_result = {}

        self.update_log = scoring.RunLog("motif-score-" + seqtype)
        self.motif_log = scoring.RunLog("motif-motif-" + seqtype)

        # precompute the sequences for all genes that are referenced in the
        # input ratios, they are used as a basis to compute the background
        # distribution for every cluster
        self.used_seqs = organism.sequences_for_genes_scan(
            sorted(matrix.row_names), seqtype=self.seqtype)
        logging.info("used sequences for motifing retrieved")
        logging.info("building reverse map...")
        start_time = util.current_millis()
        self.reverse_map = self.__build_reverse_map(matrix)
        logging.info("reverse map built in %d ms.",
                     util.current_millis() - start_time)

    def run_logs(self):
        return [self.update_log, self.motif_log]

    def __build_reverse_map(self, matrix):
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
        for row_name in matrix.row_names:
            feature_id = feature_id_for(row_name)
            if feature_id != None:
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
        return self.__compute(iteration_result, True, ref_matrix)

    def __compute(self, iteration_result, force, ref_matrix=None):
        """compute method for the specified iteration
        Note: will return None if not computed yet and the result of a previous
        scoring if the function is not supposed to actually run in this iteration
        """
        iteration = iteration_result['iteration']

        if force or self.motif_in_iteration(iteration):  # meme.iter in R
            logging.info('Running Motifing...')
            # running MEME and store the result for the non-motifing iterations
            # to reuse
            # Note: currently, iteration results are only computed here
            num_motifs = self.num_motif_func(iteration,
                                             self.config_params['num_iterations'])
            self.__last_iteration_result = {}
            self.__last_pvalues = self.compute_pvalues(self.__last_iteration_result,
                                                       num_motifs)
            #for cluster in sorted(self.__last_pvalues.keys()):
            #    print "CLUSTER ", cluster
            #    for gene in sorted(self.__last_pvalues[cluster].keys()):
            #        print "%s -> %f" % (gene, self.__last_pvalues[cluster][gene])

        if self.__last_pvalues != None and (
            force or self.update_in_iteration(iteration)):  # mot.iter in R
            logging.info('Recomputing motif scores...')
            start_time = util.current_millis()
            # running the scoring itself
            remapped = {}
            for cluster in self.__last_pvalues:
                pvalues_k = self.__last_pvalues[cluster]
                pvalues_genes = {}
                for feature_id, pvalue in pvalues_k.items():
                    pvalues_genes[self.reverse_map[feature_id]] = pvalue
                remapped[cluster] = pvalues_genes

            current = util.current_millis()
            logging.info("built remapped in %f s.", (current - start_time) / 1000.0)
            start_time = current

            # convert remapped to an actual scoring matrix
            matrix = dm.DataMatrix(len(self.gene_names()), self.num_clusters(),
                                   self.gene_names())
            mvalues = matrix.values
            for row_index in xrange(matrix.num_rows()):
                row = matrix.row_names[row_index]
                for cluster in xrange(1, self.num_clusters() + 1):
                    if (cluster in remapped.keys() and
                        row in remapped[cluster].keys()):
                        mvalues[row_index][cluster - 1] = remapped[cluster][row]
            matrix.apply_log()

            current = util.current_millis()
            logging.info("built scoring matrix in %f s.", (current - start_time) / 1000.0)
            start_time = current

            matrix.fix_extreme_values()
            current = util.current_millis()
            logging.info("fixed extreme values in %f s.", (current - start_time) / 1000.0)

            self.__last_computed_result = matrix

        self.update_log.log(self.update_in_iteration(iteration),
                            self.scaling(iteration))
        self.motif_log.log(self.motif_in_iteration(iteration),
                           self.scaling(iteration))

        if 'motifs' not in iteration_result:
            iteration_result['motifs'] = {}
        iteration_result['motifs'][self.seqtype] = self.__last_iteration_result

        # for stats, support multiple sequence types
        if 'motif-pvalue' not in iteration_result:
            iteration_result['motif-pvalue'] = {}
        iteration_result['motif-pvalue'][self.seqtype] = compute_mean_score(
            self.__last_computed_result,
            self.membership(),
            self.organism)

        return self.__last_computed_result

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
        global MOTIF_PARAMS
        def apply_sequence_filters(seqs, feature_ids):
            """apply all filters in the filters list in order"""
            for sequence_filter in self.__sequence_filters:
                seqs = sequence_filter(seqs, feature_ids)
            return seqs

        cluster_pvalues = {}
        min_cluster_rows_allowed = self.config_params['memb.min_cluster_rows_allowed']
        max_cluster_rows_allowed = self.config_params['memb.max_cluster_rows_allowed']
        use_multiprocessing = self.config_params[
            scoring.KEY_MULTIPROCESSING]

        # create parameters
        params = []
        for cluster in xrange(1, self.num_clusters() + 1):
            genes = sorted(self.rows_for_cluster(cluster))
            feature_ids = self.organism.feature_ids_for(genes)
            seqs = self.organism.sequences_for_genes_search(
                genes, seqtype=self.seqtype)
            seqs = apply_sequence_filters(seqs, feature_ids)
            if len(seqs) == 0:
                logging.warn('Cluster %i with %i genes: no sequences!' \
                                %(cluster,len(seqs)) )

            # Pass the previous run's seed if possible
            if (self.__last_run_results != None and
                cluster in self.__last_run_results.keys() and
                self.__last_run_results[cluster] != None and
                self.__last_run_results[cluster].motif_infos != None):
                previous_motif_infos = self.__last_run_results[cluster].motif_infos
            else:
                previous_motif_infos = None

            params.append(ComputeScoreParams(cluster,
                                             feature_ids,
                                             seqs,
                                             self.used_seqs,
                                             self.meme_runner(),
                                             min_cluster_rows_allowed,
                                             max_cluster_rows_allowed,
                                             num_motifs,
                                             previous_motif_infos))

        # create motif result map if necessary
        for cluster in xrange(1, self.num_clusters() + 1):
            if not cluster in iteration_result:
                iteration_result[cluster] = { }

        # compute and store motif results
        MOTIF_PARAMS = params
        self.__last_run_results = {}
        if use_multiprocessing:
            pool = mp.Pool()
            results = pool.map(compute_cluster_score, xrange(1, self.num_clusters() + 1))

            for cluster in xrange(1, self.num_clusters() + 1):
                pvalues, run_result = results[cluster - 1]
                cluster_pvalues[cluster] = pvalues
                self.__last_run_results[cluster] = run_result
                iteration_result[cluster]['motif-info'] = meme_json(run_result)
                iteration_result[cluster]['pvalues'] = pvalues
            pool.close()
        else:
            for cluster in xrange(1, self.num_clusters() + 1):
                pvalues, run_result = compute_cluster_score(cluster)
                cluster_pvalues[cluster] = pvalues
                self.__last_run_results[cluster] = run_result
                iteration_result[cluster]['motif-info'] = meme_json(run_result)
                iteration_result[cluster]['pvalues'] = pvalues

        # cleanup
        MOTIF_PARAMS = None
        return cluster_pvalues


def meme_json(run_result):
    result = []
    if run_result != None:
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
                           'annotations': motif_annot})
    return result


class ComputeScoreParams:
    """representation of parameters to compute motif scores"""
    def __init__(self, cluster, feature_ids, seqs, used_seqs, meme_runner,
                 min_cluster_rows, max_cluster_rows,
                 num_motifs, previous_motif_infos=None):
        """constructor"""
        self.cluster = cluster
        self.feature_ids = feature_ids  # to preserve the order
        self.seqs = seqs                # sequences in the cluster
        self.used_seqs = used_seqs      # all available sequences in the organism
        self.meme_runner = meme_runner
        self.min_cluster_rows = min_cluster_rows
        self.max_cluster_rows = max_cluster_rows
        self.num_motifs = num_motifs
        self.previous_motif_infos = previous_motif_infos


def compute_cluster_score(cluster):
    """This function computes the MEME score for a cluster"""
    global MOTIF_PARAMS
    params = MOTIF_PARAMS[cluster - 1]
    pvalues = {}
    run_result = None
    nseqs = len(params.seqs)
    logging.info('%d: computing cluster score for %s seqs...', params.cluster, nseqs)
    if (nseqs >= params.min_cluster_rows
        and nseqs <= params.max_cluster_rows):
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

    def __init__(self, organism, membership, matrix,
                 meme_suite,
                 seqtype='upstream',
                 sequence_filters=[],
                 scaling_func=None,
                 num_motif_func=None,
                 update_in_iteration=scoring.default_motif_iterations,
                 motif_in_iteration=scoring.default_meme_iterations,
                 config_params=None):
        """creates a ScoringFunction"""
        MotifScoringFunctionBase.__init__(self, organism, membership,
                                          matrix, meme_suite, seqtype,
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

    def __init__(self, organism, membership, matrix,
                 meme_suite, seqtype,
                 sequence_filters=[],
                 scaling_func=None,
                 num_motif_func=None,
                 update_in_iteration=None,
                 motif_in_iteration=None,
                 config_params=None):
        """creates a scoring function"""
        MotifScoringFunctionBase.__init__(self, organism, membership, matrix,
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

    def __init__(self, meme_suite):
        """create a runner object"""
        self.meme_suite = meme_suite

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
        logging.info("# PSSMS created: %d", len(pssms))
        logging.info("run MAST on '%s'", meme_outfile)

        motif_infos = []
        for i in xrange(len(pssms)):
            pssm = pssms[i]
            motif_infos.append(meme.MemeMotifInfo(pssm.values(), i + 1,
                                                  pssm.sequence_length(),
                                                  len(pssm.sites()),
                                                  None, pssm.evalue(),
                                                  pssm.sites()))

        try:
            mast_out = self.meme_suite.mast(
                meme_outfile, dbfile,
                self.meme_suite.global_background_file())
            pe_values, annotations = meme.read_mast_output(mast_out,
                                                           params.seqs.keys())
            return meme.MemeRunResult(pe_values, annotations, motif_infos)
        except:
            return meme.MemeRunResult([], {}, [])
