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
            startPos = - distance[0]
            chars[ startPos : startPos + 4] = "NNNN"
            seqs[feature_id] = "".join(chars)
        return seqs
    return remove_atgs_filter


class MinPValueFilter:
    """A minimum p-value filter. Implemented as a class, so multiprocessing
    can deal with it"""
    def __init__(self, min_value):
        self.__min_value = min_value

    def __call__(self, values):
        """all values that are below min_value are set to the
        minimum value above min_value"""
        allowed_vals = [value for key, value in values.items()
                        if value > self.__min_value]
        result = {}
        if len(allowed_vals) == 0:
            logging.warn("all values are below the threshold -> no result !!")
        else:
            min_allowed = min(allowed_vals)
            for key, value in values.items():
                if value < self.__min_value:
                    result[key] = min_allowed
                else:
                    result[key] = value
        return result


class MotifScoringFunctionBase(scoring.ScoringFunctionBase):
    """Base class for motif scoring functions that use MEME"""

    def __init__(self, organism, membership, matrix,
                 meme_suite, seqtype,
                 sequence_filters=[],
                 pvalue_filter=None,
                 scaling_func=None,
                 run_in_iteration=lambda iteration: True,
                 config_params=None):
        """creates a ScoringFunction"""
        scoring.ScoringFunctionBase.__init__(self, membership,
                                             matrix, scaling_func,
                                             config_params)
        # attributes accessible by subclasses
        self.organism = organism
        self.meme_suite = meme_suite
        self.seqtype = seqtype
        self.run_in_iteration = run_in_iteration
        self.__sequence_filters = sequence_filters
        self.__pvalue_filter = pvalue_filter

        # precompute the sequences for all genes that are referenced in the
        # input ratios, they are used as a basis to compute the background
        # distribution for every cluster
        self.used_seqs = organism.sequences_for_genes_scan(
            sorted(matrix.row_names()), seqtype=self.seqtype)
        logging.info("used sequences for motifing retrieved")
        logging.info("building reverse map...")
        start_time = util.current_millis()
        self.reverse_map = self.__build_reverse_map(matrix)
        logging.info("reverse map built in %d ms.",
                     util.current_millis() - start_time)

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
        for row_name in matrix.row_names():
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
        """compute method, iteration is the 0-based iteration number"""
        iteration = iteration_result['iteration']
        if self.run_in_iteration(iteration):
            global_start_time = util.current_millis()
            pvalues = self.compute_pvalues(iteration_result)
            remapped = {}
            for cluster in pvalues:
                pvalues_k = pvalues[cluster]
                pvalues_genes = {}
                for feature_id, pvalue in pvalues_k.items():
                    pvalues_genes[self.reverse_map[feature_id]] = pvalue
                remapped[cluster] = pvalues_genes

            # convert remapped to an actual scoring matrix
            matrix = dm.DataMatrix(len(self.gene_names()), self.num_clusters(),
                                   self.gene_names())
            for row_index in xrange(matrix.num_rows()):
                row = matrix.row_name(row_index)
                for cluster in xrange(1, self.num_clusters() + 1):
                    if (cluster in remapped.keys() and
                        row in remapped[cluster].keys()):
                        matrix[row_index][cluster - 1] = remapped[cluster][row]
            global_elapsed = util.current_millis() - global_start_time
            logging.info("GLOBAL MOTIF TIME: %d seconds",
                         (global_elapsed / 1000.0))
            return matrix
        else:
            return None

    def compute_pvalues(self, iteration_result):
        """Compute motif scores. In order to influence the sequences
        that go into meme, the user can specify a list of sequence filter
        functions that have the signature
        (seqs, feature_ids, distance) -> seqs
        These filters are applied in the order they appear in the list
        """
        def apply_sequence_filters(seqs, feature_ids):
            """apply all filters in the filters list in order"""
            for sequence_filter in self.__sequence_filters:
                seqs = sequence_filter(seqs, feature_ids)
            return seqs

        cluster_pvalues = {}
        min_cluster_rows_allowed = self.config_params[
            scoring.KEY_MOTIF_MIN_CLUSTER_ROWS_ALLOWED]
        max_cluster_rows_allowed = self.config_params[
            scoring.KEY_MOTIF_MAX_CLUSTER_ROWS_ALLOWED]
        use_multiprocessing = self.config_params[
            scoring.KEY_MULTIPROCESSING]

        # create parameters
        params = []
        for cluster in xrange(1, self.num_clusters() + 1):
            genes = sorted(self.rows_for_cluster(cluster))
            feature_ids = self.organism.feature_ids_for(genes)
            seqs = self.organism.sequences_for_genes_search(
                genes, seqtype=self.seqtype)
            # should I stay or should I go? Frank Schmitz
            seqs = apply_sequence_filters(seqs, feature_ids)

            params.append(ComputeScoreParams(cluster, genes, seqs,
                                             self.used_seqs,
                                             self.meme_runner(),
                                             self.__pvalue_filter,
                                             min_cluster_rows_allowed,
                                             max_cluster_rows_allowed))

        # create motif result map if necessary
        if not 'motifs' in iteration_result:
            iteration_result['motifs'] = {}
        for cluster in xrange(1, self.num_clusters() + 1):
            if not cluster in iteration_result['motifs']:
                iteration_result['motifs'][cluster] = { }
            if not self.seqtype in iteration_result['motifs'][cluster]:
                iteration_result['motifs'][cluster][self.seqtype] = { }

        # compute and store motif results
        if use_multiprocessing:
            pool = mp.Pool()
            results = pool.map(compute_cluster_score, params)

            for cluster in xrange(1, self.num_clusters() + 1):
                pvalues, run_result = results[cluster - 1]
                cluster_pvalues[cluster] = pvalues
                iteration_result['motifs'][cluster][self.seqtype] = meme_json(run_result)
            pool.close()
        else:
            for cluster in xrange(1, self.num_clusters() + 1):
                pvalues, run_result = compute_cluster_score(params[cluster - 1])
                cluster_pvalues[cluster] = pvalues
                iteration_result['motifs'][cluster][self.seqtype] = meme_json(run_result)

        return cluster_pvalues


def meme_json(run_result):
    result = []
    if run_result != None:
        motif_annotations = {}  # map motif_num -> [annotations]
        for gene in run_result.annotations:
            for annotation in run_result.annotations[gene]:
                motif_num = annotation[2]
                if motif_num not in motif_annotations:
                    motif_annotations[motif_num] = []

                motif_annotations[motif_num].append(
                    {'gene': gene, 'position': annotation[0], 'pvalue': annotation[1]})

        for motif_info in run_result.motif_infos:
            motif_num = motif_info.motif_num()
            motif_annot = None
            if motif_num in motif_annotations:
                motif_annot = motif_annotations[motif_num]
            result.append({'motif_num': motif_num,
                           'pssm': motif_info.pssm(),
                           'evalue': motif_info.evalue(),
                           'annotations': motif_annot})
    return result


class ComputeScoreParams:
    """representation of parameters to compute motif scores"""
    def __init__(self, cluster, genes, seqs, used_seqs, meme_runner,
                 pvalue_filter, min_cluster_rows, max_cluster_rows):
        """constructor"""
        self.cluster = cluster
        self.genes = genes
        self.seqs = seqs
        self.used_seqs = used_seqs
        self.meme_runner = meme_runner
        self.pvalue_filter = pvalue_filter
        self.min_cluster_rows = min_cluster_rows
        self.max_cluster_rows = max_cluster_rows


def compute_cluster_score(params):
    """This function computes the MEME score for a cluster"""
    pvalues = {}
    run_result = None
    nseqs = len(params.seqs)
    logging.info('computing cluster score for %s seqs...' %nseqs)
    if (nseqs >= params.min_cluster_rows
        and nseqs <= params.max_cluster_rows):
        run_result = params.meme_runner(params.seqs, params.used_seqs)
        pe_values = run_result.pe_values
        for feature_id, pvalue, evalue in pe_values:
            pvalues[feature_id] = np.log(pvalue)
        if params.pvalue_filter != None:
            pvalues = params.pvalue_filter(pvalues)
    else:
        logging.info("# seqs (= %d) outside of defined limits, "
                     "skipping cluster %d", len(params.seqs), params.cluster)
    return pvalues, run_result


class MemeScoringFunction(MotifScoringFunctionBase):
    """Scoring function for motifs"""

    def __init__(self, organism, membership, matrix,
                 meme_suite,
                 seqtype='Promoter',
                 sequence_filters=[],
                 pvalue_filter=None,
                 scaling_func=None,
                 run_in_iteration=scoring.default_motif_iterations,
                 config_params=None):
        """creates a ScoringFunction"""
        MotifScoringFunctionBase.__init__(self, organism, membership,
                                          matrix, meme_suite, seqtype,
                                          sequence_filters, pvalue_filter,
                                          scaling_func, run_in_iteration,
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
                 sequence_filters=[], pvalue_filter=None,
                 scaling_func=None,
                 run_in_iteration=scoring.default_motif_iterations,
                 config_params=None):
        """creates a scoring function"""
        MotifScoringFunctionBase.__init__(self, organism, membership, matrix,
                                          meme_suite, seqtype,
                                          sequence_filters, pvalue_filter,
                                          scaling_func, run_in_iteration, config_params)

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

    def __call__(self, seqs, all_seqs):
        """call the runner like a function"""
        with tempfile.NamedTemporaryFile(prefix='weeder.fasta',
                                         delete=False) as outfile:
            filename = outfile.name
            logging.info("Run Weeder on FASTA file: '%s'", filename)
            st.write_sequences_to_fasta_file(outfile, seqs.items())

        pssms = weeder.run_weeder(filename)
        meme_outfile = '%s.meme' % filename
        dbfile = self.meme_suite.make_sequence_file(
            [(feature_id, locseq[1])
             for feature_id, locseq in all_seqs.items()])
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
                                                           seqs.keys())
            return meme.MemeRunResult(pe_values, annotations, motif_infos)
        except:
            return meme.MemeRunResult([], {}, [])


