"""motif.py - cMonkey motif related processing
This module captures the motif-specific scoring component
of cMonkey.

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import logging
import numpy
import scoring
import datamatrix as dm
import weeder
import meme
import tempfile
import seqtools as st


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


def make_min_value_filter(min_value):
    """A function generator which creates a value filter"""

    def min_value_filter(values):
        """all values that are below min_value are set to the
        minimum value above min_value"""
        allowed_vals = [value for key, value in values.items()
                        if value > min_value]
        result = {}
        if len(allowed_vals) == 0:
            logging.warn("all values are below the threshold -> no result !!")
        else:
            min_allowed = min(allowed_vals)
            for key, value in values.items():
                if value < min_value:
                    result[key] = min_allowed
                else:
                    result[key] = value
        return result
    return min_value_filter


class MotifScoringFunctionBase(scoring.ScoringFunctionBase):
    """Base class for motif scoring functions that use MEME"""

    def __init__(self, organism, membership, matrix,
                 meme_suite, seqtype,
                 sequence_filters=[],
                 pvalue_filter=None,
                 weight_func=None,
                 interval=0,
                 config_params=None):
        """creates a ScoringFunction"""
        scoring.ScoringFunctionBase.__init__(self, membership,
                                             matrix, weight_func,
                                             config_params)
        # attributes accessible by subclasses
        self.organism = organism
        self.meme_suite = meme_suite
        self.seqtype = seqtype
        self.interval = interval
        self.__sequence_filters = sequence_filters
        self.__pvalue_filter = pvalue_filter

        # precompute the sequences for all genes that are referenced in the
        # input ratios, they are used as a basis to compute the background
        # distribution for every cluster
        self.used_seqs = organism.sequences_for_genes_scan(
            sorted(matrix.row_names()), seqtype=self.seqtype)
        logging.info("used sequences for motifing retrieved")
        logging.info("building reverse map...")
        self.reverse_map = self.__build_reverse_map(matrix)
        logging.info("reverse map built.")

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
        for row_name in matrix.row_names():
            feature_id = feature_id_for(row_name)
            if feature_id != None:
                result[feature_id] = row_name
        return result

    def name(self):
        """returns the name of this scoring function"""
        return "Motif"""

    def compute(self, iteration, ref_matrix=None):
        """compute method, iteration is the 0-based iteration number"""
        if (self.interval == 0 or
            (iteration > 0 and (iteration % self.interval == 0))):
            print "RUN MOTIF SCORING IN ITERATION ", iteration
            pvalues = self.compute_pvalues(iteration)
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
            for row in matrix.row_names():
                row_index = matrix.row_names().index(row)
                for cluster in range(1, self.num_clusters() + 1):
                    if (cluster in remapped.keys() and
                        row in remapped[cluster].keys()):
                        matrix[row_index][cluster - 1] = remapped[cluster][row]
            return matrix
        else:
            return None

    def compute_pvalues(self, iteration):
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
        meme_run_results = {}
        min_cluster_rows_allowed = self.config_params[
            scoring.KEY_MOTIF_MIN_CLUSTER_ROWS_ALLOWED]
        max_cluster_rows_allowed = self.config_params[
            scoring.KEY_MOTIF_MAX_CLUSTER_ROWS_ALLOWED]

        for cluster in range(1, self.num_clusters() + 1):
            logging.info("compute motif scores for cluster %d", cluster)
            genes = sorted(self.rows_for_cluster(cluster))
            feature_ids = self.organism.feature_ids_for(genes)
            seqs = self.organism.sequences_for_genes_search(
                genes, seqtype=self.seqtype)
            seqs = apply_sequence_filters(seqs, feature_ids)
            if (len(seqs) >= min_cluster_rows_allowed
                and len(seqs) <= max_cluster_rows_allowed):
                meme_run_results[cluster] = self.compute_meme_run(seqs)

                pvalues = {}
                pe_values = meme_run_results[cluster].pe_values
                for feature_id, pvalue, evalue in pe_values:
                    pvalues[feature_id] = numpy.log(pvalue)
                if self.__pvalue_filter != None:
                    pvalues = self.__pvalue_filter(pvalues)
                cluster_pvalues[cluster] = pvalues

                for motif_info in meme_run_results[cluster].motif_infos:
                    print ("consensus: ", motif_info.consensus_string(),
                           " evalue: ", motif_info.evalue())
            else:
                logging.info("# seqs (= %d) outside of defined limits, "
                             "skipping cluster %d", len(seqs), cluster)

        return cluster_pvalues


class MemeScoringFunction(MotifScoringFunctionBase):
    """Scoring function for motifs"""

    def __init__(self, organism, membership, matrix,
                 meme_suite,
                 seqtype='upstream',
                 sequence_filters=[],
                 pvalue_filter=None,
                 weight_func=None,
                 interval=0,
                 config_params=None):
        """creates a ScoringFunction"""
        MotifScoringFunctionBase.__init__(self, organism, membership,
                                          matrix, meme_suite, seqtype,
                                          sequence_filters, pvalue_filter,
                                          weight_func, interval,
                                          config_params)

    def name(self):
        """returns the name of this scoring function"""
        return "MEME"

    def compute_meme_run(self, seqs):
        """performs a MEME run"""
        return self.meme_suite.run_meme(seqs, self.used_seqs)


class WeederScoringFunction(MotifScoringFunctionBase):
    """Motif scoring function that runs Weeder instead of MEME"""

    def __init__(self, organism, membership, matrix,
                 meme_suite, seqtype,
                 sequence_filters=[], pvalue_filter=None,
                 weight_func=None, interval=0):
        """creates a scoring function"""
        MotifScoringFunctionBase.__init__(self, organism, membership, matrix,
                                          meme_suite, seqtype,
                                          sequence_filters, pvalue_filter,
                                          weight_func, interval, None)

    def name(self):
        """returns the name of this scoring function"""
        return "Weeder"

    def compute_meme_run(self, seqs):
        """computes the values of a meme run"""
        with tempfile.NamedTemporaryFile(prefix='weeder.fasta',
                                         delete=False) as outfile:
            filename = outfile.name
            logging.info("Run Weeder on FASTA file: '%s'", filename)
            st.write_sequences_to_fasta_file(outfile, seqs.items())

        pssms = weeder.run_weeder(filename)
        meme_outfile = '%s.meme' % filename
        dbfile = self.meme_suite.make_sequence_file(
            [(feature_id, locseq[1])
             for feature_id, locseq in self.used_seqs.items()])
        logging.info("# PSSMS created: %d", len(pssms))
        logging.info("run MAST on '%s'", meme_outfile)
        try:
            mast_out = self.meme_suite.mast(
                meme_outfile, dbfile,
                self.meme_suite.global_background_file())
            pe_values, annotations = meme.read_mast_output(mast_out,
                                                           seqs.keys())
            return meme.MemeRunResult(pe_values, annotations, [])
        except:
            return meme.MemeRunResult([], {}, [])
