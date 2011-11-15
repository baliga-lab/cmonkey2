"""set_enrichment.py - cMonkey set_enrichment scoring.

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import util
import scoring
import numpy as np
import rpy2.robjects as robjects


class EnrichmentSet:
    def __init__(self, cutoff):
        self.genes = []
        self.values = []
        self.cutoff = cutoff
        self.__genes_above_cutoff = None

    def add(self, elem, value):
        self.genes.append(elem)
        self.values.append(value)

    def genes_above_cutoff(self):
        if self.__genes_above_cutoff == None:
            self.__genes_above_cutoff = []
            for index in range(0, len(self.genes)):
                if self.cutoff == 'discrete':
                    self.__genes_above_cutoff.append(self.genes[index])
                elif self.values[index] >= self.cutoff:
                    self.__genes_above_cutoff.append(self.genes[index])
        return self.__genes_above_cutoff

    def __repr__(self):
        return "Cutoff = %s, elems: %s" % (self.cutoff, self.genes)


class SetType:
    def __init__(self, name, sets):
        self.name = name
        self.sets = sets
        self.__genes = None

    def genes(self):
        if self.__genes == None:
            self.__genes = set()
            for enrichment_set in self.sets.values():
                for gene in enrichment_set.genes:
                    self.__genes.add(gene)
        return self.__genes

    def __repr__(self):
        result = "SetType['%s'] = {\n" % self.name
        for key, value in self.sets.items():
            result +=  "%s -> %s\n\n" % (key, value)
        result += "}"
        return result

    @classmethod
    def read_csv(cls, name, infile, cutoff=None, sep=','):
        """reads a set from a CSV file"""
        dfile = util.DelimitedFile.read(infile, sep)
        sets = {}
        for line in dfile.lines():
            if line[0] not in sets:
                sets[line[0]] = EnrichmentSet('discrete')
            sets[line[0]].add(line[1].upper(), 1)
        return SetType(name, sets)


class ScoringFunction(scoring.ScoringFunctionBase):
    """Network scoring function"""

    def __init__(self, membership, matrix, set_types,
                 weight_func=None,
                 interval=0):
        """Create scoring function instance"""
        scoring.ScoringFunctionBase.__init__(self, membership,
                                             matrix, weight_func)
        self.__interval = interval
        self.__set_types = set_types

    def compute(self, iteration):
        if (self.__interval == 0 or
            (iteration > 0 and (iteration % self.__interval == 0))):
            result = {}
            for set_type in self.__set_types:

                for cluster in range(1, self.num_clusters() + 1):
                    cluster_rows = self.rows_for_cluster(cluster)
                    cluster_genes = [gene for gene in cluster_rows
                                     if gene in set_type.genes()]
                    phyper_q = []
                    phyper_m = []
                    phyper_n = []
                    phyper_k = []

                    for set_name, eset in set_type.sets.items():
                        set_genes = eset.genes_above_cutoff()
                        intersect = set(cluster_genes).intersection(set_genes)
                        num_overlap = len(intersect)
                        phyper_q.append(len(intersect))
                        phyper_m.append(len(set_genes))

                    num_sets = len(set_type.sets)
                    phyper_n = (np.array([len(set_type.genes())
                                         for _ in range(num_sets)]) -
                                np.array(phyper_m))
                    phyper_n = [value for value in phyper_n]
                    phyper_k = [len(cluster_genes) for _ in range(num_sets)]

                    print "Q = ", phyper_q
                    print "M = ", phyper_m
                    print "N = ", phyper_n
                    print "K = ", phyper_k
                    hyper_values = phyper(phyper_q, phyper_m, phyper_n, phyper_k)
                    print "HYPER: ", hyper_values
                    
                    #print "SETTYPE = %s, GENES CLUSTER %d = %s" % (set_type.name, cluster, str(cluster_genes))
                    # TODO: apply cutoff
            return None
        else:
            return None

def phyper(q, m, n, k):
    r_phyper = robjects.r['phyper']
    #kwargs = {'lower.tail': False}
    return r_phyper(q, m, n, k, **kwargs)
    return r_phyper(q, m, n, k)
