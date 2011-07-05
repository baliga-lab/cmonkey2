"""cMonkey top-level module"""

class CMonkey:
    """
    The cMonkey object controls the overall execution of the cMonkey
    algorithm.
    This top-level class takes configuration and inputs to provide
    them for the actual execution
    configuration is a dictionary with the following keys:
    'organism'         - organism to use for computation
    'cog.organism'     - organism for NCBI COG
    'rsat.species'     - species for RSAT
    'num_iterations'   - # of iterations
    'clusters_per_row' - # of clusters/row
    ...
    """
    CMONKEY_VERSION = '4.0'
    RSAT_URLS = [
        'http://rsat.ccb.sickkids.ca/', 'http://rsat.ulb.ac.be/rsat/', 
        'http://embnet.ccg.unam.mx/rsa-tools'
        ]
    KEGG_FTP = 'ftp://ftp.genome.jp/pub/kegg/genes/taxonomy'

    def __init__(self, ratio_matrices, configuration=None):        
        """create a cMonkey object"""
        self.run_finished = False
        self.ratio_matrices = ratio_matrices
        if not configuration:
            self.configuration = CMonkey.make_default_config()
        else:
            self.configuration = configuration

    @classmethod
    def make_default_config(cls):
        """creates default configuration"""
        return {'organism': 'hpy', 'cog.organism': '?', 'rsat.species': '?',
                'num_iterations': 2000, 'clusters_per_row': 2,
                'operon.shift': True, 'background.order': 3,
                'recalc.background': True, 'row.scaling': 6,
                'seed.method.rows': 'kmeans', 'seed.method.cols': 'best',
                'post.adjust': True}

    def run(self):
        """start a run"""
        self.run_finished = True

class Membership:
    """Algorithms for cluster membership"""
    def __init__(self):
        """this class has only class methods"""
        pass

    @classmethod
    def map_to_is_member_matrix(cls, membership_matrix, kcluster):
        """maps a matrix containing row/column numbers to a true/false
        matrix by checking all values i in the range [1, kcluster] for
        containment in each row of the membership matrix.
        Example: mat =  [1 2] kcluster: 3
                        [2 3]
                 result = [t f] => 1 is in     [1 2], but not in [2 3]
                          [t t] => 2 is in     [1 2], and        [2 3]
                          [f t] => 3 is not in [1 2], but in     [2 3]
        """
        result = []
        for i in range(1, kcluster + 1):
            result_row = []
            result.append(result_row)
            for matrix_row in membership_matrix:
                result_row.append(i in matrix_row)

        return result
            

###########################################################################
##### CLUSTER SEEDING METHODS
####################################
class ClusterSeedingKMeans:
    """k-means cluster seeding"""

    def __init__(self, kcluster):
        """initialize seeding method"""
        self.kcluster = kcluster

    def run(self):
        """run a seed"""
        pass

class ClusterSeedingBest:
    """k-means cluster seeding"""

    def __init__(self, cluster_size):
        """initialize seeding method"""
        self.kcluster = kcluster

    def run(self):
        """run a seed"""
        pass
