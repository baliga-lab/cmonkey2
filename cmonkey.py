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
        self.init_configuration(configuration)

        # run data
        self.stats       = None
        self.row_scores  = None
        self.col_scores  = None
        self.mot_scores  = None
        self.net_scores  = None
        self.r_scores    = None
        self.row_scaling = [0 for _ in range(self.num_iterations())]
        self.mot_scaling = [0 for _ in range(self.num_iterations())]
        self.net_scaling = [0 for _ in range(self.num_iterations())]
        self.fuzzy_index = [0 for _ in range(self.num_iterations())]


    def init_configuration(self, configuration):
        """creates default configuration"""
        if not configuration:
            configuration = {}

        # meaningful defaults
        configuration.setdefault('organism',          'hpy')
        configuration.setdefault('num_iterations',    2000)
        configuration.setdefault('clusters_per_row',  2)
        configuration.setdefault('operon.shift',      True)
        configuration.setdefault('background.order',  3)
        configuration.setdefault('recalc.background', True)
        configuration.setdefault('row.scaling',       6)
        configuration.setdefault('seed_method.rows',  'kmeans')
        configuration.setdefault('seed_method.cols',  'best')
        configuration.setdefault('post.adjust',       True)
        configuration.setdefault('verbose',           True)
        self.configuration = configuration        

    def is_verbose(self):
        """determine whether we are running in verbose mode"""
        return self.configuration['verbose']

    def num_iterations(self):
        """configured number of iterations"""
        return self.configuration.get('num_iterations', 2000)

    def run(self):
        """start a run"""
        self.seed_clusters()
        current_iteration = 0

        if self.is_verbose():
            print "# iterations: %d" % self.num_iterations()

        while current_iteration < self.num_iterations():
            self.iterate()
            current_iteration += 1

        self.run_finished = True

    def seed_clusters(self):
        """seed clusters using the selected row and column methods"""
        pass
    
    def iterate(self):
        """iteration step in cMonkey
        This is run over and over until no improvements can be achieved"""
        self.compute_all_scores()
        self.combine_scores()
        self.fuzzify_scores()       

    def compute_all_scores(self):
        """compute scores on microarray data and clusters"""
        self.compute_microarray_scores()
        self.compute_cluster_scores()

    def compute_microarray_scores(self):
        """compute scores on microarray data"""
        self.compute_row_scores()
        self.compute_column_scores()

    def compute_row_scores(self):
        """compute row scores on microarray data"""
        pass

    def compute_column_scores(self):
        """compute column scores on microarray data"""
        pass

    def compute_cluster_scores(self):
        """compute scores for clusters"""
        self.compute_meme_scores()
        self.compute_mot_scores()
        self.compute_net_scores()        

    def compute_meme_scores(self):
        """compute meme scores on clusters"""
        pass

    def compute_mot_scores(self):
        """computes mot.scores, using the output from the meme scoring"""
        pass

    def compute_net_scores(self):
        """compute net.scores from STRING, add weighted scores for other
        networks if they exist"""
        pass
    def combine_scores(self):
        """combine the computed scores"""
        pass

    def fuzzify_scores(self):
        """fuzzify scores a bit for stochasticity
        fuzz should be between 0.2 and 0 (decreasing with iter)"""
        pass

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
