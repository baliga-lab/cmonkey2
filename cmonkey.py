"""cMonkey top-level module"""
class CMonkey:
    """
    The cMonkey object controls the overall execution of the cMonkey
    algorithm.
    This top-level class takes configuration and inputs to provide
    them for the actual execution
    """
    def __init__(self, ratios, organism='hpy'):
        """create a cMonkey object"""
        self.run_finished = False
        self.organism = organism

    def run(self):
        """start a run"""
        self.run_finished = True
