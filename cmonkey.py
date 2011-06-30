# The cMonkey object controls the overall execution of the cMonkey
# algorithm
# This top-level class takes configuration and inputs to provide
# them for the actual execution
class cMonkey:
  def __init__(self, ratios, organism='hpy'):
    self.run_finished = False
    self.organism = organism

  def run(self):
    # TODO: invoke algorithm
    self.run_finished = True
