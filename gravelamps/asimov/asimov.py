'''
Gravelamps-Asimov Pipeline Integration

The following is based upon the generic work performed by Daniel Williams as well as Asimov
documentation to create a setup for Gravelamps as a workable pipeline within the Asimov
automation framework, allowing event handling to be automated

Written by Daniel Williams
           Mick Wright
'''

from asimov.pipeline import Pipeline

class Gravelamps(Pipeline):
    '''
    Gravelamps Pipeline Integration --- based primarily on the bilby Pipeline
    '''

    def __init__(self, production):
        super().__init__(production)

    def detect_completion(self):
        pass

    def collect_assets(self):
        pass

    def collect_logs(self):
        pass

    def build_dag(self):
        '''
        Build DAG using gravelamps_inference
        '''

    def submit_dag(self):
        pass

    def before_submit(self, dryrun=False):
        pass

    def after_completion(self):
        pass
