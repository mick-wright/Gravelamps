'''
Gravelamps-Asimov Pipeline Integration

The following is based upon the generic work performed by Daniel Williams as well as Asimov
documentation to create a setup for Gravelamps as a workable pipeline within the Asimov
automation framework, allowing event handling to be automated

Written by Daniel Williams
           Mick Wright
'''

import configparser
import glob
import os
import re
import subprocess

from asimov import logging
from asimov.pipeline import Pipeline, PipelineException, PipelineLogger

class Gravelamps(Pipeline):
    '''
    Gravelamps Pipeline Integration --- based primarily on the bilby Pipeline
    '''

    def __init__(self, production, category=None):
        super(Gravelamps, self).__init__(production, category)
        self.logger = logging.AsimovLogger(event=production.event)
        if not production.pipeline.lower() == "gravelamps":
            raise PipelineException

    def detect_completion(self):
        '''
        Check for the production of the bilby posterior file to signal that the job has completed
        '''
        results_dir = glob.glob(f"{self.production.rundir}/result")
        if len(results_dir) > 0:
            if len(glob.glob(os.path.join(results_dir[0], "*merge_result.json"))) > 0:
                return True
        return False

    def build_dag(self):
        '''
        Build DAG using gravelamps_inference
        '''

        ini = f"{self.production.name}.ini"
        command = ["gravelamps_inference", ini]

        print(" ".join(command))
        self.logger.info(" ".join(command))
        with subprocess.Popen(command,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.STDOUT) as pipe:
            out, err = pipe.communicate()

            if err or "To submit, use the following command:" not in str(out):
                self.production.status = "stuck"
                raise PipelineException

        return PipelineLogger(message=out,
                              production=self.production.name)

    def submit_dag(self):
        '''
        Submit DAG to condor cluster
        '''
        self.before_submit()

        try:
            dag_filename = "gravelamps_inference.dag"
            command = ["condor_submit_dag",
                       "-batch-name",
                       f"gravelamps/{self.production.event.name}/{self.production.name}",
                       os.path.join(self.production.rundir, "submit", dag_filename)]
            dagman = subprocess.Popen(command,
                                      stdout=subprocess.PIPE,
                                      stderr=subprocess.STDOUT)
        except FileNotFoundError as error:
            raise PipelineException("Condor does not appear to be installed on this system.\n"
                                    f'''Command to be run: {" ".join(command)}.''') from error

        stdout, stderr = dagman.communicate()

        if "submitted to cluster" in str(stdout):
            cluster = re.search("submitted to cluster ([\d]+)", str(stdout)).groups()[0]
            self.production.status = "running"
            self.production.job_id = int(cluster)
            return cluster, PipelineLogger(stdout)

        raise PipelineException(f"The DAG file could not be submitted.\n\n{stdout}\n\n{stderr}",
                                issue=self.production.event.issue_object,
                                production=self.production.name)

    @classmethod
    def read_ini(cls, filepath):
        '''
        Read and parse Gravelamps configuration file.
        Gravelamps configuration files are INI compliant, with dedicated and important sections
        Individual options can be repeated between sections.
        '''

        config_parser = configparser.ConfigParser()
        config_parser.read(filepath)

        return config_parser
