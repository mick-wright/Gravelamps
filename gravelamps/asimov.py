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
import subprocess
import time

import git.exc
from liquid import Liquid

from asimov import config
from asimov.pipeline import Pipeline, PipelineException, PipelineLogger
from asimov.utils import update

class Gravelamps(Pipeline):
    '''
    Gravelamps Pipeline Integration --- based primarily on the bilby Pipeline
    '''

    def __init__(self, production, category=None):
        super(Gravelamps, self).__init__(production, category)
        self.logger.info("Using the Gravelamps pipeline")

        if not production.pipeline.lower() == "gravelamps":
            raise PipelineException

    def detect_completion(self):
        '''
        Check for the production of the bilby posterior file to signal that the job has completed
        '''

        self.logger.info("Checking if the bilby parameter estimation has completed")
        results_dir = glob.glob(f"{self.production.rundir}/result")

        if len(results_dir) > 0:
            result_files = glob.glob(
                os.path.join(results_dir[0], "*merge*_result.hdf5"))
            result_files = glob.glob(
                os.path.join(results_dir[0], "*merge*_result.json"))
            self.logger.debug(f"Result Files: {result_files}")

            if len(result_files) > 0:
                self.logger.info("Result file(s) found, job completed")
                return True
            self.logger.info("Result file not found")
            return False
        self.logger.info("Result directory not found")
        return False

    def before_submit(self):
        '''
        Pre-submission hook, preserves relative paths
        '''

        self.logger.info("Running pre-submit hook")

        submission_files = glob.glob(f"{self.production.rundir}/submit/*.sub*")
        for sub_file in submission_files:
            if "dag" in sub_file:
                continue

            with open(sub_file, "r", encoding="utf-8") as file:
                original = file.read()
            with open(sub_file, "w", encoding="utf-8") as file:
                self.logger.info(f"Adding preserve relative_paths to {sub_file}")
                file.write("preserve_relative_paths = True\n" + original)

    def _determine_prior(self):
        '''
        Determines the correct choice of prior file for this production
        '''

        self.logger.info("Determining production prior file")

        if "prior file" in self.production.meta:
            self.logger.info("Production specifies prior file directly")
            self.logger.info("f{self.production.meta['prior file']}")
            return self.production.meta["prior file"]

        template = None

        if "event type" in self.production.meta:
            event_type = self.production.meta["event type"].lower()
        else:
            event_type = "bbh"
            self.production.meta["event type"] = event_type

            if self.production.event.issue_object:
                self.production.event.issue_object.update_data()

        if template is None:
            template_filename = f"{event_type}.prior.template"
            self.logger.info("[Gravelamps] Constructing a prior from template")

            try:
                template = os.path.join(
                    config.get("Gravelamps", "priors"), template_filename)
            except (configparser.NoOptionError, configparser.NoSectionError):
                from pkg_resources import resource_filename

                template = resource_filename("asimov", f"priors/{template_filename}")

        priors = {}
        priors = update(priors, self.production.event.ledger.data["priors"])
        priors = update(priors, self.production.event.meta["priors"])
        priors = update(priors, self.production.meta["priors"])

        liq = Liquid(template)
        rendered = liq.render(priors=priors, config=config)

        prior_name = f"{self.production.name}.prior"
        prior_file = os.path.join(os.getcwd(), prior_name)

        self.logger.info(f"Saving the new prior file as {prior_file}")

        with open(prior_file, "w", encoding="utf-8") as new_prior:
            new_prior.write(rendered)

        repo = self.production.event.repository
        try:
            repo.add_file(prior_file,
                          os.path.join(config.get("general", "calibration_dictionary"),
                                       prior_name))
            os.remove(prior_file)
        except git.exc.GitCommandError:
            pass

        return os.path.join(self.production.event.repository.dictionary,
                            config.get("general", "calibration_dictionary"),
                            prior_name)

    def build_dag(self, psds=None, user=None, clobber_psd=None, dryrun=False):
        '''
        Construct a DAG file in order to submit a production to the condor scheduler
        using gravelamps_inference.

        Parameters
        ----------
        production : str
            Production name
        psds : dict, optional
            The PSDs which should be used for this DAG. If no PSDs are provided
            the PSD files specified in the ini file will be used instead
        user : str
            The user accounting tag which should be used to run the job
        dryrun : bool
            If set to true the commands will not be run, but will be printed
            to standard output. Defaults to False

        Raises
        ------
        PipelineException
            Raised if the construction of the DAG fails
        '''

        cwd = os.getcwd()
        self.logger.info(f"Working in {cwd}")

        self._determine_prior()

        if self.production.event.repository:
            ini =\
                self.production.event.repository.find_prods(self.production.name, self.category)[0]
            ini = os.path.join(cwd, ini)
        else:
            ini = f"{self.production.name}.ini"

        if self.production.rundir:
            rundir = self.production.rundir
        else:
            rundir = os.path.join(os.path.expanduser("~"),
                                  self.production.event.name,
                                  self.production.name)
            self.production.rundir = rundir

        if "job label" in self.production.meta:
            job_label = self.production.meta["job label"]
        else:
            job_label = self.production.name

        command = [
            os.path.join(config.get("pipelines", "environment"), "bin", "gravelamps_inference"),
            ini]

        if dryrun:
            print(" ".join(command))
        else:
            self.logger.info(" ".join(command))
            pipe = subprocess.Popen(command,
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.STDOUT)
            out, err = pipe.communicate()
            self.logger.info(out)

            if err or "To submit, use the following command:" not in str(out):
                self.production.status = "stuck"
                self.logger.error(err)
                raise PipelineException(
                    f"DAG file could not be created.\n{command}\n{out}\n{err}",
                    production=self.production.name)

            time.sleep(10)
            return PipelineLogger(message=out, production=self.production.name)

    def submit_dag(self):

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
