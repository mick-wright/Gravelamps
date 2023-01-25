"""Asimov Pipeline Integration

Following is the implementation of Gravelamps as an Asimov pipeline based on the generic
instructions provided by Daniel Williams within the Asimov documentation. This sets up the
configuration required to automate Gravelamps running within the Asimov framework, allowing event
handling to be automated.

Written by Daniel Williams,
           Mick Wright.
"""

import configparser
import glob
import os
import re
import subprocess
import time

import git.exc
from liquid import Liquid

from asimov import config
from asimov.pipeline import Pipeline, PipelineException, PipelineLogger, PESummaryPipeline
from asimov.utils import update

class Gravelamps(Pipeline):
    """
    Gravelamps specific Pipeline configuration.

    Based primarily upon Asimov's built-in Bilby Pipeline class. This handles building and
    submtiting individual event runs from a properly configured ledger.

    Methods
    -------
    detect_completion
        Assess if job has completed
    before_submit
        Pre submission hook
    build_dag
        Build Gravelamps DAG
    submit_dag
        Submits Gravelamps DAG to HTCondor
    collect_assets
        Collect result assets
    samples
        Collect result sample files for PESummary
    after_completion
        Post completion hook to run PESummary
    collect_logs
        Collect logs into dictionary
    check_progress
        Checks job progress

    See Also
    --------
    Asimov Documentation: for greater explanation of the concepts behind the class
    construction.
    """

    def __init__(self, production, category=None):
        super(Gravelamps, self).__init__(production, category)
        self.logger.info("Using the Gravelamps pipeline")

        if not production.pipeline.lower() == "gravelamps":
            raise PipelineException

    def detect_completion(self):
        """
        Assess if job has completed.

        The Gravelamps DAG's final job is always the bilby_pipe DAG. To assess if the DAG has
        completed therefore, the function checks for the existance of the final result file in the
        ouput directory to assert the completion of the job.

        Returns
        -------
        bool
            Job completion status---true if complete, false otherwise.
        """

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
        """
        Pre submission hook.

        The hook at present adds the preserve relative file path argument to the condor submission
        file.

        Notes
        -----
        The hook currently adds the results directory from bilby_pipe to the individual submission
        files that transfer input files. This is to deal with a current ongoing issue in bilby_pipe
        that is due to be fixed in the next release, and will be modified after this occurs.
        """

        self.logger.info("Running pre-submit hook")

        submission_files = glob.glob(f"{self.event.name}/{self.production.rundir}/submit/*.submit")
        result_dir = f"{self.production.rundir}/result"

        for sub_file in submission_files:
            if "dag" in sub_file:
                continue

            with open(sub_file, "r", encoding="utf-8") as file:
                original = file.readlines()
                #Implementing Daniel's fix --- replace when bilby_pipe fixes
                for idx, line in enumerate(original):
                    if line.startswith("transfer_input_files"):
                        original[idx] = line.replace("\n", f", {result_dir}\n")

            with open(sub_file, "w", encoding="utf-8") as file:
                self.logger.info(f"Adding preserve relative_paths to {sub_file}")
                file.write("preserve_relative_paths = True\n" + ''.join(original))

    def _determine_prior(self):
        """
        Determines the correct choice of prior file for this production
        """

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
        """
        Build Gravelamps DAG.

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
        """

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

            else:
                time.sleep(10)
                return PipelineLogger(message=out, production=self.production.name)

    def submit_dag(self, dryrun=False):
        """
        Submits DAG file to the condor cluster

        Parameters
        ----------
        dryrun : bool
            If set to true the DAG will not be submitted but all commands will be
            printed to standard output instead. Defaults to False

        Returns
        -------
        int
            The cluster ID assigned to the running DAG file
        PipelineLogger
            The pipeline logger message.

        Raises
        ------
        PipelineException
            This will be raised if the pipeline fails to submit the job
        """

        cwd = os.getcwd()
        self.logger.info(f"Working in {cwd}")

        self.before_submit()

        try:
            dag_filename = "gravelamps_inference.dag"

            command = ["condor_submit_dag",
                       "-batch-name",
                       f"gravelamps/{self.production.event.name}/{self.production.name}",
                       os.path.join(self.production.rundir, "submit", dag_filename)]

            if dryrun:
                print(" ".join(command))
            else:
                dagman = subprocess.Popen(command,
                                          stdout=subprocess.PIPE,
                                          stderr=subprocess.STDOUT)

                self.logger.info(" ".join(command))

                stdout, stderr = dagman.communicate()

                if "submitted to cluster" in str(stdout):
                    cluster = re.search(r"submitted to cluster ([\d]+)", str(stdout)).groups()[0]
                    self.logger.info(f"Submitted successfully. Running with job ID {int(cluster)}")

                    self.production.status = "running"
                    self.production.job_id = int(cluster)

                    return cluster, PipelineLogger(stdout)
                else:
                    self.logger.error("Could not submit the job to the cluster")
                    self.logger.info(stdout)
                    self.logger.error(stderr)

                    raise PipelineException("The DAG file could not be submitted")

        except FileNotFoundError as error:
            self.logger.exception(error)
            raise PipelineException("It appears HTCondor is not installed on this system.\n"
                                    f'''Command that would have been run: {" ".join(command)}.''')\
                                    from error

    def collect_assets(self):
        """
        Collect result assets.

        The current result assests are deemed to be the samples produced by the nested sampling
        run.
        """

        return {"samples": self.samples()}

    def samples(self, absolute=False):
        """
        Collect the combined samples file for PESummary

        Parameters
        ----------
        absolute : bool
            Flag to return the absolute or relative filepath

        Returns
        -------
        sample_files : str
            Path to the combined sample file
        """

        if absolute:
            rundir = os.path.abspath(self.production.rundir)
        else:
            rundir = self.production.rundir
        self.logger.info(f"Rundir for samples: {rundir}")

        sample_files = glob.glob(os.path.join(rundir, "result", "*_merge*_result.hdf5"))\
                       + glob.glob(os.path.join(rundir, "result", "*_merge*_result.json"))

        return sample_files

    def after_completion(self):
        """
        Post pipeline completion hook to run PESummary.

        The hook runs the PESummary pipeline to produce post-completion output.
        """
        post_pipeline = PESummaryPipeline(production=self.production)
        self.logger.info("Job has completed. Running PESummary")
        cluster = post_pipeline.submit_dag()
        self.production.meta["job id"] = int(cluster)
        self.production.status = "processing"
        self.production.event.update_data()

    def collect_logs(self):
        """
        Collect all of the log files which have been produced by this production and return
        their contents as a dictionary

        Returns
        -------
        messages : dict
            Dictionary containing the log file content or notification that the file could not
            be opened
        """

        logs = glob.glob(f"{self.production.rundir}/submit/*.err")\
               + glob.glob(f"{self.production.rundir}/log*/*.err")\
               + glob.glob(f"{self.production.rundir}/*/*.out")\
               + glob.glob(f"{self.production.rundir}/*.log")

        messages = {}

        for log in logs:
            try:
                with open(log, "r", encoding="utf-8") as log_f:
                    message = log_f.read()
                    message = message.split("\n")
                    messages[log.split("/")[-1]] = "\n".join(message[-100:])
            except FileNotFoundError:
                messages[log.split("/")[-1]] = "There was a problem opening this log file"

        return messages

    def check_progress(self):
        """
        Checks job progress.

        The job progress is checked up on by finding the number of iterations and the current value
        of the dlogz for the sampling runs. This combined information can be used to obtain a rough
        estimate of how far through the job the run is. This is returned in dictionary format.

        Returns
        -------
        messages : dict
            Dictionary containing job progress in the form of the number of iterations and current
            dlogz value. Will contain a message noting if the log file for the job could not be
            opened.
        """

        logs = glob.glob(f"{self.production.rundir}/log_data_analysis/*.out")

        messages = {}
        for log in logs:
            try:
                with open(log, "r", encoding="utf-8") as log_f:
                    message = log_f.read()
                    message = message.split("\n")[-1]
                    pat = re.compile(r"([\d]+)it")
                    iterations = pat.search(message)
                    pat = re.compile(r"dlogz:([\d]*\.[d]*)")
                    dlogz = pat.search(message)
                    if iterations:
                        messages[log.split("/")[-1]] = (iterations.group(),
                                                        dlogz.group())
            except FileNotFoundError:
                messages[log.split("/")[-1]] = "There was a problem opening this log file."

        return messages

    @classmethod
    def read_ini(cls, filepath):
        """
        Read and parse Gravelamps configuration file.

        Gravelamps configuration files are INI compliant, with dedicated and important sections
        Individual options can be repeated between sections.

        Returns
        -------
        config_parser : ConfigParser
            Object containing Gravelamps configuration settings based on the INI structure.
        """

        config_parser = configparser.ConfigParser()
        config_parser.read(filepath)

        return config_parser

    def html(self):
        out = ""
        return out

    def resurrect(self):
        """
        Attempt to ressurect a failed job.

        A failed job will be resurrected a maximum of five times assuming that a rescue DAG
        has been produced. 
        """

        try:
            count = self.production.meta["resurrections"]
        except KeyError:
            count = 0

        if count < 5 and\
           len(glob.glob(os.path.join(self.production.rundir, "submit", "*.rescue"))) > 0:
            count += 1
            self.submit_dag()
