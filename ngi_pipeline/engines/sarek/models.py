import csv
import os
import re
from string import Template

from ngi_pipeline.engines.sarek.database import CharonConnector, TrackingConnector
from ngi_pipeline.engines.sarek.exceptions import BestPracticeAnalysisNotRecognized, ReferenceGenomeNotRecognized, \
    SampleNotValidForAnalysisError
from ngi_pipeline.engines.sarek.process import ProcessConnector, ProcessRunning, ProcessExitStatusSuccessful, \
    ProcessExitStatusFailed
from ngi_pipeline.utils.filesystem import safe_makedir


class SarekAnalysis(object):
    """
    Base class for the SarekAnalysis engine. This class contains the necessary methods for configuring and launching
    an analysis with the Sarek engine. However, some methods are not implemented (they are "abstract") and are
    expected to be implemented in subclasses, providing interfaces to the specialized analysis modes (e.g. Germline or
    Somatic).
    """

    DEFAULT_CONFIG = {
        "profile": "standard",
        "tools": None,
        "config": None,
        "sarek_path": "/lupus/ngi/production/latest/sw/sarek/",
        "nf_path": "/lupus/ngi/production/latest/sw/nextflow/nextflow",
        "containerPath": "/lupus/ngi/production/latest/resources/singularity/sarek"
    }

    def __init__(
            self,
            reference_genome,
            config,
            log,
            charon_connector=None,
            tracking_connector=None,
            process_connector=None):
        """
        Create an instance of SarekAnalysis.

        :param reference_genome: a type of ReferenceGenome indicating the reference genome to use
        :param config: a dict object with configuration options
        :param log: a log handle to use for logging
        :param charon_connector: a CharonConnector instance to use for the database connection. If not specified, a new
        connector will be created
        :param tracking_connector: a TrackingConnector instance to use for connections to the local tracking database.
        If not specified, a new connector will be created
        :param process_connector: a ProcessConnector instance to use for starting the analysis. If not specified, a new
        connector for local execution will be created
        """
        self.reference_genome = reference_genome
        self.config = config
        self.log = log
        self.sarek_config = self.configure_analysis(genome=self.reference_genome)
        self.charon_connector = charon_connector or CharonConnector(self.config, self.log)
        self.tracking_connector = tracking_connector or TrackingConnector(self.config, self.log)
        self.process_connector = process_connector or ProcessConnector(cwd=os.curdir)

    def __repr__(self):
        # returns the name of the instance type, e.g. "SarekAnalysis" or "SarekGermlineAnalysis"
        return type(self).__name__

    def configure_analysis(self, config=None, **opts):
        """
        Put together the Sarek config dict based on the default parameters in the class, updated with any options passed
        as well as the "sarek" section in a supplied config

        :param config: a config dict. If specified, any content stored under the "sarek" key will be included in the
        returned dict
        :param opts: additional options can be specified in the call and will be included in the returned dict
        :return: a config dict based on the default values and updated with the passed parameters
        """
        config = config or self.config
        sarek_config = self.DEFAULT_CONFIG.copy()
        sarek_config.update(opts)
        sarek_config.update(config.get("sarek", {}))
        return sarek_config

    @staticmethod
    def get_analysis_type_for_workflow(workflow):
        """
        Gets the type of the SarekAnalysis instance to use for a supplied workflow name.

        :param workflow: name of the workflow
        :return: a SarekAnalysis type appropriate for the workflow or None if no appropriate type exists
        """
        if workflow == "SarekAnalysis":
            return SarekAnalysis
        if workflow == "SarekGermlineAnalysis":
            return SarekGermlineAnalysis
        return None

    @staticmethod
    def get_analysis_instance_for_project(
            projectid,
            config,
            log,
            charon_connector=None,
            tracking_connector=None,
            process_connector=None):
        """
        Factory method returning a SarekAnalysis subclass instance corresponding to the best practice analysis specified
        for the supplied project.

        :param projectid: the projectid of the project to get an analysis instance for
        :param config: a config dict
        :param log: a log handle
        :param charon_connector: a connector instance to the charon database. If None, the default connector will be
        used
        :param tracking_connector: a TrackingConnector instance to use for connections to the local tracking database.
        If not specified, a new connector will be created
        :param process_connector: a ProcessConnector instance to use for starting the analysis. If not specified, a new
        connector for local execution will be created
        :return: an instance of a SarekAnalysis subclass
        """

        charon_connector = charon_connector or CharonConnector(config, log)
        tracking_connector = tracking_connector or TrackingConnector(config, log)
        process_connector = process_connector or ProcessConnector(cwd=os.curdir)

        # fetch the best practice analysis specified in Charon. This is expected to be a string formatted as:
        # ENGINE_MODE_GENOME, e.g. "Sarek_Germline_GRCh38"
        best_practice_analysis = charon_connector.best_practice_analysis(projectid)

        reference_genome = ReferenceGenome.get_instance(best_practice_analysis.split("_")[2])
        sarek_analysis_type = best_practice_analysis.split("_")[1].lower()
        if sarek_analysis_type == "germline":
            return SarekGermlineAnalysis(
                reference_genome,
                config,
                log,
                charon_connector,
                tracking_connector,
                process_connector)
        elif sarek_analysis_type == "somatic":
            raise NotImplementedError(
                "best-practice.analysis for {} is not implemented".format(best_practice_analysis))
        raise BestPracticeAnalysisNotRecognized(best_practice_analysis)

    def status_should_be_started(
            self,
            status,
            restart_failed_jobs=False,
            restart_finished_jobs=False,
            restart_running_jobs=False):
        """
        Takes a status string (e.g. the analysis_status or alignment_status as stored in Charon) and decides whether
        analysis should be started based on the status, taking the value of the restart flags into account.

        :param status: the status string as stored in Charon (e.g. "UNDER_ANALYSIS", "NOT RUNNING" etc.)
        :param restart_failed_jobs: if True, jobs marked as failed are ok to start (default is False)
        :param restart_finished_jobs: if True, jobs marked as finished are ok to start (default is False)
        :param restart_running_jobs: if True, jobs marked as running are ok to start (default is False)
        :return: True if the analysis is ok to start or False otherwise
        """
        def _charon_status_list_from_process_status(process_status):
            analysis_status = self.charon_connector.analysis_status_from_process_status(process_status)
            alignment_status = self.charon_connector.alignment_status_from_analysis_status(analysis_status)
            return analysis_status, alignment_status

        if status in _charon_status_list_from_process_status(ProcessRunning):
            return restart_running_jobs
        if status in _charon_status_list_from_process_status(ProcessExitStatusSuccessful):
            return restart_finished_jobs
        if status in _charon_status_list_from_process_status(ProcessExitStatusFailed):
            return restart_failed_jobs
        return True

    def sample_should_be_started(self,
                                 projectid,
                                 sampleid,
                                 restart_options):
        """
        Decides whether the analysis for a sample should be started based on the analysis status recorded in Charon
        and taking the value of the restart flags into account.

        :param projectid: the project id for the sample
        :param sampleid: the sample id
        :param restart_options: a dict with the restart options to take into account when deciding whether to start the
        sample
        :return: True if the analysis for the sample is ok to start or False otherwise
        """
        analysis_status = self.charon_connector.sample_analysis_status(projectid, sampleid)
        should_be_started = self.status_should_be_started(analysis_status, **restart_options)
        self.log.info(
            "{} - {}: sample analysis status is '{}' -> sample will{} be included".format(
                projectid,
                sampleid,
                analysis_status,
                "" if should_be_started else " NOT"))
        return should_be_started

    def libprep_should_be_started(self,
                                  projectid,
                                  sampleid,
                                  libprepid,
                                  start_failed_libpreps=False):
        """
        Decides whether the analysis for a libprep should be started based on the QC status recorded in Charon
        and taking the value of the start flags into account.

        :param projectid: the project id for the sample
        :param sampleid: the sample id
        :param libprepid: the libprep id
        :param start_failed_libpreps: if True, failed libpreps will be included in the analysis (default is False)
        :return: True if the analysis for the libprep is ok to start or False otherwise
        """
        libprep_qc_status = self.charon_connector.libprep_qc_status(projectid, sampleid, libprepid)
        should_be_started = libprep_qc_status != "FAILED" or start_failed_libpreps
        self.log.info(
            "{} - {} - {}: libprep QC is '{}' -> libprep will{} be included".format(
                projectid,
                sampleid,
                libprepid,
                libprep_qc_status,
                "" if should_be_started else " NOT"))
        return should_be_started

    def seqrun_should_be_started(self,
                                 projectid,
                                 sampleid,
                                 libprepid,
                                 seqrunid,
                                 restart_options):
        """
        Decides whether the analysis for a seqrun should be started based on the alignment status recorded in Charon
        and taking the value of the restart flags into account.

        :param projectid: the project id for the sample
        :param sampleid: the sample id
        :param libprepid: the libprep id
        :param seqrunid: the seqrun id
        :param restart_options: a dict with the restart options to take into account when deciding whether to start the
        seqrun
        :return: True if the analysis for the seqrun is ok to start or False otherwise
        """
        seqrun_alignment_status = self.charon_connector.seqrun_alignment_status(
            projectid, sampleid, libprepid, seqrunid)
        should_be_started = self.status_should_be_started(seqrun_alignment_status, **restart_options)
        self.log.info(
            "{} - {} - {} - {}: seqrun alignment status is '{}' -> seqrun will{} be included".format(
                projectid,
                sampleid,
                libprepid,
                seqrunid,
                seqrun_alignment_status,
                "" if should_be_started else " NOT"))
        return should_be_started

    def analyze_sample(self, sample_object, analysis_object):
        """
        Start the analysis for the supplied NGISample object and with the analysis details contained within the supplied
        NGIAnalysis object. If analysis is successfully started, will record the analysis in the local tracking
        database.

        Before starting, the status of the sample will be checked against the restart options in the analysis object.

        :raises: a SampleNotValidForAnalysisError if the sample is not eligible for analysis based on its status and the
        analysis options in the NGIAnalysis object
        :param sample_object: a NGISample object representing the sample to start analysis for
        :param analysis_object: a NGIAnalysis object containing the details for the analysis
        :return: None
        """
        projectid = analysis_object.project.project_id
        restart_options = {
            "restart_failed_jobs": analysis_object.restart_failed_jobs,
            "restart_finished_jobs": analysis_object.restart_finished_jobs,
            "restart_running_jobs": analysis_object.restart_running_jobs}

        if not self.sample_should_be_started(projectid, sample_object.name, restart_options):
            raise SampleNotValidForAnalysisError(projectid, sample_object.name, "nothing to analyze")

        args_to_collect_paths = [
            analysis_object.project.base_path,
            analysis_object.project.dirname,
            sample_object.name]

        # get the paths needed for the analysis
        sample_tsv_file = self.create_tsv_file(sample_object, analysis_object)
        sample_output_dir = self.sample_analysis_path(*args_to_collect_paths)
        sample_exit_code_path = self.sample_analysis_exit_code_path(*args_to_collect_paths)

        # get the command line to use for the analysis and execute it using the process connector
        cmd = self.command_line(sample_tsv_file, sample_output_dir)
        pid = self.process_connector.execute_process(
            cmd,
            working_dir=sample_output_dir,
            exit_code_path=sample_exit_code_path,
            job_name="{}-{}-{}".format(
                analysis_object.project.name,
                sample_object.name,
                str(self)))
        self.log.info("launched '{}', with {}, pid: {}".format(cmd, type(self.process_connector), pid))

        # record the analysis details in the local tracking database
        self.tracking_connector.record_process_sample(
            analysis_object.project.name,
            sample_object.name,
            analysis_object.project.base_path,
            str(self),
            "sarek",
            pid,
            type(self.process_connector))

    def command_line(self, sample_tsv_file, sample_output_dir):
        raise NotImplementedError("command_line should be implemented in the subclasses")

    @classmethod
    def sample_data_path(cls, base_path, projectid, sampleid):
        return os.path.join(base_path, "DATA", projectid, sampleid)

    @classmethod
    def sample_seqrun_path(cls, base_path, projectid, sampleid, libprepid, seqrunid):
        return os.path.join(cls.sample_data_path(base_path, projectid, sampleid), libprepid, seqrunid)

    @classmethod
    def sample_analysis_path(cls, base_path, projectid, sampleid):
        return os.path.join(base_path, "ANALYSIS", projectid, sampleid, cls.__name__)

    @classmethod
    def sample_analysis_exit_code_path(cls, base_path, projectid, sampleid):
        return os.path.join(
            cls.sample_analysis_path(base_path, projectid, sampleid),
            "{}-{}-{}.exit_code".format(projectid, sampleid, cls.__name__))

    @classmethod
    def sample_analysis_tsv_file(cls, base_path, projectid, sampleid):
        return os.path.join(
            cls.sample_analysis_path(base_path, projectid, sampleid),
            "{}-{}-{}.tsv".format(projectid, sampleid, cls.__name__))

    def generate_tsv_file_contents(self, sample_object, analysis_object):
        raise NotImplementedError("creation of sample tsv file contents should be implemented by subclasses")

    def create_tsv_file(self, sample_object, analysis_object):
        """
        Create a tsv file containing the information needed by Sarek for starting the analysis. Will decide the path to
        the tsv file based on the sample amd project information. If the path does not exist, it will be created.

        :raises: a SampleNotValidForAnalysisError if no libpreps or seqruns for the sample were eligible for analysis
        :param sample_object: a NGISample object representing the sample to create the tsv file for
        :param analysis_object: a NGIAnalysis object containing the details for the analysis
        :return: the path to the created tsv file
        """
        rows = self.generate_tsv_file_contents(sample_object, analysis_object)
        if not rows:
            raise SampleNotValidForAnalysisError(
                analysis_object.project.project_id,
                sample_object.name,
                "no libpreps or seqruns to analyze")

        tsv_file = self.sample_analysis_tsv_file(
            analysis_object.project.base_path,
            analysis_object.project.project_id,
            sample_object.name)
        safe_makedir(os.path.dirname(tsv_file))
        with open(tsv_file, "w") as fh:
            writer = csv.writer(fh, dialect=csv.excel_tab)
            writer.writerows(rows)
        return tsv_file


class SarekGermlineAnalysis(SarekAnalysis):
    """
    Class representing the Sarek Germline analysis mode. This inherits the SarekAnalysis class but any mode-specific
    methods and configurations are overriding the base class equivalents.
    """

    def command_line(self, sample_tsv_file, sample_output_dir):
        """
        Creates the command line for launching the sample analysis and returns it as a string.

        :param sample_tsv_file: the path to the tsv file needed by Sarek for the analysis
        :param sample_output_dir: the path to the folder where the sample analysis will be run
        :return: the command line as a string
        """
        # get the path to the nextflow executable and the sarek main script. If not present in the config, rely on the
        # environment to be aware of them
        step_args = [
            self.sarek_config.get("nf_path", "nextflow"),
            self.sarek_config.get("sarek_path", "sarek")]

        # each analysis step is represented by a SarekWorkflowStep instance
        processing_steps = [
            SarekPreprocessingStep(*step_args, sample=sample_tsv_file, outDir=sample_output_dir, **self.sarek_config),
            SarekGermlineVCStep(*step_args, outDir=sample_output_dir, **self.sarek_config),
            SarekAnnotateStep(*step_args, outDir=sample_output_dir, **self.sarek_config),
            SarekMultiQCStep(*step_args, outDir=sample_output_dir, **self.sarek_config)
        ]
        # create the command line by chaining the command lines from each processing step
        return " && ".join(
            map(lambda step: step.command_line(), processing_steps))

    def generate_tsv_file_contents(self, sample_object, analysis_object):
        """
        Create the contents of the tsv file used by Sarek for the analysis. This will check the libpreps and seqruns
        and decide whether to include them in the analysis based on QC flag and alignment status, respectively.

        :param sample_object: a NGISample object representing the sample to create the tsv contents for
        :param analysis_object: a NGIAnalysis object containing the details for the analysis
        :return: a list of lists representing tsv entries where the outer list is the rows and each element is a list of
        the fields that should be tab-separated in the tsv file
        """
        rows = []
        projectid = analysis_object.project.project_id
        restart_options = {
            "restart_failed_jobs": analysis_object.restart_failed_jobs,
            "restart_finished_jobs": analysis_object.restart_finished_jobs,
            "restart_running_jobs": analysis_object.restart_running_jobs}

        patientid = sample_object.name
        gender = "ZZ"
        status = 0
        sampleid = patientid
        # filter out libpreps that are not eligible for analysis
        for libprep in filter(
                lambda lp: self.libprep_should_be_started(projectid, sampleid, lp.name),
                sample_object):
            libprepid = libprep.name
            # filter out seqruns that are not eligible for analysis
            for seqrun in filter(
                    lambda sr: self.seqrun_should_be_started(
                        projectid, sampleid, libprepid, sr.name, restart_options),
                    libprep):
                # the runfolder is represented with a Runfolder object
                runfolder = Runfolder(
                    SarekAnalysis.sample_seqrun_path(
                        analysis_object.project.base_path,
                        projectid,
                        sampleid,
                        libprepid,
                        seqrun.name))
                flowcellid = runfolder.flowcell_id
                # iterate over the fastq file pairs belonging to the seqrun
                for sample_fastq_file_pair in SarekGermlineAnalysis._sample_fastq_file_pair(
                        map(lambda f: SampleFastq(os.path.join(runfolder.path, f)), seqrun.fastq_files)):
                    runid = "{}.{}.{}".format(
                        flowcellid, sample_fastq_file_pair[0].lane_number, sample_fastq_file_pair[0].sample_number)
                    rows.append(
                        [patientid, gender, status, sampleid, runid] +
                        map(
                            lambda x: x.path, sample_fastq_file_pair))
        return rows

    @staticmethod
    def fastq_files_from_tsv_file(tsv_file):
        """
        Get the path to the fastq files listed in a tsv file

        :param tsv_file: the path to the tsv file
        :return: a list of fastq file paths
        """
        fastq_files = []
        with open(tsv_file) as fh:
            reader = csv.reader(fh, dialect=csv.excel_tab)
            for sample in reader:
                fastq_files.extend(sample[5:])
        return fastq_files

    @staticmethod
    def _sample_fastq_file_pair(sample_fastqs):
        """
        Take a list of SampleFastq objects and return a generator where each element is a list containing a fastq file
        pair (if paired-end, for single-end the list will just contain one entry), with first element being R1 and
        second element being R2.

        :param sample_fastqs: a list of SampleFastq objects in any order
        :return: a generator of lists of fastq file pairs
        """
        sorted_fastqs = sorted(sample_fastqs, key=lambda f: (f.sample_number, f.lane_number, f.read_number))
        n = 0
        while n < len(sorted_fastqs):
            n += 1
            if int(sorted_fastqs[n-1].read_number) > 1:
                continue
            fastqs_to_yield = [sorted_fastqs[n-1]]
            try:
                if int(sorted_fastqs[n].read_number) > 1:
                    fastqs_to_yield.append(sorted_fastqs[n])
                    n += 1
            except IndexError:
                pass
            yield fastqs_to_yield


class SampleFastq(object):
    """
    SampleFastq represents a fastq file and takes care of identifying properties such as sample name, lane, read etc.
    by splitting the path based on an expected regexp.
    """
    # the regexp used for splitting the file name
    FASTQ_REGEXP = r'^(.+)_([^_]+)_L00(\d)_R(\d)[^\.]*(\..*)$'

    def __init__(self, path):
        """
        Creates a SampleFastq object representing a fastq file.

        :param path: the path to the fastq file, where the file name will be used to identify properties
        """
        self.path = path
        self.dirname = os.path.dirname(self.path)
        self.filename = os.path.basename(self.path)
        self.sample_name, \
            self.sample_number, \
            self.lane_number, \
            self.read_number, \
            self.file_extension = self.split_filename(self.filename)

    def split_filename(self, filename):
        """
        Split the filename according to the regexp and return the identified groups. If the regexp was not matched,
        None will be returned for all properties.

        :param filename: the filename to split by the regexp
        :return: a list of the values captured by the regexp or a list of None values if the regexp did not capture
        anything
        """
        match = re.match(self.FASTQ_REGEXP, filename)
        return match.groups() if match is not None else [None for i in range(5)]


class Runfolder(object):
    """
    Runfolder represents a runfolder and identifies properties such as run date, instrument id etc. based on the
    name element of the path and a regexp used to identify them.
    """
    # the regexp used to split the runfolder name
    RUNFOLDER_REGEXP = r'^(\d{6})_([^_]+)_(\d+)_([AB])(\w+)$'

    def __init__(self, path):
        """
        Creates a Runfolder object representing a runfolder directory.

        :param path: the path to the runfolder where the name will be used to identify properties
        """
        self.path = path
        self.dirname = os.path.dirname(self.path)
        self.runfolder_name = os.path.basename(self.path)
        self.run_date, \
            self.instrument_id, \
            self.run_number, \
            self.flowcell_position, \
            self.flowcell_id = self.split_runfolder_name(self.runfolder_name)

    def split_runfolder_name(self, runfolder_name):
        """
        Split the runfolder name according to the regexp and return the identified groups. If the regexp was not
        matched, None will be returned for all properties.

        :param runfolder_name: the runfolder name to split by the regexp
        :return: a list of the values captured by the regexp or a list of None values if the regexp did not capture
        anything
        """
        match = re.match(self.RUNFOLDER_REGEXP, runfolder_name)
        return match.groups() if match is not None else [None for i in range(5)]


class ReferenceGenome(object):
    """
    The ReferenceGenome class represents the reference genome used by Sarek for the analysis. It has a factory method
    for returning the correct instance from a string representation.
    """
    NAME = None

    def __repr__(self):
        return self.NAME

    @staticmethod
    def get_instance(genome_name):
        """
        Factory method to get a ReferenceGenome instance representing the genome from a string.

        :raises: ReferenceGenomeNotRecognized if the genome name was not recognized
        :param genome_name: the name of the reference as a string, e.g. "GRCh37" or "GRCh38"
        :return: a ReferenceGenome instance
        """
        try:
            return filter(
                lambda ref: str(ref).lower() == genome_name.lower(),
                (GRCh37(), GRCh38())
            )[0]
        except IndexError:
            raise ReferenceGenomeNotRecognized(genome_name)


class GRCh37(ReferenceGenome):
    """Class representing the GRCh37 reference genome"""
    NAME = "GRCh37"


class GRCh38(ReferenceGenome):
    """Class representing the GRCh38 reference genome"""
    NAME = "GRCh38"


class SarekWorkflowStep(object):
    """
    The SarekWorkflowStep class represents an analysis step in the Sarek workflow. Primarily, it provides a method for
    creating the step-specific command line.
    """

    available_tools = []

    def __init__(
            self,
            path_to_nextflow,
            path_to_sarek,
            **sarek_args):
        """
        Create a SarekWorkflowStep instance according to the passed parameters.

        :param path_to_nextflow: path to the nextflow executable
        :param path_to_sarek: path to the main Sarek folder
        :param sarek_args: additional Sarek parameters to be included on the command line
        """
        self.nf_path = path_to_nextflow
        self.sarek_path = path_to_sarek
        # create a dict with parameters based on the passed key=value arguments
        self.sarek_args = {k: v for k, v in sarek_args.items() if k not in ["nf_path", "sarek_path"]}
        # add/filter a tools parameter against the valid tools for the workflow step
        self.sarek_args["tools"] = self.valid_tools(self.sarek_args.get("tools", []))
        # expand any parameters passed as list items into a ","-separated string
        self.sarek_args = {k: v if type(v) is not list else ",".join(v) for k, v in self.sarek_args.items()}

    def _append_argument(self, base_string, name, hyphen="--"):
        """
        Append an argument with a placeholder for the value to the supplied string in a format suitable for the
        string.Template constructor. If no value exists for the argument name among this workflow step's config
        parameters, the supplied string is returned untouched.

        Example: step._append_argument("echo", "hello", "") should return "echo hello ${hello}", provided the step
        instance has a "hello" key in the step.sarek_args dict.

        :param base_string: the string to append an argument to
        :param name: the argument name to add a placeholder for
        :param hyphen: the hyphen style to prefix the argument name with (default "--")
        :return: the supplied string with an appended argument name and placeholder
        """
        # NOTE: a numeric value of 0 will be excluded (as will a boolean value of False)!
        if not self.sarek_args.get(name):
            return base_string
        return "{0} {2}{1} ${{{1}}}".format(base_string, name, hyphen)

    def command_line(self):
        """
        Generate the command line for launching this analysis workflow step. The command line will be built using the
        Sarek arguments passed to the step's constructor and returned as a string.

        :return: the command line for the workflow step as a string
        """
        single_hyphen_args = ["config", "profile"]
        template_string = "${nf_path} run ${sarek_step_path}"
        for argument_name in single_hyphen_args:
            template_string = self._append_argument(template_string, argument_name, hyphen="-")
        for argument_name in filter(lambda n: n not in single_hyphen_args, self.sarek_args.keys()):
            template_string = self._append_argument(template_string, argument_name, hyphen="--")
        command_line = Template(template_string).substitute(
            nf_path=self.nf_path,
            sarek_step_path=os.path.join(self.sarek_path, self.sarek_step()),
            **self.sarek_args)
        return command_line

    def sarek_step(self):
        raise NotImplementedError("The Sarek workflow step definition for {} has not been defined".format(type(self)))

    def valid_tools(self, tools):
        """
        Filter a list of tools against the list of available tools for the analysis step.

        :param tools: a list of tool names
        :return: a list of tool names valid for this workflow step
        """
        return list(filter(lambda t: t in self.available_tools, tools))


class SarekPreprocessingStep(SarekWorkflowStep):

    def sarek_step(self):
        return "main.nf"


class SarekGermlineVCStep(SarekWorkflowStep):

    available_tools = [
        "haplotypecaller",
        "strelka",
        "manta"
    ]

    def sarek_step(self):
        return "germlineVC.nf"


class SarekAnnotateStep(SarekWorkflowStep):

    available_tools = [
        "snpeff",
        "vep"
    ]

    def sarek_step(self):
        return "annotate.nf"


class SarekMultiQCStep(SarekWorkflowStep):

    def sarek_step(self):
        return "runMultiQC.nf"

