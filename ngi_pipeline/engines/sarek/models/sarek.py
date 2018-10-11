import csv
import os

from ngi_pipeline.engines.sarek.database import CharonConnector, TrackingConnector
from ngi_pipeline.engines.sarek.exceptions import BestPracticeAnalysisNotRecognized, SampleNotValidForAnalysisError
from ngi_pipeline.engines.sarek.models.resources import ReferenceGenome
from ngi_pipeline.engines.sarek.models.sample import SarekAnalysisSample
from ngi_pipeline.engines.sarek.models.workflow import SarekPreprocessingStep, SarekGermlineVCStep, SarekAnnotateStep, \
    SarekMultiQCStep
from ngi_pipeline.engines.sarek.parsers import ParserIntegrator
from ngi_pipeline.engines.sarek.process import ProcessConnector, ProcessRunning, ProcessExitStatusSuccessful, \
    ProcessExitStatusFailed
from ngi_pipeline.utils.filesystem import safe_makedir


class SarekAnalysis:
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
        "containerPath": "/lupus/ngi/production/latest/resources/containers/sarek"
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
    def get_analysis_instance_for_workflow(
            workflow,
            config,
            log,
            reference_genome=None,
            charon_connector=None,
            tracking_connector=None,
            process_connector=None):
        """
        Factory method returning a SarekAnalysis subclass instance corresponding to the best practice analysis specified
        by the supplied workflow name.

        :param workflow: name of the workflow
        :param config: a config dict
        :param log: a log handle
        :param reference_genome: an instance of ReferenceGenome. If omitted, it will be None
        :param charon_connector: a connector instance to the charon database. If None, the default connector will be
        used
        :param tracking_connector: a TrackingConnector instance to use for connections to the local tracking database.
        If not specified, a new connector will be created
        :param process_connector: a ProcessConnector instance to use for starting the analysis. If not specified, a new
        connector for local execution will be created
        :return: an instance of a SarekAnalysis subclass
        """
        instance_type = SarekAnalysis.get_analysis_type_for_workflow(workflow)
        return instance_type(
            reference_genome,
            config,
            log,
            charon_connector=charon_connector,
            tracking_connector=tracking_connector,
            process_connector=process_connector)

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

        analysis_sample = SarekAnalysisSample(
            analysis_object.project,
            sample_object,
            self,
            restart_options={
                "restart_failed_jobs": analysis_object.restart_failed_jobs,
                "restart_finished_jobs": analysis_object.restart_finished_jobs,
                "restart_running_jobs": analysis_object.restart_running_jobs})

        if not self.sample_should_be_started(
                analysis_sample.projectid, analysis_sample.sampleid, analysis_sample.restart_options):
            raise SampleNotValidForAnalysisError(
                analysis_sample.projectid, analysis_sample.sampleid, "nothing to analyze")

        # get the paths needed for the analysis
        self.create_tsv_file(analysis_sample)

        # get the command line to use for the analysis and execute it using the process connector
        cmd = self.command_line(analysis_sample)
        pid = self.process_connector.execute_process(
            cmd,
            working_dir=analysis_sample.sample_analysis_path(),
            exit_code_path=analysis_sample.sample_analysis_exit_code_path(),
            job_name="{}-{}-{}".format(
                analysis_object.project.name,
                sample_object.name,
                str(self)))
        self.log.info("launched '{}', with {}, pid: {}".format(cmd, type(self.process_connector), pid))

        # record the analysis details in the local tracking database
        self.tracking_connector.record_process_sample(
            analysis_sample.projectid,
            analysis_sample.sampleid,
            analysis_sample.project_base_path,
            str(self),
            "sarek",
            pid,
            type(self.process_connector))

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

    def command_line(self, analysis_sample):
        raise NotImplementedError("command_line should be implemented in the subclasses")

    def generate_tsv_file_contents(self, analysis_sample):
        raise NotImplementedError("creation of sample tsv file contents should be implemented by subclasses")

    def collect_analysis_metrics(self, analysis_sample):
        raise NotImplementedError("collection of analysis results should be implemented by subclasses")

    def create_tsv_file(self, analysis_sample):
        """
        Create a tsv file containing the information needed by Sarek for starting the analysis. Will decide the path to
        the tsv file based on the sample amd project information. If the path does not exist, it will be created.

        :raises: a SampleNotValidForAnalysisError if no libpreps or seqruns for the sample were eligible for analysis
        :param sample_object: a NGISample object representing the sample to create the tsv file for
        :param analysis_object: a NGIAnalysis object containing the details for the analysis
        :return: the path to the created tsv file
        """
        rows = self.generate_tsv_file_contents(analysis_sample)
        if not rows:
            raise SampleNotValidForAnalysisError(
                analysis_sample.projectid,
                analysis_sample.sampleid,
                "no libpreps or seqruns to analyze")

        tsv_file = analysis_sample.sample_analysis_tsv_file()
        safe_makedir(os.path.dirname(tsv_file))
        with open(tsv_file, "w") as fh:
            writer = csv.writer(fh, dialect=csv.excel_tab)
            writer.writerows(rows)
        return tsv_file

    @staticmethod
    def _sample_paths(project_base_path, projectid, sampleid, subroot=None):
        try:
            return os.path.join(project_base_path, subroot, projectid, sampleid)
        except AttributeError:
            return os.path.join(project_base_path, projectid, sampleid)

    @classmethod
    def sample_data_path(cls, *args):
        return cls._sample_paths(*args, subroot="DATA")

    @classmethod
    def sample_analysis_path(cls, *args):
        return os.path.join(cls._sample_paths(*args, subroot="ANALYSIS"), cls.__name__)

    @classmethod
    def _sample_analysis_file(cls, project_base_path, projectid, sampleid, extension):
        return os.path.join(
            cls.sample_analysis_path(project_base_path, projectid, sampleid),
            "{}-{}-{}.{}".format(
                projectid,
                sampleid,
                cls.__name__,
                extension))

    @classmethod
    def sample_analysis_exit_code_path(cls, *args):
        return cls._sample_analysis_file(*args, extension="exit_code")

    @classmethod
    def sample_analysis_tsv_file(cls, *args):
        return cls._sample_analysis_file(*args, extension="tsv")


class SarekGermlineAnalysis(SarekAnalysis):
    """
    Class representing the Sarek Germline analysis mode. This inherits the SarekAnalysis class but any mode-specific
    methods and configurations are overriding the base class equivalents.
    """

    def processing_steps(self, analysis_sample):
        # get the path to the nextflow executable and the sarek main script. If not present in the config, rely on the
        # environment to be aware of them
        step_args = [
            self.sarek_config.get("nf_path", "nextflow"),
            self.sarek_config.get("sarek_path", "sarek")]
        local_sarek_config = {"outDir": analysis_sample.sample_analysis_path()}
        local_sarek_config.update(self.sarek_config)
        return [SarekPreprocessingStep(
            *step_args, sample=analysis_sample.sample_analysis_tsv_file(), **local_sarek_config)] + map(
            lambda step: step(*step_args, **local_sarek_config),
            [SarekGermlineVCStep, SarekAnnotateStep, SarekMultiQCStep])

    def command_line(self, analysis_sample):
        """
        Creates the command line for launching the sample analysis and returns it as a string.

        :param sample_tsv_file: the path to the tsv file needed by Sarek for the analysis
        :param sample_output_dir: the path to the folder where the sample analysis will be run
        :return: the command line as a string
        """
        # each analysis step is represented by a SarekWorkflowStep instance
        # create the command line by chaining the command lines from each processing step
        return " && ".join(
            map(lambda step: step.command_line(), self.processing_steps(analysis_sample)))

    def generate_tsv_file_contents(self, analysis_sample):
        """
        Create the contents of the tsv file used by Sarek for the analysis. This will check the libpreps and seqruns
        and decide whether to include them in the analysis based on QC flag and alignment status, respectively.

        :param sample_object: a NGISample object representing the sample to create the tsv contents for
        :param analysis_object: a NGIAnalysis object containing the details for the analysis
        :return: a list of lists representing tsv entries where the outer list is the rows and each element is a list of
        the fields that should be tab-separated in the tsv file
        """
        rows = []
        patientid = analysis_sample.sampleid
        gender = "ZZ"
        status = 0
        sampleid = patientid
        for sample_fastq in analysis_sample.runid_and_fastq_files_for_sample():
            rows.append([patientid, gender, status, sampleid] + sample_fastq)
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
    def sample_fastq_file_pair(sample_fastqs):
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

    def collect_analysis_metrics(self, analysis_sample):
        results_parser = ParserIntegrator()
        for processing_step in self.processing_steps(analysis_sample):
            for parser_type, results_file in processing_step.report_files(analysis_sample):
                results_parser.add_parser(parser_type(results_file))
        return {
            metric: results_parser.query_parsers("get_{}".format(metric)).pop()
            for metric in ["percent_duplication", "autosomal_coverage", "total_reads"]}

