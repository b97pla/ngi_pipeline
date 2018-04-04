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
    Base class for the SarekAnalysis engine
    """

    DEFAULT_CONFIG = {
        "profile": "standard",
        "tools": None,
        "sarek_path": "/lupus/ngi/staging/latest/sw/sarek/",
        "nf_path": "/lupus/ngi/staging/latest/sw/nextflow/nextflow"
    }

    def __init__(
            self,
            reference_genome,
            config,
            log,
            charon_connector=None,
            tracking_connector=None,
            process_connector=None):
        self.reference_genome = reference_genome
        self.config = config
        self.log = log
        self.profile, self.tools, self.nf_path, self.sarek_path = self.configure_analysis()
        self.charon_connector = charon_connector or CharonConnector(self.config, self.log)
        self.tracking_connector = tracking_connector or TrackingConnector(self.config, self.log)
        self.process_connector = process_connector or ProcessConnector(cwd=os.curdir)

    def __repr__(self):
        return type(self).__name__

    def configure_analysis(self, config=None):
        config = config or self.config
        sarek_config = self.DEFAULT_CONFIG.copy()
        sarek_config.update(config.get("sarek", {}))
        profile = sarek_config.get("profile")
        tools = sarek_config.get("tools")
        nf_path = sarek_config.get("nf_path")
        sarek_path = sarek_config.get("sarek_path")
        return profile, tools, nf_path, sarek_path

    @staticmethod
    def get_analysis_type_for_workflow(workflow):
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
        Factory method returning a SarekAnalysis subclass corresponding to the best practice analysis specified for
        the supplied project.

        :param projectid: the projectid of the project to get an analysis instance for
        :param config: a config dict
        :param log: a log handle
        :param charon_connector: a connector instance to the charon database. If None, the default connector will be
        used
        :return: an instance of a SarekAnalysis subclass
        """

        charon_connector = charon_connector or CharonConnector(config, log)
        tracking_connector = tracking_connector or TrackingConnector(config, log)
        process_connector = process_connector or ProcessConnector(cwd=os.curdir)

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

        sample_tsv_file = self.create_tsv_file(sample_object, analysis_object)
        sample_output_dir = self.sample_analysis_path(*args_to_collect_paths)
        sample_exit_code_path = self.sample_analysis_exit_code_path(*args_to_collect_paths)

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
        self.tracking_connector.record_process_sample(
            analysis_object.project.name,
            sample_object.name,
            analysis_object.project.base_path,
            str(self),
            "sarek",
            pid,
            TrackingConnector.pidfield_from_process_connector_type(
                type(self.process_connector)))

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
    Class representing the Sarek Germline analysis type
    """

    def command_line(self, sample_tsv_file, sample_output_dir):
        step_args = [
            self.nf_path,
            self.sarek_path]
        step_kwargs = {
            "profile": self.profile,
            "sample": sample_tsv_file,
            "outdir": sample_output_dir,
            "genome": self.reference_genome,
            "tools": self.tools
        }
        processing_steps = [
            SarekPreprocessingStep(*step_args, **step_kwargs),
            SarekGermlineVCStep(*step_args, **step_kwargs),
            SarekMultiQCStep(*step_args, **step_kwargs)
        ]
        return " && ".join(
            map(lambda step: step.command_line(), processing_steps))

    def generate_tsv_file_contents(self, sample_object, analysis_object):
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
        for libprep in filter(
                lambda lp: self.libprep_should_be_started(projectid, sampleid, lp.name),
                sample_object):
            libprepid = libprep.name
            for seqrun in filter(
                    lambda sr: self.seqrun_should_be_started(
                        projectid, sampleid, libprepid, sr.name, restart_options),
                    libprep):
                runfolder = Runfolder(
                    SarekAnalysis.sample_seqrun_path(
                        analysis_object.project.base_path,
                        projectid,
                        sampleid,
                        libprepid,
                        seqrun.name))
                flowcellid = runfolder.flowcell_id
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
        fastq_files = []
        with open(tsv_file) as fh:
            reader = csv.reader(fh, dialect=csv.excel_tab)
            for sample in reader:
                fastq_files.extend(sample[5:])
        return fastq_files

    @staticmethod
    def _sample_fastq_file_pair(sample_fastqs):
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

    FASTQ_REGEXP = r'^(.+)_([^_]+)_L00(\d)_R(\d)[^\.]*(\..*)$'

    def __init__(self, path):
        self.path = path
        self.dirname = os.path.dirname(self.path)
        self.filename = os.path.basename(self.path)
        self.sample_name, \
            self.sample_number, \
            self.lane_number, \
            self.read_number, \
            self.file_extension = self.split_filename(self.filename)

    def split_filename(self, filename):
        match = re.match(self.FASTQ_REGEXP, filename)
        return match.groups() if match is not None else [None for i in range(5)]


class Runfolder(object):

    RUNFOLDER_REGEXP = r'^(\d{6})_([^_]+)_(\d+)_([AB])(\w+)$'

    def __init__(self, path):
        self.path = path
        self.dirname = os.path.dirname(self.path)
        self.runfolder_name = os.path.basename(self.path)
        self.run_date, \
            self.instrument_id, \
            self.run_number, \
            self.flowcell_position, \
            self.flowcell_id = self.split_runfolder_name(self.runfolder_name)

    def split_runfolder_name(self, runfolder_name):
        match = re.match(self.RUNFOLDER_REGEXP, runfolder_name)
        return match.groups() if match is not None else [None for i in range(5)]


class ReferenceGenome(object):
    NAME = None

    def __repr__(self):
        return self.NAME

    @staticmethod
    def get_instance(genome_name):
        if genome_name.lower() == "grch37":
            return GRCh37()
        if genome_name.lower() == "grch38":
            return GRCh38()
        raise ReferenceGenomeNotRecognized(genome_name)


class GRCh37(ReferenceGenome):
    NAME = "GRCh37"


class GRCh38(ReferenceGenome):
    NAME = "GRCh38"


class SarekWorkflowStep(object):

    def __init__(self, nf_path, sarek_path, profile=None, sample=None, step=None, tools=None, genome=None, outdir=None):
        self.nf_path = nf_path
        self.sarek_path = sarek_path
        self.profile = profile
        self.sample = sample
        self.step = step
        self.tools = tools
        self.genome = genome
        self.outDir = outdir

    def _append_argument(self, base_string, name, hyphen="--"):
        if getattr(self, name) is None:
            return base_string
        return "{0} {2}{1} ${{{1}}}".format(base_string, name, hyphen)

    def command_line(self):
        template_string = self._append_argument("${nf_path} run ${sarek_step_path}", "profile", hyphen="-")
        for argument_name in ["sample", "step", "tools", "genome", "outDir"]:
            template_string = self._append_argument(template_string, argument_name)
        command_line = Template(template_string).substitute(
            nf_path=self.nf_path,
            sarek_step_path=os.path.join(self.sarek_path, self.sarek_step()),
            profile=self.profile,
            sample=self.sample,
            step=self.step,
            tools=",".join(self.tools or []),
            genome=self.genome,
            outDir=self.outDir)
        return command_line

    @classmethod
    def sarek_step(cls):
        raise NotImplementedError("The Sarek workflow step definition for {} has not been defined".format(cls))


class SarekPreprocessingStep(SarekWorkflowStep):

    @classmethod
    def sarek_step(cls):
        return "main.nf"


class SarekGermlineVCStep(SarekWorkflowStep):

    @classmethod
    def sarek_step(cls):
        return "germlineVC.nf"


class SarekMultiQCStep(SarekWorkflowStep):

    @classmethod
    def sarek_step(cls):
        return "runMultiQC.nf"

