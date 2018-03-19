import os

from ngi_pipeline.engines.sarek.database import CharonConnector, TrackingConnector
from ngi_pipeline.engines.sarek.exceptions import BestPracticeAnalysisNotRecognized, ReferenceGenomeNotRecognized
from ngi_pipeline.engines.sarek.process import ProcessConnector


class SarekAnalysis(object):
    """
    Base class for the SarekAnalysis engine
    """

    DEFAULT_CONFIG = {
        "profile": "standard",
        "tools": "haplotypecaller",
        "sarek_path": "/lupus/ngi/staging/wildwest/ngi2016001/private/CAW",
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
        sarek_config = config.get("sarek", {})
        profile = sarek_config.get("profile", self.DEFAULT_CONFIG.get("profile"))
        tools = sarek_config.get("tools", self.DEFAULT_CONFIG.get("tools"))
        nf_path = sarek_config.get("nf_path", self.DEFAULT_CONFIG.get("nf_path"))
        sarek_path = sarek_config.get("sarek_path", self.DEFAULT_CONFIG.get("sarek_path"))
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

    def sample_should_be_started(self,
                                 projectid,
                                 sampleid,
                                 restart_failed_jobs=False,
                                 restart_finished_jobs=False,
                                 restart_running_jobs=False):
        analysis_status = self.charon_connector.sample_analysis_status(projectid, sampleid)
        if analysis_status == "FAILED":
            return restart_failed_jobs
        if analysis_status == "ANALYZED":
            return restart_finished_jobs
        if analysis_status == "UNDER_ANALYSIS":
            return restart_running_jobs
        return True

    def analyze_sample(self, sample_object, analysis_object):
        args_to_collect_paths = [
            analysis_object.project.base_path,
            analysis_object.project.dirname,
            sample_object.dirname]
        sample_input_dir = self.sample_data_path(*args_to_collect_paths)
        sample_output_dir = self.sample_analysis_path(*args_to_collect_paths)
        sample_exit_code_path = self.sample_analysis_exit_code_path(*args_to_collect_paths)

        cmd = self.command_line(sample_input_dir, sample_output_dir)
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

    def command_line(self, sample_input_dir, sample_output_dir):
        raise NotImplementedError("command_line should be implemented in the subclasses")

    @classmethod
    def sample_data_path(cls, base_path, projectid, sampleid):
        return os.path.join(base_path, "DATA", projectid, sampleid)

    @classmethod
    def sample_analysis_path(cls, base_path, projectid, sampleid):
        return os.path.join(base_path, "ANALYSIS", projectid, sampleid, cls.__name__)

    @classmethod
    def sample_analysis_exit_code_path(cls, base_path, projectid, sampleid):
        return os.path.join(
            cls.sample_analysis_path(base_path, projectid, sampleid),
            "{}-{}-{}.exit_code".format(projectid, sampleid, cls.__name__))


class SarekGermlineAnalysis(SarekAnalysis):
    """
    Class representing the Sarek Germline analysis type
    """

    def command_line(self, sample_input_dir, sample_output_dir):
        cmd_template = \
            "{nf_path} run {sarek_path} " \
            "--tools={tools} " \
            "--genome={genome} " \
            "--sampleDir={sample_dir} " \
            "--outDir={output_dir} " \
            "-profile {profile}"
        return "sleep 60"
        return cmd_template.format(
            nf_path = self.nf_path,
            sarek_path=self.sarek_path,
            tools=self.tools,
            genome=str(self.reference_genome),
            sample_dir=sample_input_dir,
            output_dir=sample_output_dir,
            profile=self.profile)


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
