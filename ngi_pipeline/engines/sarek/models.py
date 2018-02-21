import os

from ngi_pipeline.engines.sarek.database import CharonConnector
from ngi_pipeline.engines.sarek.exceptions import BestPracticeAnalysisNotRecognized, ReferenceGenomeNotRecognized
from ngi_pipeline.utils.filesystem import execute_command_line, safe_makedir


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

    def __init__(self, reference_genome, config, log, charon_connector=None, **kwargs):
        self.reference_genome = reference_genome
        self.config = config
        self.log = log
        self.profile, self.tools, self.nf_path, self.sarek_path = self.configure_analysis()
        self.charon_connector = charon_connector or CharonConnector(self.config, self.log)

    def __repr__(self):
        return "{}-{}".format(type(self).__name__, str(self.reference_genome))

    def configure_analysis(self, config=None):
        config = config or self.config
        sarek_config = config.get("sarek", {})
        profile = sarek_config.get("profile", self.DEFAULT_CONFIG.get("profile"))
        tools = sarek_config.get("tools", self.DEFAULT_CONFIG.get("tools"))
        nf_path = sarek_config.get("nf_path", self.DEFAULT_CONFIG.get("nf_path"))
        sarek_path = sarek_config.get("sarek_path", self.DEFAULT_CONFIG.get("sarek_path"))
        return profile, tools, nf_path, sarek_path

    @staticmethod
    def get_analysis_instance_for_project(projectid, config, log, charon_connector=None):
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

        if not charon_connector:
            charon_connector = CharonConnector(config, log)

        best_practice_analysis = charon_connector.best_practice_analysis(projectid)

        reference_genome = ReferenceGenome.get_instance(best_practice_analysis.split("_")[2])
        sarek_analysis_type = best_practice_analysis.split("_")[1].lower()
        if sarek_analysis_type == "germline":
            return SarekGermlineAnalysis(reference_genome, config, log, charon_connector=charon_connector)
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
        sample_input_dir = os.path.join(
            analysis_object.project.base_path,
            "DATA",
            analysis_object.project.dirname,
            sample_object.dirname)
        sample_output_dir = os.path.join(
            analysis_object.project.base_path,
            "ANALYSIS",
            analysis_object.project.dirname,
            sample_object.dirname,
            str(self))
        cmd = self.command_line(sample_input_dir, sample_output_dir)
        safe_makedir(sample_output_dir)
        proc = execute_command_line(cmd, shell=False, cwd=sample_output_dir)
        self.log.info("launched '{}', pid: {}".format(cmd, proc.pid))

    def command_line(self, sample_input_dir, sample_output_dir):
        raise NotImplementedError("command_line should be implemented in the subclasses")


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
