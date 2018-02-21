import mock
import os
import unittest

from ngi_pipeline.conductor.classes import NGIAnalysis, NGIProject, NGISample, NGILibraryPrep, NGISeqRun
from ngi_pipeline.log.loggers import minimal_logger


class TestLaunchers(unittest.TestCase):

    CONFIG = {
        "analysis": {
            "best_practice_analysis": {
                "sarek_germline_grch38": {
                    "analysis_engine": "ngi_pipeline.engines.sarek"
                }
            }
        }
    }

    @staticmethod
    def get_NGIObject(cls, n):
        name = "{}_{}".format(cls.__name__, n)
        return cls(name, "{}".format(name))

    @staticmethod
    def get_NGISeqRun(n):
        srun = TestLaunchers.get_NGIObject(NGISeqRun, "{}_{}_{}".format(n, n, n))
        srun.add_fastq_files(["{}_L001_R{}_001.fastq.gz".format(srun.name, i) for i in range(1, 3)])
        return srun

    @staticmethod
    def get_NGILibraryPrep(n):
        lprep = TestLaunchers.get_NGIObject(NGILibraryPrep, n)
        lprep._subitems = {srun.name: srun for srun in map(TestLaunchers.get_NGISeqRun, range(1, 3))}
        return lprep

    @staticmethod
    def get_NGISample(n):
        sample = TestLaunchers.get_NGIObject(NGISample, n)
        sample._subitems = {lprep.name: lprep for lprep in map(TestLaunchers.get_NGILibraryPrep, range(1, 3))}
        return sample

    @staticmethod
    def get_NGIProject(n):
        name = "{}_{}".format(NGIProject.__name__, n)
        project = NGIProject(
            name, "{}".format(name), "{}".format(name), os.path.join("/path", "to", name, "base"))
        project._subitems = {sample.name: sample for sample in map(TestLaunchers.get_NGISample, range(1, 3))}
        return project

    @staticmethod
    def get_NGIAnalysis(best_practice_analysis="sarek_germline_grch38", config=None, log=None):
        project = TestLaunchers.get_NGIProject(1)
        with mock.patch("ngi_pipeline.conductor.classes.CharonSession", autospec=True) as CharonSessionMock:
            charon_mock = CharonSessionMock.return_value
            charon_mock.project_get.return_value = {"best_practice_analysis": best_practice_analysis}
            return NGIAnalysis(
                project,
                config=(config or TestLaunchers.CONFIG),
                log=(log or minimal_logger(__name__, to_file=False, debug=True)))

    def setUp(self):
        self.config = TestLaunchers.CONFIG
        self.log = minimal_logger(__name__, to_file=False, debug=True)
        self.analysis_object = TestLaunchers.get_NGIAnalysis(config=self.config, log=self.log)

    def test_analyze(self):
        with mock.patch("ngi_pipeline.engines.sarek.models.SarekAnalysis", autospec=True) as SarekAnalysisMock:
            pass
