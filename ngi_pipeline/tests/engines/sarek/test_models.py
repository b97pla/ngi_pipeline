import mock
import unittest

from ngi_pipeline.engines.sarek.exceptions import BestPracticeAnalysisNotRecognized, ReferenceGenomeNotRecognized
from ngi_pipeline.engines.sarek.models import SarekAnalysis, SarekGermlineAnalysis, ReferenceGenome
from ngi_pipeline.log.loggers import minimal_logger
from ngi_pipeline.tests.engines.sarek.test_launchers import TestLaunchers


@mock.patch("ngi_pipeline.engines.sarek.models.ReferenceGenome", autospec=True)
@mock.patch("ngi_pipeline.engines.sarek.database.CharonConnector", autospec=True)
class TestSarekAnalysis(unittest.TestCase):

    CONFIG = {}

    def setUp(self):
        self.log = minimal_logger(__name__, to_file=False, debug=True)
        self.config = TestSarekAnalysis.CONFIG

    def test_get_analysis_instance_for_project_germline(self, charon_connector_mock, reference_genome_mock):
        reference_genome_mock.get_instance.return_value = "this-is-a-reference-genome"
        charon_connector = charon_connector_mock.return_value
        charon_connector.best_practice_analysis.return_value = "something_germline_something"
        expected_analysis_class = SarekGermlineAnalysis
        observed_analysis_instance = SarekAnalysis.get_analysis_instance_for_project(
            "this-is-a-project-id",
            self.config,
            self.log,
            charon_connector=charon_connector)
        self.assertEqual(expected_analysis_class, type(observed_analysis_instance))

    def test_get_analysis_instance_for_project_somatic(self, charon_connector_mock, reference_genome_mock):
        reference_genome_mock.get_instance.return_value = "this-is-a-reference-genome"
        charon_connector = charon_connector_mock.return_value
        charon_connector.best_practice_analysis.return_value = "something_somatic_something"
        with self.assertRaises(NotImplementedError):
            SarekAnalysis.get_analysis_instance_for_project(
                "this-is-a-project-id",
                self.config,
                self.log,
                charon_connector=charon_connector)

    def test_get_analysis_instance_for_project_unknown(self, charon_connector_mock, reference_genome_mock):
        reference_genome_mock.get_instance.return_value = "this-is-a-reference-genome"
        charon_connector = charon_connector_mock.return_value
        charon_connector.best_practice_analysis.return_value = "something_unknown_something"
        with self.assertRaises(BestPracticeAnalysisNotRecognized):
            SarekAnalysis.get_analysis_instance_for_project(
                "this-is-a-project-id",
                self.config,
                self.log,
                charon_connector=charon_connector)

    def test_sample_should_be_started(self, charon_connector_mock, reference_genome_mock):

        def _echo_status(status, *args):
            return status

        charon_connector = charon_connector_mock.return_value
        charon_connector.sample_analysis_status.side_effect = _echo_status
        sarek_analysis = SarekAnalysis(
            reference_genome_mock.return_value, self.config, self.log, charon_connector=charon_connector)
        self.assertTrue(
            sarek_analysis.sample_should_be_started(
                "FAILED", "this-is-a-sample-id", restart_failed_jobs=True))
        self.assertFalse(
            sarek_analysis.sample_should_be_started(
                "FAILED", "this-is-a-sample-id", restart_failed_jobs=False))
        self.assertTrue(
            sarek_analysis.sample_should_be_started(
                "ANALYZED", "this-is-a-sample-id", restart_finished_jobs=True))
        self.assertFalse(
            sarek_analysis.sample_should_be_started(
                "ANALYZED", "this-is-a-sample-id", restart_finished_jobs=False))
        self.assertTrue(
            sarek_analysis.sample_should_be_started(
                "UNDER_ANALYSIS", "this-is-a-sample-id", restart_running_jobs=True))
        self.assertFalse(
            sarek_analysis.sample_should_be_started(
                "UNDER_ANALYSIS", "this-is-a-sample-id", restart_running_jobs=False))
        self.assertTrue(
            sarek_analysis.sample_should_be_started(
                "TO_BE_ANALYZED", "this-is-a-sample-id"))

    def test_configure_analysis(self, charon_connector_mock, reference_genome_mock):
        config = {
            "profile": "this-is-a-profile",
            "tools": "this-is-the-tools"}

        sarek_analysis = SarekAnalysis(
            reference_genome_mock.return_value,
            {"sarek": config},
            self.log,
            charon_connector=charon_connector_mock.return_value)

        self.assertEqual(config["profile"], sarek_analysis.profile)
        self.assertEqual(config["tools"], sarek_analysis.tools)
        self.assertEqual(SarekAnalysis.DEFAULT_CONFIG["nf_path"], sarek_analysis.nf_path)
        self.assertEqual(SarekAnalysis.DEFAULT_CONFIG["sarek_path"], sarek_analysis.sarek_path)

        config = {
            "profile": "this-is-another-profile",
            "nf_path": "this-is-the-path-to-nextflow",
            "sarek_path": "this-is-the-path-to-sarek"}
        profile, tools, nf_path, sarek_path = sarek_analysis.configure_analysis(config={"sarek": config})
        self.assertEqual(config["profile"], profile)
        self.assertEqual(config["nf_path"], nf_path)
        self.assertEqual(config["sarek_path"], sarek_path)
        self.assertEqual(SarekAnalysis.DEFAULT_CONFIG["tools"], tools)

    def test_command_line(self, charon_connector_mock, reference_genome_mock):

        sarek_analysis = SarekAnalysis(
            reference_genome_mock.return_value,
            self.config,
            self.log,
            charon_connector=charon_connector_mock.return_value)

        with self.assertRaises(NotImplementedError):
            sarek_analysis.command_line(None, None)


@mock.patch("ngi_pipeline.engines.sarek.models.ReferenceGenome", autospec=True)
@mock.patch("ngi_pipeline.engines.sarek.database.CharonConnector", autospec=True)
class TestSarekGermlineAnalysis(unittest.TestCase):

    CONFIG = {
        "profile": "this-is-a-profile",
        "tools": "this-is-the-tools",
        "nf_path": "this-is-the-nextflow-path",
        "sarek_path": "this-is-the-sarek-path",
        "sampleDir": "this-is-the-sample-dir",
        "outDir": "this-is-the-output-dir",
        "genome": "this-is-the-genome"}

    def setUp(self):
        self.log = minimal_logger(__name__, to_file=False, debug=True)
        self.config = TestSarekGermlineAnalysis.CONFIG

    def get_instance(self, charon_connector_mock, reference_genome_mock):
        reference_genome = reference_genome_mock.return_value
        reference_genome.__str__.return_value = self.config["genome"]
        return SarekGermlineAnalysis(
            reference_genome,
            {"sarek": self.config},
            self.log,
            charon_connector=charon_connector_mock)

    def test_command_line(self, charon_connector_mock, reference_genome_mock):
        sarek_analysis = self.get_instance(charon_connector_mock, reference_genome_mock)
        observed_cmd = sarek_analysis.command_line(self.config["sampleDir"], self.config["outDir"])

        self.assertIn("-profile {}".format(self.config["profile"]), observed_cmd)
        self.assertTrue(observed_cmd.startswith("{} run {}".format(self.config["nf_path"], self.config["sarek_path"])))
        for key in filter(lambda k: k not in ["profile", "nf_path", "sarek_path"], self.config.keys()):
            self.assertIn("--{}={}".format(key, self.config[key]), observed_cmd)

    def test_analyze_sample(self, charon_connector_mock, reference_genome_mock):
        with mock.patch("ngi_pipeline.engines.sarek.models.execute_command_line") as execute_mock, \
                mock.patch("ngi_pipeline.engines.sarek.models.safe_makedir") as makedir_mock:
            analysis_obj = TestLaunchers.get_NGIAnalysis(log=self.log)
            sarek_analysis = self.get_instance(charon_connector_mock, reference_genome_mock)
            for sample_obj in analysis_obj.project:
                sarek_analysis.analyze_sample(sample_obj, analysis_obj)
                makedir_mock.assert_called_once()
                execute_mock.assert_called_once()
                makedir_mock.reset_mock()
                execute_mock.reset_mock()


class TestReferenceGenome(unittest.TestCase):

    def test_get_instance(self):
        self.assertEqual("GRCh37", str(ReferenceGenome.get_instance("GRCh37")))
        self.assertEqual("GRCh38", str(ReferenceGenome.get_instance("GRCh38")))
        with self.assertRaises(ReferenceGenomeNotRecognized):
            ReferenceGenome.get_instance("unknown")
