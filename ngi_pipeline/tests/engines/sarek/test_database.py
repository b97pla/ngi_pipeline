import mock
import unittest

from ngi_pipeline.engines.sarek.database import CharonConnector
from ngi_pipeline.engines.sarek.exceptions import BestPracticeAnalysisNotSpecifiedError, \
    SampleAnalysisStatusNotFoundError
from ngi_pipeline.log.loggers import minimal_logger


class TestCharonConnector(unittest.TestCase):

    CONFIG = {}

    def setUp(self):
        log = minimal_logger(__name__, to_file=False, debug=True)
        config = TestCharonConnector.CONFIG
        with mock.patch("ngi_pipeline.database.classes.CharonSession", autospec=True) as CharonSessionMock:
            charon_session = CharonSessionMock.return_value
            self.charon_connector = CharonConnector(config, log, charon_session=charon_session)

    def test_best_practice_analysis(self):
        self.charon_connector.charon_session.project_get.return_value = {}
        with self.assertRaises(BestPracticeAnalysisNotSpecifiedError):
            self.charon_connector.best_practice_analysis("this-is-a-project-id")

        expected_bp_analysis = "this-is-a-best-practice-analysis"
        self.charon_connector.charon_session.project_get.return_value = {
            "best_practice_analysis": expected_bp_analysis}
        observed_bp_analysis = self.charon_connector.best_practice_analysis("this-is-a-project-id")

        self.assertEqual(expected_bp_analysis, observed_bp_analysis)

    def test_sample_analysis_status(self):
        self.charon_connector.charon_session.sample_get.return_value = {}
        with self.assertRaises(SampleAnalysisStatusNotFoundError):
            self.charon_connector.sample_analysis_status("this-is-a-project-id", "this-is-a-sample-id")

        expected_sample_analysis_status = "this-is-a-sample-analysis-status"
        self.charon_connector.charon_session.sample_get.return_value = {
            "analysis_status": expected_sample_analysis_status}
        observed_sample_analysis_status = self.charon_connector.sample_analysis_status(
            "this-is-a-project-id", "this-is-a-sample-id")

        self.assertEqual(expected_sample_analysis_status, observed_sample_analysis_status)
