import mock
import unittest

from ngi_pipeline.engines.sarek.database import CharonConnector, CharonError, TrackingConnector
from ngi_pipeline.engines.sarek.exceptions import BestPracticeAnalysisNotSpecifiedError, \
    SampleAnalysisStatusNotFoundError, SampleLookupError, AnalysisStatusForProcessStatusNotFoundError, \
    AlignmentStatusForAnalysisStatusNotFoundError, SampleAnalysisStatusNotSetError
from ngi_pipeline.engines.sarek.process import ProcessStopped
from ngi_pipeline.log.loggers import minimal_logger


@mock.patch("ngi_pipeline.database.classes.CharonSession", autospec=True)
class TestCharonConnector(unittest.TestCase):

    CONFIG = {}

    def setUp(self):
        self.log = minimal_logger(__name__, to_file=False, debug=True)
        self.config = TestCharonConnector.CONFIG
        self.project_id = "this-is-a-project-id"
        self.sample_id = "this-is-a-sample-id"
        self.libprep_id = "this-is-a-libprep-id"
        self.libpreps = [
            {"libprepid": "this-is-a-libprep-1"},
            {"libprepid": "this-is-a-libprep-2"},
            {"libprepid": "this-is-a-libprep-3"}]
        self.seqruns = [
            {"seqrunid": "this-is-a-seqrun-1"},
            {"seqrunid": "this-is-a-seqrun-2"},
            {"seqrunid": "this-is-a-seqrun-3"}]

    def _get_charon_connector(self, charon_session):
        self.charon_connector = CharonConnector(self.config, self.log, charon_session=charon_session)

    def test_best_practice_analysis(self, charon_session_mock):
        self._get_charon_connector(charon_session_mock.return_value)
        self.charon_connector.charon_session.project_get.return_value = {}
        with self.assertRaises(BestPracticeAnalysisNotSpecifiedError):
            self.charon_connector.best_practice_analysis(self.project_id)

        expected_bp_analysis = "this-is-a-best-practice-analysis"
        self.charon_connector.charon_session.project_get.return_value = {
            "best_practice_analysis": expected_bp_analysis}
        observed_bp_analysis = self.charon_connector.best_practice_analysis(self.project_id)

        self.assertEqual(expected_bp_analysis, observed_bp_analysis)

    def test_sample_analysis_status(self, charon_session_mock):
        self._get_charon_connector(charon_session_mock.return_value)
        self.charon_connector.charon_session.sample_get.return_value = {}
        with self.assertRaises(SampleAnalysisStatusNotFoundError):
            self.charon_connector.sample_analysis_status(self.project_id, self.sample_id)

        expected_sample_analysis_status = "this-is-a-sample-analysis-status"
        self.charon_connector.charon_session.sample_get.return_value = {
            "analysis_status": expected_sample_analysis_status}
        observed_sample_analysis_status = self.charon_connector.sample_analysis_status(
            self.project_id, self.sample_id)

        self.assertEqual(expected_sample_analysis_status, observed_sample_analysis_status)

    def libpreps_seqruns_helper(self, test_fn, test_args, expected_list, element_key):

        # simplest case, no restricting
        self.assertListEqual(
            expected_list,
            test_fn(*test_args))

        # filtering for a non-existing element should exclude everything
        self.assertListEqual(
            [],
            test_fn(*test_args, restrict_to=["this-element-does-not-exist"]))

        # filtering for existing elements should return only those
        self.assertListEqual(
            expected_list[0:2],
            test_fn(*test_args, restrict_to=map(lambda x: x[element_key], expected_list[0:2])))

        self.assertListEqual(
            expected_list[0:1],
            test_fn(*test_args, restrict_to=map(lambda x: x[element_key], expected_list[0:1])))

    def test_sample_libpreps(self, charon_session_mock):
        self._get_charon_connector(charon_session_mock.return_value)
        self.charon_connector.charon_session.sample_get_libpreps.return_value = {}

        with self.assertRaises(SampleLookupError) as sle:
            self.charon_connector.sample_libpreps(self.project_id, self.sample_id)

        expected_libpreps = self.libpreps
        self.charon_connector.charon_session.sample_get_libpreps.return_value = {"libpreps": expected_libpreps}
        self.libpreps_seqruns_helper(
            self.charon_connector.sample_libpreps, [self.project_id, self.sample_id], expected_libpreps, "libprepid")

    def test_libprep_seqruns(self, charon_session_mock):
        self._get_charon_connector(charon_session_mock.return_value)
        self.charon_connector.charon_session.libprep_get_seqruns.return_value = {}

        with self.assertRaises(SampleLookupError) as sle:
            self.charon_connector.libprep_seqruns(self.project_id, self.sample_id, self.libprep_id)

        expected_seqruns = self.seqruns
        self.charon_connector.charon_session.libprep_get_seqruns.return_value = {"seqruns": expected_seqruns}
        self.libpreps_seqruns_helper(
            self.charon_connector.libprep_seqruns,
            [self.project_id, self.sample_id, self.libprep_id],
            expected_seqruns,
            "seqrunid")

    def test_analysis_status_from_process_status(self, charon_session_mock):
        self._get_charon_connector(charon_session_mock.return_value)
        with self.assertRaises(AnalysisStatusForProcessStatusNotFoundError) as e:
            self.charon_connector.analysis_status_from_process_status(ProcessStopped)

        expected_status = CharonConnector._ANALYSIS_STATUS_FROM_PROCESS_STATUS.values()
        self.assertListEqual(
            expected_status,
            map(
                lambda p: self.charon_connector.analysis_status_from_process_status(p),
                CharonConnector._ANALYSIS_STATUS_FROM_PROCESS_STATUS.keys()))

    def test_alignment_status_from_analysis_status(self, charon_session_mock):
        self._get_charon_connector(charon_session_mock.return_value)
        with self.assertRaises(AlignmentStatusForAnalysisStatusNotFoundError) as e:
            self.charon_connector.alignment_status_from_analysis_status("this-status-does-not-exist")

        expected_status = CharonConnector._ALIGNMENT_STATUS_FROM_ANALYSIS_STATUS.values()
        self.assertListEqual(
            expected_status,
            map(
                lambda p: self.charon_connector.alignment_status_from_analysis_status(p),
                CharonConnector._ALIGNMENT_STATUS_FROM_ANALYSIS_STATUS.keys()))

    def test_set_sample_analysis_status(self, charon_session_mock):
        self._get_charon_connector(charon_session_mock.return_value)
        expected_return_value = "this-is-the-return-value"
        self.charon_connector.charon_session.sample_update.return_value = expected_return_value
        self.assertEqual(
            expected_return_value,
            self.charon_connector.set_sample_analysis_status(
                self.project_id, self.sample_id, "this-is-the-analysis-status"))

        self.charon_connector.charon_session.sample_update.side_effect = CharonError("raised CharonError")
        with self.assertRaises(SampleAnalysisStatusNotSetError) as e:
            self.charon_connector.set_sample_analysis_status(
                self.project_id, self.sample_id, "this-is-the-analysis-status")

        # set up for a recursive update
        self.charon_connector.charon_session.sample_update.side_effect = None
        self.charon_connector.charon_session.sample_get_libpreps.return_value = {"libpreps": self.libpreps}
        self.charon_connector.charon_session.libprep_get_seqruns.return_value = {"seqruns": self.seqruns}
        expected_libpreps = self.libpreps[-1].values()
        expected_seqruns = self.seqruns[1].values()
        self.charon_connector.set_sample_analysis_status(
            self.project_id,
            self.sample_id,
            CharonConnector._ALIGNMENT_STATUS_FROM_ANALYSIS_STATUS.keys()[0],
            recurse=True,
            restrict_to_libpreps=expected_libpreps,
            restrict_to_seqruns={lp.values()[0]: expected_seqruns for lp in self.libpreps})
        self.charon_connector.charon_session.seqrun_update.assert_called_once()


class TestTrackingConnector(unittest.TestCase):

    def test_pidfield_from_process_connector_type(self):
        self.assertListEqual(
            TrackingConnector.PIDFIELD_FROM_PROCESS_CONNECTOR_TYPE.values(),
            map(
                lambda c: TrackingConnector.pidfield_from_process_connector_type(c),
                TrackingConnector.PIDFIELD_FROM_PROCESS_CONNECTOR_TYPE.keys()))
