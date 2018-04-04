import mock
import os
import unittest

from ngi_pipeline.conductor.classes import NGIProject
from ngi_pipeline.engines.sarek import local_process_tracking
from ngi_pipeline.engines.sarek.database import TrackingConnector
from ngi_pipeline.engines.sarek.process import ProcessRunning, ProcessStopped
from ngi_pipeline.tests.engines.sarek.test_launchers import TestLaunchers


@mock.patch("ngi_pipeline.engines.sarek.database.TrackingConnector", autospec=True)
@mock.patch("ngi_pipeline.engines.sarek.database.CharonConnector", autospec=True)
class TestLocalProcessTracking(unittest.TestCase):

    def test_remove_analysis(self, _, tracking_connector_mock):
        tracking_connector = tracking_connector_mock.return_value
        local_process_tracking.remove_analysis(
            "this-is-the-analysis", ProcessRunning, tracking_connector)
        tracking_connector.remove_analysis.assert_not_called()
        local_process_tracking.remove_analysis(
            "this-is-the-analysis", ProcessRunning, tracking_connector, force=True)
        tracking_connector.remove_analysis.assert_called_once()
        tracking_connector.remove_analysis.reset_mock()
        local_process_tracking.remove_analysis(
            "this-is-the-analysis", ProcessStopped, tracking_connector)
        tracking_connector.remove_analysis.assert_called_once()

    @mock.patch("ngi_pipeline.engines.sarek.models.SarekAnalysis.sample_analysis_exit_code_path")
    def helper_analysis_status(self, process_mock, analysis, exit_code_mock):
        analysis.project_id = "this-is-a-project-id"
        analysis.sample_id = "this-is-a-sample-id"
        analysis.project_base_path = "this-is-a-project-base-path"
        analysis.workflow = "SarekGermlineAnalysis"
        expected_exit_code_path = "this-is-a-path"
        exit_code_mock.return_value = expected_exit_code_path
        expected_status = ProcessRunning
        process_mock.return_value = expected_status
        self.assertEqual(
            expected_status,
            local_process_tracking.get_analysis_status(analysis))
        process_mock.assert_called_once()

    @mock.patch(
        "ngi_pipeline.engines.sarek.local_process_tracking.ProcessStatus.get_type_from_processid_and_exit_code_path")
    def test_get_analysis_status_process(self, process_mock, *mocks):
        expected_process_id = 12345
        expected_analysis = TrackingConnector._SampleAnalysis(process_id=expected_process_id)
        self.helper_analysis_status(process_mock, expected_analysis)

    @mock.patch(
        "ngi_pipeline.engines.sarek.local_process_tracking.JobStatus.get_type_from_processid_and_exit_code_path")
    def test_get_analysis_status_job(self, process_mock, *mocks):
        expected_process_id = 98765
        expected_analysis = TrackingConnector._SampleAnalysis(slurm_job_id=expected_process_id)
        self.helper_analysis_status(process_mock, expected_analysis)

    def test_update_charon_with_local_jobs_status(self, charon_connector_mock, tracking_connector_mock):
        expected_analyses = (
            TrackingConnector._SampleAnalysis(
                project_id="this-is-a-project-id-1",
                sample_id="this-is-a-sample-id",
                engine="this-is-the-engine",
                workflow="this-is-the-workflow",
                process_id=12345),
            TrackingConnector._SampleAnalysis(
                project_id="this-is-a-project-id-2",
                sample_id="this-is-a-sample-id",
                engine="this-is-the-engine",
                workflow="this-is-the-workflow",
                slurm_job_id=98765)
        )
        project_obj = NGIProject(None, None, None, None)
        project_obj.add_sample("this-is-a-sample-id", "this-is-a-sample-id")
        tracking_connector = tracking_connector_mock.return_value
        tracking_connector.tracked_analyses.return_value = expected_analyses
        with mock.patch.object(local_process_tracking, 'get_analysis_status') as analysis_status_mock, \
                mock.patch.object(local_process_tracking, 'remove_analysis') as remove_mock, \
                mock.patch.object(local_process_tracking, 'project_from_analysis') as project_mock:
            project_mock.return_value = project_obj
            local_process_tracking.update_charon_with_local_jobs_status(
                config={"this-is-not-an-empty-config": "nope"},
                log=None,
                tracking_connector=tracking_connector,
                charon_connector=charon_connector_mock.return_value)
            analysis_status_mock.assert_called()
            remove_mock.assert_called()

    def test__project_from_fastq_file_paths(self, *mocks):
        expected_project_obj = TestLaunchers.get_NGIProject("1")
        fastq_files = []
        for sample_obj in expected_project_obj:
            for libprep_obj in sample_obj:
                for seqrun_obj in libprep_obj:
                    fastq_files.extend(
                        map(
                            lambda f: os.path.join(
                                expected_project_obj.base_path,
                                "DATA",
                                expected_project_obj.dirname,
                                sample_obj.dirname,
                                libprep_obj.dirname,
                                seqrun_obj.dirname,
                                f),
                            seqrun_obj.fastq_files))
        observed_project_obj = local_process_tracking._project_from_fastq_file_paths(fastq_files)
        self.assertEqual(expected_project_obj, observed_project_obj)
