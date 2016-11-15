import mock
import unittest
import tempfile
import os

from ngi_pipeline.conductor.classes import NGIProject
from ngi_pipeline.database.classes import CharonError
from ngi_pipeline.engines.piper_ngi.launchers import analyze
from ngi_pipeline.tests.engines.piper_ngi.test_utils import ProjectData


class TestAnalyze(unittest.TestCase):

    def setUp(self):
        self.project_data = ProjectData()

    @classmethod
    def setUpClass(cls):
        # set up test data
        cls.engine_name = "piper_ngi"
        cls.workflow_name = "merge_process_variantcall"
        cls.proj_data = {
            "sthlm": {
                "proj_name": "Y.Mom_15_01",
                "proj_id": "P1155",
                "libprep_id": "P1155_prepA",
                "seqrun_id": "P1155_seqrunA",
                "sample_name": "P1155_101",
                "sequencing_facility": "NGI-S",
                "basepath": tempfile.mkdtemp()},
            "upps": {
                "proj_name": "OK-9999",
                "proj_id": "OK-9999",
                "libprep_id": "OK-9999.SampleA.v1",
                "seqrun_id": "161019_ST-1234_0123_ATEST123CXX",
                "sample_name": "SampleA",
                "sequencing_facility": "NGI-U",
                "basepath": tempfile.mkdtemp()}}
        cls.proj_obj = {
            site: NGIProject(
                name=data["proj_name"],
                dirname=data["proj_name"],
                project_id=data["proj_id"],
                base_path=data["basepath"]) for site, data in cls.proj_data.items()}
        cls.sample_obj = {
            site: obj.add_sample(
                name=cls.proj_data[site]["sample_name"],
                dirname=cls.proj_data[site]["sample_name"]) for site, obj in cls.proj_obj.items()}

        # create a mock that can replace calls to charon
        cls.charon_mock = {}
        for site, data in cls.proj_data.items():
            cls.charon_mock[site] = mock.Mock()
            cls.charon_mock[site].sample_get_libpreps = mock.Mock(
                return_value={
                    'libpreps': [{
                        'qc': 'PASSED',
                        'libprepid': data["libprep_id"]}]})
            cls.charon_mock[site].libprep_get_seqruns = mock.Mock(
                return_value={
                    'seqruns': [{
                        'seqrunid': data["seqrun_id"]}]})
            cls.charon_mock[site].seqrun_get = mock.Mock(
                return_value={
                    'alignment_status': ''})
            cls.charon_mock[site].project_get = mock.Mock(
                return_value={
                    'sequencing_facility': data["sequencing_facility"]})

    @mock.patch('ngi_pipeline.engines.piper_ngi.launchers.check_for_preexisting_sample_runs', side_effect=RuntimeError)
    def test_running_sample(self, *args):
        """ A sample that shouldn't be restarted should raise a RuntimeError """
        with self.assertRaises(RuntimeError):
            analyze("some", "value")

    @mock.patch('ngi_pipeline.engines.piper_ngi.launchers.check_for_preexisting_sample_runs', return_value=None)
    def test_exec_mode(self, *args):
        """ An unknown exec_mode should not be accepted """
        with self.assertRaises(ValueError):
            analyze("some", "value", exec_mode="this_is_not_a_permitted_value")

    @mock.patch('ngi_pipeline.engines.piper_ngi.launchers.load_modules', side_effect=RuntimeError)
    @mock.patch('ngi_pipeline.engines.piper_ngi.launchers.check_for_preexisting_sample_runs', return_value=None)
    def test_load_modules(self, *args):
        """ Failing to load modules should raise an execption """
        with self.assertRaises(RuntimeError):
            analyze("some", "value", exec_mode="local")

    @mock.patch('ngi_pipeline.engines.piper_ngi.launchers.check_for_preexisting_sample_runs', return_value=None)
    @mock.patch('ngi_pipeline.engines.piper_ngi.launchers.find_previous_genotype_analyses')
    def test_genotype_charon_error(self, gt_analyses_mock, *args):
        """ A Charon error for a sample should result in skipping of that sample """
        charon_mock = mock.Mock()
        charon_mock.sample_get = mock.Mock(side_effect=CharonError("mocked error"))
        with mock.patch('ngi_pipeline.engines.piper_ngi.launchers.CharonSession', return_value=charon_mock):
            for site in self.proj_obj.keys():
                self.assertIsNone(analyze(self.proj_obj[site], self.sample_obj[site], level="genotype"))
                charon_mock.sample_get.assert_called_once_with(
                    projectid=self.proj_obj[site].project_id,
                    sampleid=self.sample_obj[site].name)
                gt_analyses_mock.assert_not_called()
                charon_mock.sample_get.reset_mock()

    @mock.patch('ngi_pipeline.engines.piper_ngi.launchers.check_for_preexisting_sample_runs', return_value=None)
    @mock.patch('ngi_pipeline.engines.piper_ngi.launchers.find_previous_genotype_analyses')
    @mock.patch('ngi_pipeline.engines.piper_ngi.launchers.is_sample_analysis_running_local', side_effect=AssertionError)
    def test_genotype_restarting(self, run_local, gt_analyses_mock, *args):
        """ Whether to restart a sample or skip it should be handled appropriately """
        charon_mock = mock.Mock()
        charon_mock.sample_get = mock.Mock()
        with mock.patch('ngi_pipeline.engines.piper_ngi.launchers.CharonSession', return_value=charon_mock):
            for site in self.proj_obj.keys():
                # iterate over combinations that should mean that a sample is skipped or not skipped
                for prev_analyses, genotype_status, should_skip in (
                        (False, {"genotype_status": "DONE"}, True),
                        (True, {"genotype_status": "DONE"}, True),
                        (True, {"genotype_status": "FAILED"}, True),
                        (False, {"genotype_status": "FAILED"}, False)):
                    # configure the mocks to behave as desired
                    gt_analyses_mock.return_value = prev_analyses
                    charon_mock.sample_get.return_value = genotype_status
                    # also verify the force behavior
                    for force in (True, False):
                        if should_skip and not force:
                            # assert that the sample is skipped
                            try:
                                self.assertIsNone(analyze(
                                    self.proj_obj[site],
                                    self.sample_obj[site],
                                    level="genotype",
                                    restart_finished_jobs=force,
                                    restart_running_jobs=False))
                            except AssertionError as ae:
                                self.fail("sample was not skipped")
                        else:
                            with self.assertRaises(AssertionError):
                                analyze(
                                    self.proj_obj[site],
                                    self.sample_obj[site],
                                    level="genotype",
                                    restart_finished_jobs=force,
                                    restart_running_jobs=False)

