import mock
import unittest
import tempfile
import os

from ngi_pipeline.engines.piper_ngi.launchers import analyze


class TestAnalyze(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.sthlm_data = {
            "proj_name": "Y.Mom_15_01",
            "proj_id": "P1155",
            "libprep_id": "P1155_prepA",
            "seqrun_id": "P1155_seqrunA",
            "sample_name": "P1155_101"}
        cls.upps_data = {
            "proj_name": "OK-9999",
            "proj_id": "OK-9999",
            "libprep_id": "OK-9999.SampleA.v1",
            "seqrun_id": "161019_ST-1234_0123_ATEST123CXX",
            "sample_name": "SampleA"}
        cls.engine_name = "piper_ngi"
        cls.proj_basepath = tempfile.mkdtemp()
        cls.workflow_name = "merge_process_variantcall"
        cls.xml_path = os.path.join(cls.proj_basepath, "some_config.xml")
        cls.exit_file = os.path.join(cls.proj_basepath, "some_file.exit")
        cls.project_obj = NGIProject(name=cls.proj_name,
                                     dirname=cls.proj_name,
                                     project_id=cls.proj_id,
                                     base_path=cls.proj_basepath)
        cls.sample_obj = cls.project_obj.add_sample(name=cls.sample_name,
                                                    dirname=cls.sample_name)
        # create a mock that can replace calls to charon
        cls.charon_mock = mock.Mock()
        cls.charon_mock.sample_get_libpreps = mock.Mock(return_value = {
            'libpreps': [{'qc': 'PASSED', 'libprepid': cls.libprep_id}]})
        cls.charon_mock.libprep_get_seqruns = mock.Mock(return_value = {
            'seqruns': [{'seqrunid': cls.seqrun_id}]})
        cls.charon_mock.seqrun_get = mock.Mock(return_value = {
            'alignment_status': ''})
        cls.charon_mock.project_get = mock.Mock(return_value={
            'sequencing_facility': 'Unknown'})

    def test_analyze(self):
        with mock.patch('ngi_pipeline.engines.piper_ngi.utils.CharonSession',
                        spec=False,
                        return_value=self.charon_mock) as dbmock:
            analyze(self.project_obj, self.sample_obj)
        return True
