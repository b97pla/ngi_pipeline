import datetime
import mock
import os
import random
import tempfile
import unittest

from ngi_pipeline.utils.parsers import get_flowcell_id_from_dirtree, parse_lane_from_filename, \
                                       find_fastq_read_pairs, find_fastq_read_pairs_from_dir, \
                                       determine_library_prep_from_samplesheet, _get_and_trim_field_value, \
                                        _get_libprepid_from_description, get_sample_numbers_from_samplesheet
from ngi_pipeline.tests import generate_test_data as gtd

class TestCommon(unittest.TestCase):

    def test_get_flowcell_id_from_dirtree(self):
        date = datetime.date.today().strftime("%y%m%d")
        flowcell_id = gtd.generate_flowcell_id()
        run_id = gtd.generate_identifier(date=date, flowcell_id=flowcell_id)
        sthlm_project_path = "/proj/a2010002/analysis/{project_name}/{sample_name}/{identifier}".format(
                project_name = gtd.generate_project_name(),
                sample_name = gtd.generate_sample_name(),
                identifier = run_id)
        uppsala_project_path = "/proj/a2010002/analysis/{run_id}/Sample_{sample_name}".format(
                run_id = run_id,
                sample_name = gtd.generate_sample_name())
        # Test Sthlm format
        self.assertEqual(get_flowcell_id_from_dirtree(sthlm_project_path), flowcell_id)
        # Test Uppsala format
        self.assertEqual(get_flowcell_id_from_dirtree(uppsala_project_path), flowcell_id)

    def test_parse_lane_from_filename(self):
        lane = random.randint(1,8)
        file_name = gtd.generate_sample_file_name(lane=lane)
        self.assertEqual(lane, parse_lane_from_filename(file_name))

    def test_find_fastq_read_pairs(self):
        # Test list functionality
        file_list = [ "P123_456_AAAAAA_L001_R1_001.fastq.gz",
                      "P123_456_AAAAAA_L001_R2_001.fastq.gz",]
        expected_output = {"P123_456_AAAAAA_L001_": sorted(file_list) }
        self.assertEqual(expected_output, find_fastq_read_pairs(file_list))

    def test_find_fastq_read_pairs_from_dir(self):
        tmp_dir = tempfile.mkdtemp(suffix="_AFCID")
        file_list = map(lambda x: os.path.join(tmp_dir, x),
                    [ "P123_456_AAAAAA_L001_R1_001.fastq.gz",
                      "P123_456_AAAAAA_L001_R2_001.fastq.gz",])
        for file_name in file_list:
            # Touch the file
            open(os.path.join(tmp_dir, file_name), 'w').close()
        expected_output = {"P123_456_AAAAAA_L001_AFCID": file_list }
        produced_output = find_fastq_read_pairs_from_dir(tmp_dir)
        # This makes the test pass but is annoying
        produced_output = { k: sorted(v) for k,v in produced_output.items() }
        self.assertEqual(expected_output, produced_output)

    def test_determine_library_prep_from_samplesheet(self):
        ss_v18 = tempfile.mkstemp()[1]
        ss_v25 = tempfile.mkstemp()[1]
        samplesheet_v18_text = [
            "FCID,Lane,SampleID,SampleRef,Index,Description,Control,Recipe,Operator,SampleProject",
            "C45KVANXX,1,P1139_147,hg19,GTAGAGGA-TAGATCGC,Y__Mom_15_01,N,,Your Mother,Y__Mom_15_01",
            "C45KVANXX,2,P1139_145,hg19,AAGAGGCA-TAGATCGC,Y__Mom_15_01,N,,Your Mother,Y__Mom_15_01"]
        with open(ss_v18, 'w') as f:
            f.write("\n".join(samplesheet_v18_text))
        with self.assertRaises(ValueError): # Shouldn't find anything - raises ValueError
            libprep = determine_library_prep_from_samplesheet(samplesheet_path=ss_v18,
                                                              project_id="Y__Mom_15_01",
                                                              sample_id="P1139_147",
                                                              lane_num=1)
        samplesheet_v25_text = [
            "[Header],,,,,,,,",
            "Date,2014-02-12,,,,,,,",
            "[Data],,,,,,,,",
            "Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,Sample_Project,Description",
            "1,Sample_CEP-NA10860-PCR-free,CEP-NA10860-PCR-free,,,,AGTTCC,YM01,LIBRARY_NAME:CEP_Pool8",
            "2,Sample_CEP-NA11992-PCR-free,CEP-NA11992-PCR-free,,,,ATGTCA,YM01,LIBRARY_NAME:CEP_Pool9"]
        with open(ss_v25, 'w') as f:
            f.write("\n".join(samplesheet_v25_text))
        determined_libprep = determine_library_prep_from_samplesheet(samplesheet_path=ss_v25,
                                                                     project_id="YM01",
                                                                     sample_id="CEP-NA10860-PCR-free",
                                                                     lane_num=1)
        self.assertEqual("CEP_Pool8", determined_libprep)
        with self.assertRaises(ValueError):
            determined_libprep = determine_library_prep_from_samplesheet(samplesheet_path=ss_v25,
                                                                         project_id="YM01",
                                                                         sample_id="Sample_CEP-NA10860-PCR-free,",
                                                                         lane_num=2)

    def test__get_and_trim_field_value(self):
        test_dict = {
            "key1": "this-is-value1",
            "key2": "this-is-value2",
            "key3": "this-is-value3"
        }
        test_fields = sorted(test_dict.keys())
        self.assertEqual(_get_and_trim_field_value(test_dict, test_fields), test_dict[test_fields[0]])
        self.assertEqual(
            _get_and_trim_field_value(test_dict, ["key-does-not-match"] + test_fields),
            test_dict[test_fields[0]])
        self.assertEqual(
            _get_and_trim_field_value(test_dict, test_fields, "-is-"),
            test_dict[test_fields[0]].replace("-is-", ""))
        self.assertEqual(
            _get_and_trim_field_value(test_dict, test_fields, "string-does-not-match"),
            test_dict[test_fields[0]])
        self.assertIsNone(_get_and_trim_field_value(test_dict, ["key-does-not-match"]))

    def test__get_libprepid_from_description(self):
        self.assertIsNone(_get_libprepid_from_description("this;is;a;description"))
        self.assertEqual(_get_libprepid_from_description("this;is;a;description;LIBRARY_NAME:"), "")
        self.assertEqual(
            _get_libprepid_from_description("this:is;a:description;LIBRARY_NAME:this-is-the-libprepid;suffix"),
            "this-is-the-libprepid")

    def test_get_sample_numbers_from_samplesheet(self):
        samplesheet_header = "Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,Sample_Project," \
                             "Description".split(",")
        samplesheet_samples = [
            "1,Sample_CEP-NA10860-PCR-free,CEP-NA10860-PCR-free,,,,AGTTCC,YM01,LIBRARY_NAME:CEP_Pool8".split(","),
            "2,Sample_CEP-NA11992-PCR-free,CEP-NA11992-PCR-free,,,,ATGTCA,YM01,LIBRARY_NAME:CEP_Pool9".split(","),
            "3,Sample_CEP-NA10860-PCR-free,CEP-NA10860-PCR-free,,,,AGTTCC,YM01,LIBRARY_NAME:CEP_Pool8".split(","),
            "4,Sample_CEP-NA10860-PCR-free,CEP-NA10860-PCR-free,,,,GGTTCC,YM01,LIBRARY_NAME:CEP_Pool1".split(","),
            "4,Sample_CEP-NA10860-PCR-free,CEP-NA10860-PCR-free,,,,GGTTCC,YM01,LIBRARY_NAME:CEP_Pool11".split(",")
        ]
        samplesheet_rows = [dict(zip(samplesheet_header, sample)) for sample in samplesheet_samples]
        for row in samplesheet_rows[0:-1]:
            row["index2"] = "TGCGTT"

        with mock.patch('ngi_pipeline.utils.parsers.parse_samplesheet') as samplesheet_mock:
            samplesheet_mock.return_value = samplesheet_rows
            observed_sample_numbers = get_sample_numbers_from_samplesheet("samplesheet_path")
            self.assertListEqual(["S1", "S2", "S1", "S1", "S1"], [osn[0] for osn in observed_sample_numbers])
