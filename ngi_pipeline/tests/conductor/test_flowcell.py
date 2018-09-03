
import unittest

from ngi_pipeline.conductor.flowcell import match_fastq_sample_number_to_samplesheet


class TestFlowcell(unittest.TestCase):

    def test_match_fastq_sample_number_to_samplesheet(self):
        test_samples = {
            "S1": [["S1", "proj1", "sample-name", "", "", 1],
                   ["S1", "proj1", "sample-name", "", "", 2],
                   ["S1", "proj1", "sample-name", "", "", 3],
                   ["S1", "proj1", "sample-name", "", "", 4]],
            "S2": [["S2", "proj1", "sample-name", "", "", 1],
                   ["S2", "proj1", "sample-name", "", "", 2]]
        }
        test_data = [
            ["sample-name_S1_L001_R1_001.fastq.gz", test_samples["S1"]+test_samples["S2"], None, test_samples["S1"][0]],
            ["sample-name_S1_L001_R1_001.fastq.gz", test_samples["S1"]+test_samples["S2"], "not-proj1", None],
            ["sample-name_S1_L003_R1_001.fastq.gz", test_samples["S1"]+test_samples["S2"], None, test_samples["S1"][2]],
            ["sample-name_S2_L002_R1_001.fastq.gz", test_samples["S1"]+test_samples["S2"], None, test_samples["S2"][1]],
            ["not-sample-name_S2_L002_R1_001.fastq.gz", test_samples["S1"]+test_samples["S2"], None, None],
            ["sample-name_S2_L003_R1_001.fastq.gz", test_samples["S1"]+test_samples["S2"], None, None],
            ["sample-name_S3_L003_R1_001.fastq.gz", test_samples["S1"]+test_samples["S2"], None, None],
            ["sample-name_S3_L004_R1_001.fastq.gz", None, None, None],
            ["sample-name_S3_L004_R1_001.fastq.gz", [], None, None],
            ["sample-name_S3_L004_R1_001.fastq.gz", [[]], None, None],
            ["sample-name_SX_L004_R1_001.fastq.gz", test_samples["S1"] + test_samples["S2"], None, None]
        ]
        for test_item in test_data:
            self.assertEqual(
                test_item[3],
                match_fastq_sample_number_to_samplesheet(test_item[0], test_item[1], test_item[2]))
