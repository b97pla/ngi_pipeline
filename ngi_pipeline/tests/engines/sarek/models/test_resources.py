import os
import unittest

from ngi_pipeline.engines.sarek.exceptions import ReferenceGenomeNotRecognized
from ngi_pipeline.engines.sarek.models.resources import ReferenceGenome, Runfolder, SampleFastq


class TestReferenceGenome(unittest.TestCase):

    def test_get_instance(self):
        self.assertEqual("GRCh37", str(ReferenceGenome.get_instance("GRCh37")))
        self.assertEqual("GRCh38", str(ReferenceGenome.get_instance("GRCh38")))
        with self.assertRaises(ReferenceGenomeNotRecognized):
            ReferenceGenome.get_instance("unknown")


class TestSampleFastq(unittest.TestCase):

    def setUp(self):
        self.sample_name = "AB-0123_Sample-456"
        self.sample_number = "S2"
        self.lane_number = "3"
        self.read_number = "2"
        self.file_extension = ".fastq.gz"

    def get_file_name(self):
        return "_".join([
            self.sample_name,
            self.sample_number,
            "L00{}".format(self.lane_number),
            "R{}".format(self.read_number),
            "001{}".format(self.file_extension)])

    def test_split_filename(self):
        fastq_file_name = self.get_file_name()
        fastq_file_path = os.path.join("/path", "to")
        sample_fastq = SampleFastq(os.path.join(fastq_file_path, fastq_file_name))
        self.assertEqual(self.sample_name, sample_fastq.sample_name)
        self.assertEqual(self.sample_number, sample_fastq.sample_number)
        self.assertEqual(self.lane_number, sample_fastq.lane_number)
        self.assertEqual(self.read_number, sample_fastq.read_number)
        self.assertEqual(fastq_file_path, sample_fastq.dirname)
        self.assertEqual(fastq_file_name, sample_fastq.filename)


class TestRunfolder(unittest.TestCase):

    def setUp(self):
        self.run_date = "180406"
        self.instrument_id = "ST-01234"
        self.run_number = "0123"
        self.flowcell_position = "A"
        self.flowcell_id = "ABC123CXY"

    def get_runfolder_name(self):
        return "_".join([
            self.run_date,
            self.instrument_id,
            self.run_number,
            "{}{}".format(self.flowcell_position, self.flowcell_id)])

    def test_split_runfolder_name(self):
        runfolder_path = os.path.join("/path", "to")
        runfolder_name = self.get_runfolder_name()
        runfolder = Runfolder(os.path.join(runfolder_path, runfolder_name))
        self.assertEqual(runfolder_path, runfolder.dirname)
        self.assertEqual(runfolder_name, runfolder.runfolder_name)
        self.assertEqual(self.run_date, runfolder.run_date)
        self.assertEqual(self.instrument_id, runfolder.instrument_id)
        self.assertEqual(self.run_number, runfolder.run_number)
        self.assertEqual(self.flowcell_position, runfolder.flowcell_position)
        self.assertEqual(self.flowcell_id, runfolder.flowcell_id)

