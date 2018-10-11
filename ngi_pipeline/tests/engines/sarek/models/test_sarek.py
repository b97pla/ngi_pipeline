import mock
import os
import shutil
import tempfile
import unittest

from ngi_pipeline.engines.sarek.database import CharonConnector
from ngi_pipeline.engines.sarek.exceptions import BestPracticeAnalysisNotRecognized, SampleNotValidForAnalysisError
from ngi_pipeline.engines.sarek.models.resources import SampleFastq
from ngi_pipeline.engines.sarek.models.sample import SarekAnalysisSample
from ngi_pipeline.engines.sarek.models.sarek import SarekAnalysis, SarekGermlineAnalysis
from ngi_pipeline.engines.sarek.parsers import ParserIntegrator
from ngi_pipeline.log.loggers import minimal_logger
from ngi_pipeline.tests.engines.sarek.test_launchers import TestLaunchers


@mock.patch("ngi_pipeline.engines.sarek.models.sarek.ReferenceGenome", autospec=True)
@mock.patch("ngi_pipeline.engines.sarek.database.CharonConnector", autospec=True)
class TestSarekAnalysis(unittest.TestCase):

    CONFIG = {}

    def setUp(self):
        self.log = minimal_logger(__name__, to_file=False, debug=True)
        self.config = TestSarekAnalysis.CONFIG
        self.analysis_obj = TestLaunchers.get_NGIAnalysis(log=self.log)

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

    def test_status_should_be_started(self, charon_connector_mock, reference_genome_mock):

        def _analysis_status_wrapper(status, *args):
            return CharonConnector(
                self.config, self.log, charon_session="charon-session").analysis_status_from_process_status(status)

        def _alignment_status_wrapper(status, *args):
            return CharonConnector(
                self.config, self.log, charon_session="charon-session").alignment_status_from_analysis_status(status)

        charon_connector = charon_connector_mock.return_value
        charon_connector.analysis_status_from_process_status.side_effect = _analysis_status_wrapper
        charon_connector.alignment_status_from_analysis_status.side_effect = _alignment_status_wrapper
        sarek_analysis = SarekAnalysis(
            reference_genome_mock.return_value, self.config, self.log, charon_connector=charon_connector)

        for restart_mode in True, False:
            self.assertEqual(
                restart_mode,
                sarek_analysis.status_should_be_started(
                    "FAILED",
                    restart_failed_jobs=restart_mode))

        for restart_mode in True, False:
            self.assertEqual(
                restart_mode,
                sarek_analysis.status_should_be_started(
                    "ANALYZED",
                    restart_finished_jobs=restart_mode))

        for restart_mode in True, False:
            self.assertEqual(
                restart_mode,
                sarek_analysis.status_should_be_started(
                    "UNDER_ANALYSIS",
                    restart_running_jobs=restart_mode))

        self.assertTrue(
            sarek_analysis.status_should_be_started("TO_BE_ANALYZED"))

    def test_configure_analysis(self, charon_connector_mock, reference_genome_mock):
        config = {
            "profile": "this-is-a-profile",
            "tools": ["tool-A", "tool-B"]}

        sarek_analysis = SarekAnalysis(
            reference_genome_mock.return_value,
            {"sarek": config},
            self.log,
            charon_connector=charon_connector_mock.return_value)

        for key in config.keys():
            self.assertEqual(config[key], sarek_analysis.sarek_config[key])
        self.assertEqual(SarekAnalysis.DEFAULT_CONFIG["nf_path"], sarek_analysis.sarek_config["nf_path"])
        self.assertEqual(SarekAnalysis.DEFAULT_CONFIG["sarek_path"], sarek_analysis.sarek_config["sarek_path"])

        config = {
            "config": "this-is-a-site-specific-config-file",
            "profile": "this-is-another-profile",
            "nf_path": "this-is-the-path-to-nextflow",
            "sarek_path": "this-is-the-path-to-sarek"}
        sarek_config = sarek_analysis.configure_analysis(config={"sarek": config})
        for key in config.keys():
            self.assertEqual(config[key], sarek_config[key])
        self.assertEqual(SarekAnalysis.DEFAULT_CONFIG["tools"], sarek_config["tools"])

    def test_command_line(self, charon_connector_mock, reference_genome_mock):

        sarek_analysis = SarekAnalysis(
            reference_genome_mock.return_value,
            self.config,
            self.log,
            charon_connector=charon_connector_mock.return_value)

        with self.assertRaises(NotImplementedError):
            sarek_analysis.command_line(None)

    def test_generate_tsv_file_contents(self, charon_connector_mock, reference_genome_mock):

        sarek_analysis = SarekAnalysis(
            reference_genome_mock.return_value,
            self.config,
            self.log,
            charon_connector=charon_connector_mock.return_value)

        with self.assertRaises(NotImplementedError):
            sarek_analysis.generate_tsv_file_contents(None)

    def test_create_tsv_file(self, charon_connector_mock, reference_genome_mock):
        sarek_analysis = SarekAnalysis(
            reference_genome_mock.return_value,
            self.config,
            self.log,
            charon_connector=charon_connector_mock.return_value)
        sample_obj = self.analysis_obj.project.samples.values()[0]
        analysis_sample = SarekAnalysisSample(self.analysis_obj.project, sample_obj, sarek_analysis)
        with mock.patch.object(sarek_analysis, "generate_tsv_file_contents") as tsv_mock, \
                mock.patch.object(analysis_sample, "sample_analysis_tsv_file") as tsv_path_mock:
            tsv_mock.return_value = []
            with self.assertRaises(SampleNotValidForAnalysisError) as e:
                sarek_analysis.create_tsv_file(analysis_sample)
            tsv_mock.return_value = [["this", "is"], ["some", "content"]]
            tempdir = tempfile.mkdtemp(prefix="test_create_tsv_file_")
            tsv_path_mock.return_value = os.path.join(tempdir, "tsv_parent_folder", "tsv_file.tsv")
            sarek_analysis.create_tsv_file(analysis_sample)
            self.assertTrue(os.path.exists(tsv_path_mock.return_value))
            shutil.rmtree(tempdir)


@mock.patch("ngi_pipeline.engines.sarek.models.resources.ReferenceGenome", autospec=True)
@mock.patch("ngi_pipeline.engines.sarek.database.CharonConnector", autospec=True)
@mock.patch("ngi_pipeline.engines.sarek.database.TrackingConnector", autospec=True)
@mock.patch("ngi_pipeline.engines.sarek.process.ProcessConnector", autospec=True)
class TestSarekGermlineAnalysis(unittest.TestCase):

    CONFIG = {
        "profile": "this-is-a-profile",
        "tools": ["haplotypecaller", "manta"],
        "nf_path": os.path.join("/this", "is", "the", "nextflow", "path"),
        "sarek_path": os.path.join("/this", "is", "the", "path", "to", "sarek"),
        "genome": "this-is-the-genome"}

    def setUp(self):
        self.log = minimal_logger(__name__, to_file=False, debug=True)
        self.config = TestSarekGermlineAnalysis.CONFIG
        self.analysis_obj = TestLaunchers.get_NGIAnalysis(log=self.log)

    def get_instance(
            self, process_connector_mock, tracking_connector_mock, charon_connector_mock, reference_genome_mock):
        reference_genome = reference_genome_mock.return_value
        reference_genome.__str__.return_value = self.config["genome"]
        return SarekGermlineAnalysis(
            reference_genome,
            {"sarek": self.config},
            self.log,
            charon_connector=charon_connector_mock,
            tracking_connector=tracking_connector_mock,
            process_connector=process_connector_mock)

    def test_command_line(
            self, process_connector_mock, tracking_connector_mock, charon_connector_mock, reference_genome_mock):
        sarek_analysis = self.get_instance(
            process_connector_mock, tracking_connector_mock, charon_connector_mock, reference_genome_mock)
        sample_obj = self.analysis_obj.project.samples.values()[0]
        analysis_sample = SarekAnalysisSample(self.analysis_obj.project, sample_obj, sarek_analysis)
        self.config["outDir"] = analysis_sample.sample_analysis_path()
        self.config["sample"] = analysis_sample.sample_analysis_tsv_file()
        observed_cmd = sarek_analysis.command_line(analysis_sample)

        self.assertIn("-profile {}".format(self.config["profile"]), observed_cmd)
        self.assertIn("--tools {}".format(",".join(self.config["tools"])), observed_cmd)
        self.assertTrue(observed_cmd.startswith("{} run {}".format(self.config["nf_path"], self.config["sarek_path"])))
        for key in filter(lambda k: k not in ["profile", "nf_path", "sarek_path", "tools"], self.config.keys()):
            self.assertIn("--{} {}".format(key, self.config[key]), observed_cmd)

    def test_generate_tsv_file_contents(
            self, process_connector_mock, tracking_connector_mock, charon_connector_mock, reference_genome_mock):
        sarek_analysis = self.get_instance(
            process_connector_mock, tracking_connector_mock, charon_connector_mock, reference_genome_mock)
        sample_obj = self.analysis_obj.project.samples.values()[0]
        analysis_sample = SarekAnalysisSample(self.analysis_obj.project, sample_obj, sarek_analysis)
        sample_data_path = analysis_sample.sample_data_path()
        with mock.patch.object(
                sarek_analysis, "libprep_should_be_started") as libprep_mock, \
                mock.patch.object(
                    sarek_analysis, "seqrun_should_be_started") as seqrun_mock:
            for start_mode in True, False:
                libprep_mock.return_value = start_mode
                seqrun_mock.return_value = not start_mode
                self.assertListEqual(
                    [],
                    sarek_analysis.generate_tsv_file_contents(analysis_sample))
            libprep_mock.return_value = True
            seqrun_mock.return_value = True
            expected_tsv = []
            for libprepid in map(lambda lp: '{}-{}'.format(sample_obj.name, lp), ['libprep1', 'libprep2']):
                for sample_index in ['1', '2']:
                    tsv_row = [sample_obj.name, 'ZZ', 0, sample_obj.name, 'ABC00{0}CXY.1.S{0}'.format(sample_index)]
                    tsv_row.extend([
                        os.path.join(
                            sample_data_path,
                            libprepid,
                            '180411_ST-0123_001{0}_AABC00{0}CXY'.format(sample_index),
                            '{}_S{}_L001_R{}_001.fastq.gz'.format(
                                libprepid, sample_index, read_num)) for read_num in ['1', '2']])
                    expected_tsv.append(tsv_row)
            observed_tsv = sarek_analysis.generate_tsv_file_contents(analysis_sample)
            self.assertListEqual(sorted(expected_tsv), sorted(observed_tsv))

    def test__sample_fastq_file_pair_sorting(self, *args):
        sample_numbers = {
            "S03": "sample_A",
            "S1": "sample_A",
            "S12": "sample_A",
            "S2": "sample_A",
            "S3": "sample_B"}
        lane_numbers = ["2", "4", "5"]
        sample_fastq_files = {
            "1": [],
            "2": []}
        for sample_number in sorted(sample_numbers.keys()):
            sample_name = sample_numbers[sample_number]
            for lane_number in lane_numbers:
                for read_number in sorted(sample_fastq_files.keys()):
                    sample_fastq_files[read_number].append(
                        SampleFastq("{}_{}_L00{}_R{}_001.fastq.gz".format(
                            sample_name, sample_number, lane_number, read_number)))
        n = 0
        for fastq_pair in SarekGermlineAnalysis.sample_fastq_file_pair(
                sample_fastq_files["1"] + sample_fastq_files["2"]):
            for i in [0, 1]:
                self.assertEqual(sample_fastq_files[sorted(sample_fastq_files.keys())[i]][n], fastq_pair[i])
            n += 1

    def test__sample_fastq_file_pair_single_read(self, *args):
        sample_fastq_files = map(SampleFastq, [
            "sample_A_S1_L001_R1_001.fastq.gz",
            "sample_A_S1_L002_R1_001.fastq.gz",
            "sample_A_S2_L001_R1_001.fastq.gz",
            "sample_A_S2_L005_R1_001.fastq.gz",
            "sample_B_S3_L002_R1_001.fastq.gz"
        ])
        n = 0
        for fastq_pair in SarekGermlineAnalysis.sample_fastq_file_pair(sample_fastq_files):
            self.assertEqual(1, len(fastq_pair))
            self.assertEqual(sample_fastq_files[n], fastq_pair[0])
            n += 1

    def test_analyze_sample(
            self, process_connector_mock, tracking_connector_mock, charon_connector_mock, reference_genome_mock):
        sarek_analysis = self.get_instance(
            process_connector_mock, tracking_connector_mock, charon_connector_mock, reference_genome_mock)
        with mock.patch.object(
                sarek_analysis, "create_tsv_file", return_value=os.path.join("/path", "to", "tsv", "file")) as tsv_mock:
            for sample_obj in self.analysis_obj.project:
                sarek_analysis.analyze_sample(sample_obj, self.analysis_obj)

    def test_fastq_files_from_tsv_file(self, *args):
        fh, fake_tsv_file = tempfile.mkstemp(prefix="test_fastq_files_from_tsv_")
        with mock.patch("ngi_pipeline.engines.sarek.models.sarek.csv.reader", autospec=True) as csv_mock:
            expected_fastq_files = [
                os.path.join("/path", "to", "fastq", "file_{}.fastq.gz".format(i+1)) for i in range(4)]
            mocked_tsv_data = (
                [1, 2, 3, 4, 5, expected_fastq_files[0], expected_fastq_files[1]],
                [6, 7, 8, 9, 10, expected_fastq_files[2], expected_fastq_files[3]])
            csv_mock.return_value = mocked_tsv_data
            observed_fastq_files = SarekGermlineAnalysis.fastq_files_from_tsv_file(fake_tsv_file)
            self.assertListEqual(sorted(expected_fastq_files), sorted(observed_fastq_files))
        os.unlink(fake_tsv_file)

    def test_collect_analysis_metrics(
            self, process_connector_mock, tracking_connector_mock, charon_connector_mock, reference_genome_mock):
        sarek_analysis = self.get_instance(
            process_connector_mock, tracking_connector_mock, charon_connector_mock, reference_genome_mock)
        expected_metrics = {
            "percent_duplication": [25.3],
            "autosomal_coverage": [35.8],
            "total_reads": [123456789]
        }

        def _serve_metric(metric):
            return expected_metrics[metric[4:]]

        with mock.patch.object(
                sarek_analysis, "processing_steps") as processing_steps, \
                mock.patch(
                    "ngi_pipeline.engines.sarek.models.sarek.ParserIntegrator", autospec=ParserIntegrator) as \
                        parser_mock:
            processing_steps.return_value = []
            parser_instance = parser_mock.return_value
            query_mock = parser_instance.query_parsers
            query_mock.side_effect = _serve_metric
            self.assertDictEqual(expected_metrics, sarek_analysis.collect_analysis_metrics("this-is-a-sample-analysis"))
