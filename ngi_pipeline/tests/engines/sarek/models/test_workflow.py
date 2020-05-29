import mock
import os
import unittest

from ngi_pipeline.engines.sarek.exceptions import ParserException
from ngi_pipeline.engines.sarek.models.sample import SarekAnalysisSample
from ngi_pipeline.engines.sarek.models.workflow import SarekWorkflowStep, SarekMainStep


class TestSarekWorkflowStep(unittest.TestCase):

    CONFIG = {
        "profile": "this-is-a-profile",
        "tools": ["haplotypecaller", "manta"],
        "nf_path": os.path.join("/this", "is", "the", "nextflow", "path"),
        "sarek_path": os.path.join("/this", "is", "the", "path", "to", "sarek"),
        "genome": "this-is-the-genome"}

    def setUp(self):
        self.sarek_args = self.CONFIG.copy()
        self.sarek_workflow_step = SarekWorkflowStep(**self.sarek_args)

    def test__append_argument(self):
        base_string = "this-is-the-base-string"

        # test a non-existing attribute
        attribute = "this-attribute-does-not-exist"
        self.assertEqual(base_string, self.sarek_workflow_step._append_argument(base_string, attribute))

        # test a None attribute
        attribute = "none_attribute"
        setattr(self.sarek_workflow_step, attribute, None)
        self.assertEqual(base_string, self.sarek_workflow_step._append_argument(base_string, attribute))

        # test a list attribute
        attribute = "list_attribute"
        value = ["this", "is", "a", "list"]
        self.sarek_workflow_step.sarek_args[attribute] = value
        expected_result = "{0} --{1} ${{{1}}}".format(base_string, attribute)
        self.assertEqual(expected_result, self.sarek_workflow_step._append_argument(base_string, attribute))

        # test a string attribute
        attribute = "string_attribute"
        value = "this-is-a-string"
        self.sarek_workflow_step.sarek_args[attribute] = value
        expected_result = "{0} --{1} ${{{1}}}".format(base_string, attribute)
        self.assertEqual(expected_result, self.sarek_workflow_step._append_argument(base_string, attribute))

        # test a custom hyphen string
        hyphen = "xyz"
        expected_result = "{0} {2}{1} ${{{1}}}".format(base_string, attribute, hyphen)
        self.assertEqual(
            expected_result, self.sarek_workflow_step._append_argument(
                base_string, attribute, hyphen=hyphen))

    def test_sarek_step(self):
        with self.assertRaises(NotImplementedError):
            self.sarek_workflow_step.sarek_step()

    def test_command_line(self):
        sarek_step = "sarek_step"
        sarek_workflow_step = SarekWorkflowStep()
        with mock.patch.object(sarek_workflow_step, "sarek_step", return_value=sarek_step) as sarek_step_mock:
            expected_command_line = "sarek"
            self.assertEqual(expected_command_line, sarek_workflow_step.command_line())

    def test_command_line_args(self):
        sarek_step = "sarek_step"
        valid_tools = self.sarek_args["tools"]
        with mock.patch.object(self.sarek_workflow_step, "sarek_step", return_value=sarek_step):
            observed_command_line = self.sarek_workflow_step.command_line()
            self.assertIn("--tools", observed_command_line)
            self.assertIn(",".join(valid_tools), observed_command_line)
            for key, value in self.sarek_args.items():
                if key not in ["tools", "sarek_path", "nf_path"]:
                    self.assertIn("-{} {}".format(key, value), observed_command_line)


class TestSarekMainStep(unittest.TestCase):

    def test_report_files(self):
        analysis_sample = mock.Mock(spec=SarekAnalysisSample)
        analysis_sample.sampleid = "this-is-a-sample-id"
        analysis_sample.sample_analysis_path.return_value = "this-is-a-path"
        with mock.patch("os.listdir") as list_mock:
            file_list = ["file1", "file2.extension", "file3.metrics", "file4.metrics"]
            for listed_files in [file_list[0:2], file_list]:
                list_mock.return_value = listed_files
                with self.assertRaises(ParserException):
                    SarekMainStep.report_files(analysis_sample)
            list_mock.return_value = file_list[0:3]
            self.assertEqual(
                file_list[2], os.path.basename(SarekMainStep.report_files(analysis_sample)[1][1]))
