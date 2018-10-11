import mock
import os
import unittest

from ngi_pipeline.engines.sarek.models.workflow import SarekWorkflowStep


class TestSarekWorkflowStep(unittest.TestCase):

    CONFIG = {
        "profile": "this-is-a-profile",
        "tools": ["haplotypecaller", "manta"],
        "nf_path": os.path.join("/this", "is", "the", "nextflow", "path"),
        "sarek_path": os.path.join("/this", "is", "the", "path", "to", "sarek"),
        "genome": "this-is-the-genome"}

    def setUp(self):
        self.sarek_args = self.CONFIG.copy()
        self.sarek_workflow_step = SarekWorkflowStep(
            self.sarek_args["nf_path"], self.sarek_args["sarek_path"], **self.sarek_args)

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
        sarek_workflow_step = SarekWorkflowStep(self.sarek_args["nf_path"], self.sarek_args["sarek_path"])
        with mock.patch.object(sarek_workflow_step, "sarek_step", return_value=sarek_step) as sarek_step_mock:
            expected_command_line = "{} run {}".format(
                self.sarek_args["nf_path"],
                os.path.join(self.sarek_args["sarek_path"], sarek_step))
            self.assertEqual(expected_command_line, sarek_workflow_step.command_line())
            sarek_step_mock.assert_called_once()

    def test_command_line_args(self):
        sarek_step = "sarek_step"
        valid_tools = self.sarek_args["tools"]
        with mock.patch.object(self.sarek_workflow_step, "sarek_step", return_value=sarek_step):
            observed_command_line = self.sarek_workflow_step.command_line()
            self.assertNotIn("--tools", observed_command_line)
            for key, value in self.sarek_args.items():
                if key not in ["tools", "sarek_path", "nf_path"]:
                    self.assertIn("-{} {}".format(key, value), observed_command_line)

    def test_valid_tools(self):
        self.sarek_workflow_step.available_tools = []
        self.assertListEqual([], self.sarek_workflow_step.valid_tools(["A", "B", "C", "D"]))
        self.sarek_workflow_step.available_tools = ["B", "D"]
        self.assertListEqual(["B", "D"], self.sarek_workflow_step.valid_tools(["A", "B", "C", "D"]))

