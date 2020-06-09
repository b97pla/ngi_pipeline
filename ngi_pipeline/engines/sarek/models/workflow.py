import os
from string import Template

from ngi_pipeline.engines.sarek.exceptions import ParserException
from ngi_pipeline.engines.sarek.parsers import QualiMapParser, PicardMarkDuplicatesParser


class SarekWorkflowStep(object):
    """
    The SarekWorkflowStep class represents an analysis step in the Sarek workflow. Primarily, it provides a method for
    creating the step-specific command line.
    """

    available_tools = []

    def __init__(self, **sarek_args):
        """
        Create a SarekWorkflowStep instance according to the passed parameters.

        :param sarek_args: additional Sarek parameters to be included on the command line
        """
        # use a separate variable for the command to invoke sarek. This defaults to just "sarek" which is a defined
        # alias on Irma but it can be overridden through the pipeline config option "sarek_cmd"
        self.sarek_cmd = sarek_args.get("sarek_cmd", "sarek")
        # create a dict with parameters based on the passed key=value arguments
        self.sarek_args = {k: v for k, v in sarek_args.items() if k not in ["sarek_cmd"]}
        # expand any parameters passed as list items into a ","-separated string
        self.sarek_args = {k: v if type(v) is not list else ",".join(v) for k, v in self.sarek_args.items()}

    def _append_argument(self, base_string, name, hyphen="--"):
        """
        Append an argument with a placeholder for the value to the supplied string in a format suitable for the
        string.Template constructor. If no value exists for the argument name among this workflow step's config
        parameters, the supplied string is returned untouched.

        Example: step._append_argument("echo", "hello", "") should return "echo hello ${hello}", provided the step
        instance has a "hello" key in the step.sarek_args dict.

        :param base_string: the string to append an argument to
        :param name: the argument name to add a placeholder for
        :param hyphen: the hyphen style to prefix the argument name with (default "--")
        :return: the supplied string with an appended argument name and placeholder
        """
        # NOTE: a numeric value of 0 will be excluded (as will a boolean value of False)!
        if not self.sarek_args.get(name):
            return base_string
        return "{0} {2}{1} ${{{1}}}".format(base_string, name, hyphen)

    def command_line(self):
        """
        Generate the command line for launching this analysis workflow step. The command line will be built using the
        Sarek arguments passed to the step's constructor and returned as a string.

        :return: the command line for the workflow step as a string
        """
        single_hyphen_args = ["config", "profile", "resume"]
        template_string = "${sarek_cmd}"
        for argument_name in single_hyphen_args:
            template_string = self._append_argument(template_string, argument_name, hyphen="-")
        for argument_name in filter(lambda n: n not in single_hyphen_args, self.sarek_args.keys()):
            template_string = self._append_argument(template_string, argument_name, hyphen="--")
        command_line = Template(template_string).substitute(
            sarek_cmd=self.sarek_cmd,
            **self.sarek_args)
        return command_line

    def sarek_step(self):
        raise NotImplementedError("The Sarek workflow step definition for {} has not been defined".format(type(self)))

    @classmethod
    def report_files(cls, analysis_sample):
        return []


class SarekMainStep(SarekWorkflowStep):

    def sarek_step(self):
        return "main.nf"

    @classmethod
    def report_files(cls, analysis_sample):
        """
        Get a list of the report files resulting from this processing step and the associated parsers.

        :param analysis_sample: the SarekAnalysisSample that was analyzed
        :return: a list of tuples where the first element is a parser class instance and the second is the path to the
        result file that the parser instance should parse
        """
        report_dir = os.path.join(
            analysis_sample.sample_analysis_results_dir(),
            "Reports",
            analysis_sample.sampleid)
        # MarkDuplicates output files may be named differently depending on if the pipeline was started with a single
        # fastq file pair or multiple file pairs
        markdups_dir = os.path.join(report_dir, "MarkDuplicates")
        metric_files = filter(lambda f: f.endswith(".metrics"), os.listdir(markdups_dir))
        if not metric_files:
            raise ParserException(cls, "no metrics file for MarkDuplicates found for sample {} in {}".format(
                analysis_sample.sampleid, markdups_dir))
        markdups_metrics_file = metric_files.pop()
        if metric_files:
            raise ParserException(cls, "multiple metrics files for MarkDuplicates found for sample {} in {}".format(
                analysis_sample.sampleid, markdups_dir))
        return [
            [
                QualiMapParser,
                os.path.join(report_dir, "bamQC", "{}.recal".format(analysis_sample.sampleid), "genome_results.txt")],
            [
                PicardMarkDuplicatesParser,
                os.path.join(markdups_dir, markdups_metrics_file)]]
