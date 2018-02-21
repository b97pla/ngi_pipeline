

class SarekException(Exception):
    pass


class BestPracticeAnalysisNotSpecifiedError(SarekException):

    def __init__(self, projectid, reason=None):
        super(BestPracticeAnalysisNotSpecifiedError, self).__init__(
            "best-practice-analysis could not be found for project '{}'{}".format(
                projectid,
                ": {}".format(reason.__repr__()) if reason else ""))
        self.projectid = projectid
        self.reason = reason


class BestPracticeAnalysisNotRecognized(SarekException):

    def __init__(self, best_practice_analysis):
        super(BestPracticeAnalysisNotRecognized, self).__init__(
            "best-practice-analysis '{}' not recognized".format(best_practice_analysis))
        self.best_practice_analysis= best_practice_analysis


class ReferenceGenomeNotRecognized(Exception):

    def __init__(self, genome_name):
        super(ReferenceGenomeNotRecognized, self).__init__(
            "reference genome '{}' not recognized".format(genome_name))
        self.genome_name = genome_name


class SampleAnalysisStatusNotFoundError(Exception):

    def __init__(self, projectid, sampleid, reason=None):
        super(SampleAnalysisStatusNotFoundError, self).__init__(
            "sample analysis status not found for sample '{}' in project '{}'{}".format(
                sampleid, projectid, ": {}".format(reason.__repr__()) if reason else ""))
        self.projectid = projectid
        self.sampleid = sampleid
        self.reason = reason


class SampleAnalysisStatusNotSetError(Exception):

    def __init__(self, projectid, sampleid, status, reason=None):
        super(SampleAnalysisStatusNotSetError, self).__init__(
            "sample analysis status '{}' could not be set for sample '{}' in project '{}'{}".format(
                status, sampleid, projectid, ": {}".format(reason.__repr__()) if reason else ""))
        self.projectid = projectid
        self.sampleid = sampleid
        self.status = status
        self.reason = reason
