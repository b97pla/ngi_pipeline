

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


class SampleNotValidForAnalysisError(Exception):

    def __init__(self, projectid, sampleid, reason=None):
        super(SampleNotValidForAnalysisError, self).__init__(
            "requirements for starting analysis not fulfilled for sample '{}' in project '{}'{}".format(
                sampleid, projectid, ": {}".format(reason.__repr__()) if reason else ""))
        self.projectid = projectid
        self.sampleid = sampleid
        self.reason = reason


class SampleLookupError(Exception):

    def __init__(self, projectid, sampleid, reason=None):
        super(SampleLookupError, self).__init__(
            "error fetching property for sample '{}' in project '{}'{}".format(
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


class SampleAlignmentStatusNotSetError(Exception):

    def __init__(self, projectid, sampleid, libprepid, seqrunid, status, reason=None):
        super(SampleAlignmentStatusNotSetError, self).__init__(
            "seqrun alignment status '{}' could not be set for seqrun '{}' on libprep '{}' for sample '{}' in project"
            " '{}'{}".format(status, seqrunid, libprepid, sampleid, projectid, ": {}".format(
                reason.__repr__()) if reason else ""))
        self.projectid = projectid
        self.sampleid = sampleid
        self.libprepid = libprepid
        self.seqrunid = seqrunid
        self.status = status
        self.reason = reason


class SlurmStatusNotRecognizedError(Exception):

    def __init__(self, jobid, job_state):
        super(SlurmStatusNotRecognizedError, self).__init__(
            "SLURM job state '{}' for job id {} is not recognized".format(job_state, jobid))
        self.jobid = jobid
        self.job_state = job_state


class AnalysisStatusForProcessStatusNotFoundError(Exception):

    def __init__(self, process_status):
        super(AnalysisStatusForProcessStatusNotFoundError, self).__init__(
            "Charon analysis status corresponding to process status {} not found".format(process_status.__name__))
        self.process_status = process_status


class AlignmentStatusForAnalysisStatusNotFoundError(Exception):

    def __init__(self, analysis_status):
        super(AlignmentStatusForAnalysisStatusNotFoundError, self).__init__(
            "Charon alignment status corresponding to analysis status {} not found".format(analysis_status))
        self.analysis_status = analysis_status


class ParserException(Exception):

    def __init__(self, instance, message):
        super(ParserException, self).__init__(
            "{} raised an exception: {}".format(type(instance).__name__, message)
        )
        self.instance = instance
