from ngi_pipeline.database.classes import CharonError, CharonSession
from ngi_pipeline.engines.sarek.exceptions import BestPracticeAnalysisNotSpecifiedError, \
    SampleAnalysisStatusNotFoundError, SampleAnalysisStatusNotSetError


class CharonConnector:

    def __init__(self, config, log, charon_session=None):
        self.config = config
        self.log = log
        self.charon_session = charon_session or CharonSession(config=self.config)

    def best_practice_analysis(self, projectid):
        try:
            return self.charon_session.project_get(projectid)["best_practice_analysis"]
        except (KeyError, CharonError) as e:
            best_practice_analysis_exception = BestPracticeAnalysisNotSpecifiedError(projectid, reason=e)
            self.log.error(best_practice_analysis_exception)
            raise best_practice_analysis_exception

    def sample_analysis_status(self, projectid, sampleid):
        try:
            return self.charon_session.sample_get(projectid, sampleid)["analysis_status"]
        except (KeyError, CharonError) as e:
            analysis_status_exception = SampleAnalysisStatusNotFoundError(projectid, sampleid, reason=e)
            self.log.error(analysis_status_exception)
            raise analysis_status_exception

    def set_sample_analysis_status(self, projectid, sampleid, status):
        try:
            return self.charon_session.sample_update(projectid, sampleid, analysis_status=status)
        except CharonError as e:
            analysis_status_exception = SampleAnalysisStatusNotSetError(projectid, sampleid, status, reason=e)
            self.log.error(analysis_status_exception)
            raise analysis_status_exception
    