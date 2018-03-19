import contextlib

from ngi_pipeline.database.classes import CharonError, CharonSession
from ngi_pipeline.engines.sarek.process import ProcessRunning, ProcessExitStatusSuccessful, \
    ProcessExitStatusFailed, ProcessExitStatusUnknown, ProcessConnector, SlurmConnector
from ngi_pipeline.engines.piper_ngi.database import get_db_session, SampleAnalysis
from ngi_pipeline.engines.sarek.exceptions import BestPracticeAnalysisNotSpecifiedError, SampleLookupError, \
    AnalysisStatusForProcessStatusNotFoundError, SampleAnalysisStatusNotFoundError, SampleAnalysisStatusNotSetError, \
    SampleAlignmentStatusNotSetError, AlignmentStatusForAnalysisStatusNotFoundError


class CharonConnector:

    _ANALYSIS_STATUS_FROM_PROCESS_STATUS = {
        ProcessRunning: "UNDER_ANALYSIS",
        ProcessExitStatusSuccessful: "ANALYZED",
        ProcessExitStatusFailed: "FAILED",
        ProcessExitStatusUnknown: "FAILED"
    }

    _ALIGNMENT_STATUS_FROM_ANALYSIS_STATUS = {
        "TO_ANALYZE": "NOT_RUNNING",
        "UNDER_ANALYSIS": "RUNNING",
        "ANALYZED": "DONE",
        "FAILED": "FAILED"
    }

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

    def set_sample_analysis_status(
            self, projectid, sampleid, status, recurse=False, restrict_to_libpreps=None, restrict_to_seqruns=None):
        try:
            if recurse:
                # update the status for all seqruns
                for libprep in self.sample_libpreps(
                        projectid, sampleid, restrict_to=restrict_to_libpreps):
                    libprepid = libprep["libprepid"]
                    for seqrun in self.libprep_seqruns(
                            projectid, sampleid, libprepid,
                            restrict_to=restrict_to_seqruns.get(libprepid) if restrict_to_seqruns else None):
                        self.set_seqrun_alignment_status(
                            projectid,
                            sampleid,
                            libprepid,
                            seqrun["seqrunid"],
                            self.alignment_status_from_analysis_status(status))

            return self.charon_session.sample_update(projectid, sampleid, analysis_status=status)
        except CharonError as e:
            analysis_status_exception = SampleAnalysisStatusNotSetError(projectid, sampleid, status, reason=e)
            self.log.error(analysis_status_exception)
            raise analysis_status_exception

    def sample_libpreps(self, projectid, sampleid, restrict_to=None):
        try:
            return filter(
                lambda x: restrict_to is None or x["libprepid"] in restrict_to,
                self.charon_session.sample_get_libpreps(projectid, sampleid)["libpreps"])
        except (KeyError, CharonError) as e:
            sample_libpreps_exception = SampleLookupError(projectid, sampleid, reason=e)
            self.log.error(sample_libpreps_exception)
            raise sample_libpreps_exception

    def libprep_seqruns(self, projectid, sampleid, libprepid, restrict_to=None):
        try:
            return filter(
                lambda x: restrict_to is None or x["seqrunid"] in restrict_to,
                self.charon_session.libprep_get_seqruns(projectid, sampleid, libprepid)["seqruns"])
        except (KeyError, CharonError) as e:
            libprep_seqruns_exception = SampleLookupError(projectid, sampleid, reason=e)
            self.log.error(libprep_seqruns_exception)
            raise libprep_seqruns_exception

    def set_seqrun_alignment_status(self, projectid, sampleid, libprepid, seqrunid, status):
        try:
            return self.charon_session.seqrun_update(projectid, sampleid, libprepid, seqrunid, alignment_status=status)
        except CharonError as e:
            alignment_status_exception = SampleAlignmentStatusNotSetError(
                projectid, sampleid, libprepid, seqrunid, status, reason=e)
            self.log.error(alignment_status_exception)
            raise alignment_status_exception

    def analysis_status_from_process_status(self, process_status):
        try:
            return self._ANALYSIS_STATUS_FROM_PROCESS_STATUS[process_status]
        except KeyError as ke:
            charon_status_error = AnalysisStatusForProcessStatusNotFoundError(process_status)
            self.log.error(charon_status_error)
            raise charon_status_error

    def alignment_status_from_analysis_status(self, analysis_status):
        try:
            return self._ALIGNMENT_STATUS_FROM_ANALYSIS_STATUS[analysis_status]
        except KeyError as ke:
            status_error = AlignmentStatusForAnalysisStatusNotFoundError(analysis_status)
            self.log.error(status_error)
            raise status_error


class TrackingConnector:

    PIDFIELD_FROM_PROCESS_CONNECTOR_TYPE = {
        ProcessConnector: "process_id",
        SlurmConnector: "slurm_job_id"
    }

    def __init__(self, config, log, tracking_session=None):
        self.config = config
        self.log = log
        self.tracking_session = tracking_session

    class _SampleAnalysis(SampleAnalysis):
        pass

    @staticmethod
    def pidfield_from_process_connector_type(process_connector_type):
        return TrackingConnector.PIDFIELD_FROM_PROCESS_CONNECTOR_TYPE[process_connector_type]

    @contextlib.contextmanager
    def db_session(self):
        if self.tracking_session is not None:
            yield self.tracking_session
        else:
            with get_db_session(config=self.config) as db_session:
                self.tracking_session = db_session
                yield self.tracking_session

    def record_process_sample(
            self, projectid, sampleid, project_base_path, analysis_type, engine, pid, pidfield):

        db_obj = self._SampleAnalysis(
            project_id=projectid,
            project_name=projectid,
            sample_id=sampleid,
            project_base_path=project_base_path,
            workflow=analysis_type,
            engine=engine,
            **{pidfield: pid})

        with self.db_session() as db_session:
            db_session.add(db_obj)
            db_session.commit()

    def remove_analysis(self, analysis):
        with self.db_session() as db_session:
            db_session.delete(analysis)
            db_session.commit()

    def tracked_analyses(self):
        with self.db_session() as db_session:
            for analysis in db_session.query(self._SampleAnalysis)\
                    .filter(self._SampleAnalysis.engine == "sarek")\
                    .all():
                yield analysis
