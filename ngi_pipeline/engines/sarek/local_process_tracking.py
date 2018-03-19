from ngi_pipeline.engines.sarek.database import CharonConnector, TrackingConnector
from ngi_pipeline.engines.sarek.models import SarekAnalysis
from ngi_pipeline.engines.sarek.process import JobStatus, ProcessStatus, ProcessRunning
from ngi_pipeline.log.loggers import minimal_logger
from ngi_pipeline.utils.classes import with_ngi_config


@with_ngi_config
def update_charon_with_local_jobs_status(
        config=None, log=None, tracking_connector=None, charon_connector=None, **kwargs):
    log = log or minimal_logger(__name__, debug=True)

    tracking_connector = tracking_connector or TrackingConnector(config, log)
    charon_connector = charon_connector or CharonConnector(config, log)
    log.debug("updating Charon status for locally tracked jobs")
    for analysis in tracking_connector.tracked_analyses():
        log.debug("checking status for analysis of {}:{} with {}:{}, having {}".format(
            analysis.project_id,
            analysis.sample_id,
            analysis.engine,
            analysis.workflow,
            "pid {}".format(analysis.process_id) if analysis.process_id is not None else
            "sbatch job id {}".format(analysis.slurm_job_id)))
        process_status = get_analysis_status(analysis)
        log.debug(
            "{} with id {} has status {}".format(
                "process" if analysis.process_id is not None else "job",
                analysis.process_id or analysis.slurm_job_id,
                str(process_status)))
        charon_connector.set_sample_analysis_status(
            analysis.project_id,
            analysis.sample_id,
            charon_connector.analysis_status_from_process_status(process_status),
            recurse=True)
        log.debug(
            "{}removing from local tracking db".format(
                "not " if process_status is ProcessRunning else ""))
        remove_analysis(analysis, process_status, tracking_connector, force=False)


def get_analysis_status(analysis):
    status_type = ProcessStatus if analysis.process_id is not None else JobStatus
    processid_or_jobid = analysis.process_id or analysis.slurm_job_id
    exit_code_path = SarekAnalysis.get_analysis_type_for_workflow(analysis.workflow).sample_analysis_exit_code_path(
        analysis.project_base_path, analysis.project_id, analysis.sample_id)
    return status_type.get_type_from_processid_and_exit_code_path(processid_or_jobid, exit_code_path)


def remove_analysis(analysis, process_status, tracking_connector, force=False):
    # only remove the analysis if the process is not running or it should be forced
    if process_status == ProcessRunning and not force:
        return
    tracking_connector.remove_analysis(analysis)
