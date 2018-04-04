import os

from ngi_pipeline.conductor.classes import NGIProject
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

        project_obj = project_from_analysis(analysis)
        sample_obj = list(filter(lambda x: x.name == analysis.sample_id, project_obj)).pop()
        restrict_to_seqruns = {}
        for libprep_obj in sample_obj:
            restrict_to_seqruns[libprep_obj.name] = list(map(lambda x: x.name, libprep_obj))
        restrict_to_libpreps = restrict_to_seqruns.keys()

        charon_connector.set_sample_analysis_status(
            analysis.project_id,
            analysis.sample_id,
            charon_connector.analysis_status_from_process_status(process_status),
            recurse=True,
            restrict_to_libpreps=restrict_to_libpreps,
            restrict_to_seqruns=restrict_to_seqruns)
        log.debug(
            "{}removing from local tracking db".format(
                "not " if process_status is ProcessRunning else ""))
        remove_analysis(analysis, process_status, tracking_connector, force=False)


def _project_from_fastq_file_paths(fastq_file_paths):
    """
    recreate the project object from fastq files being analysed
    :param fastq_file_paths:
    :return:
    """
    project_obj = None
    for fastq_file_path in fastq_file_paths:
        seqrun_path, fastq_file_name = os.path.split(fastq_file_path)
        libprep_path, seqrun_name = os.path.split(seqrun_path)
        sample_path, libprep_name = os.path.split(libprep_path)
        project_path, sample_name = os.path.split(sample_path)
        project_data_path, project_name = os.path.split(project_path)
        project_base_path = os.path.dirname(project_data_path)

        project_obj = project_obj or NGIProject(project_name, project_name, project_name, project_base_path)
        sample_obj = project_obj.add_sample(sample_name, sample_name)
        libprep_obj = sample_obj.add_libprep(libprep_name, libprep_name)
        seqrun_obj = libprep_obj.add_seqrun(seqrun_name, seqrun_name)
        seqrun_obj.add_fastq_files(fastq_file_name)
    return project_obj


def project_from_analysis(analysis):
    analysis_type = SarekAnalysis.get_analysis_type_for_workflow(analysis.workflow)
    tsv_file_path = analysis_type.sample_analysis_tsv_file(
        analysis.project_base_path, analysis.project_id, analysis.sample_id)
    fastq_file_paths = analysis_type.fastq_files_from_tsv_file(tsv_file_path)
    return _project_from_fastq_file_paths(fastq_file_paths)


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
