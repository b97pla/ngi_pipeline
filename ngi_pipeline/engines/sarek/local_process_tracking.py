import os

from ngi_pipeline.conductor.classes import NGIProject
from ngi_pipeline.engines.sarek.database import CharonConnector, TrackingConnector
from ngi_pipeline.engines.sarek.models.sarek import SarekAnalysis
from ngi_pipeline.engines.sarek.models.sample import SarekAnalysisSample
from ngi_pipeline.engines.sarek.process import JobStatus, ProcessStatus, ProcessRunning, ProcessExitStatusSuccessful
from ngi_pipeline.log.loggers import minimal_logger
from ngi_pipeline.utils.classes import with_ngi_config


@with_ngi_config
def update_charon_with_local_jobs_status(
        config=None, log=None, tracking_connector=None, charon_connector=None, **kwargs):
    """
    Update Charon with the local changes in the SQLite tracking database.

    :param config: optional dict with configuration options. If not specified, the global configuration will be used
    instead
    :param log: optional log instance. If not specified, a new log instance will be created
    :param tracking_connector: optional connector to the tracking database. If not specified,
    a new connector will be created
    :param charon_connector: optional connector to the charon database. If not specified, a new connector will be
    created
    :param kwargs: placeholder for additional, unused options
    :return: None
    """
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

        # get an analysis instance corresponding to the Sarek workflow
        analysis_instance = SarekAnalysis.get_analysis_instance_for_workflow(
            analysis.workflow,
            config,
            log,
            charon_connector=charon_connector,
            tracking_connector=tracking_connector)

        # update the sample status
        update_charon_with_sample_local_job_status(analysis, analysis_instance, log)


def update_charon_with_sample_local_job_status(analysis, analysis_instance, log):
    # recreate a NGIProject object from the analysis
    project_obj = project_from_analysis(analysis, analysis_instance)
    # extract the sample object corresponding to the analysis entry
    sample_obj = list(filter(lambda x: x.name == analysis.sample_id, project_obj)).pop()
    analysis_sample = SarekAnalysisSample(project_obj, sample_obj, analysis_instance)

    # get the analysis status
    process_status = get_analysis_status(analysis, analysis_instance)
    log.debug(
        "{} with id {} has status {}".format(
            "process" if analysis.process_id is not None else "job",
            analysis.process_id or analysis.slurm_job_id,
            str(process_status)))

    # set the analysis status of the sample in charon, recursing to libpreps and seqruns
    report_analysis_status(analysis_sample, process_status, log)
    # set the analysis results of the sample in charon if the process has finished successfully
    report_analysis_results(analysis_sample, process_status, log)

    # remove the entry from the tracking database, but only if the process is not running
    log.debug(
        "{}removing from local tracking db".format(
            "not " if process_status is ProcessRunning else ""))
    remove_analysis(analysis, process_status, analysis_instance.tracking_connector, force=False)


def _get_libpreps_and_seqruns(analysis_sample):
    # extract the libpreps and seqruns from the sample object
    restrict_to_seqruns = {
        libprepid: analysis_sample.libprep_seqrun_ids(libprepid) for libprepid in analysis_sample.sample_libprep_ids()}
    return restrict_to_seqruns.keys(), restrict_to_seqruns


def report_analysis_status(analysis_sample, process_status, log):
    # set the status in charon, recursing into libpreps and seqruns
    restrict_to_libpreps, restrict_to_seqruns = _get_libpreps_and_seqruns(analysis_sample)
    log.debug("setting analysis status in charon to match '{}' for '{}' - '{}' - '{}' - '{}'".format(
        process_status,
        analysis_sample.projectid,
        analysis_sample.sampleid,
        ",".join(restrict_to_libpreps),
        ",".join(set([seqrun for seqruns in restrict_to_seqruns.values() for seqrun in seqruns]))
    ))
    charon_connector = analysis_sample.analysis_object.charon_connector
    charon_connector.set_sample_analysis_status(
        charon_connector.analysis_status_from_process_status(process_status),
        analysis_sample.projectid,
        analysis_sample.sampleid,
        recurse=True,
        restrict_to_libpreps=restrict_to_libpreps,
        restrict_to_seqruns=restrict_to_seqruns)


def report_analysis_results(analysis_sample, process_status, log):
    # only report the analysis results if the process exited successfully
    if process_status != ProcessExitStatusSuccessful:
        return

    analysis_metrics = analysis_sample.analysis_object.collect_analysis_metrics(analysis_sample)
    log.debug("setting analysis metrics {} in charon for sample '{}' in project '{}'".format(
        analysis_metrics, analysis_sample.sampleid, analysis_sample.projectid))
    # for convenience, map the metrics to the corresponding setter in the charon_connector
    for metric, setter in {
        "total_reads": "set_sample_total_reads",
        "autosomal_coverage": "set_sample_autosomal_coverage",
        "percent_duplication": "set_sample_duplication"
    }.items():
        # set the status in charon, skip recursing into libpreps and seqruns
        log.debug("setting {} to {} with {}".format(metric, analysis_metrics[metric], setter))
        getattr(
            analysis_sample.analysis_object.charon_connector,
            setter)(
            analysis_metrics[metric],
            analysis_sample.projectid,
            analysis_sample.sampleid)


def _project_from_fastq_file_paths(fastq_file_paths):
    """
    recreate the project object from a list of fastq file paths
    :param fastq_file_paths: list of fastq file paths, expected to be arranged in subfolders according to
    [/]path/to/project name/sample name/libprep name/seqrun name/fastq_file_name.fastq.gz

    :return: a NGIProject object recreated from the directory tree and fastq files
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


def project_from_analysis(analysis, analysis_instance):
    """
    Recreate a NGIProject object based on an analysis object fetched from the local tracking database. This will
    use the Sarek TSV file and the fastq files and directory structure.

    :param analysis: an analysis object fetched from a tracking connector
    :return: a NGIProject object recreated from the information in the analysis object
    """
    tsv_file_path = analysis_instance.sample_analysis_tsv_file(
        analysis.project_base_path, analysis.project_id, analysis.sample_id)
    fastq_file_paths = analysis_instance.fastq_files_from_tsv_file(tsv_file_path)
    return _project_from_fastq_file_paths(fastq_file_paths)


def get_analysis_status(analysis, analysis_instance):
    """
    Figure out the analysis status of an analysis object from the local tracking database. Will check the exit code
    written to an exit code file if the process is not running.

    :param analysis: an analysis object fetched from a tracking connector
    :return: a subclass of ProcessStatus representing the status of the process
    """
    status_type = ProcessStatus if analysis.process_id is not None else JobStatus
    processid_or_jobid = analysis.process_id or analysis.slurm_job_id
    exit_code_path = analysis_instance.sample_analysis_exit_code_path(
        analysis.project_base_path, analysis.project_id, analysis.sample_id)
    return status_type.get_type_from_processid_and_exit_code_path(processid_or_jobid, exit_code_path)


def remove_analysis(analysis, process_status, tracking_connector, force=False):
    """
    Remove the analysis object from the tracking database iff the process status is not ProcessRunning or force is True

    :param analysis: the analysis object to remove
    :param process_status: the process status to compare
    :param tracking_connector: the tracking connector to the tracking database
    :param force: if True, remove the entry even if the process is still running (default is False)
    :return: None
    """
    # only remove the analysis if the process is not running or it should be forced
    if process_status == ProcessRunning and not force:
        return
    tracking_connector.remove_analysis(analysis)


def parse_analysis_reports_and_update_database(projectid, config, log, charon_connector, tracking_connector):
    analysis_instance = SarekAnalysis.get_analysis_instance_for_project(
        projectid,
        config,
        log,
        charon_connector=charon_connector,
        tracking_connector=tracking_connector)
