import collections
import datetime
import fnmatch
import glob
import os
import shutil
import subprocess
import yaml

from ngi_pipeline.conductor.classes import NGIProject
from ngi_pipeline.database.classes import CharonSession
from ngi_pipeline.log.loggers import log_process_non_blocking, minimal_logger
from ngi_pipeline.utils.filesystem import execute_command_line, rotate_file, safe_makedir

LOG = minimal_logger(__name__)


def launch_piper_job(command_line, project, log_file_path=None):
    """Launch the Piper command line.

    :param str command_line: The command line to execute
    :param Project project: The Project object (needed to set the CWD)

    :returns: The subprocess.Popen object for the process
    :rtype: subprocess.Popen
    """
    working_dir = os.path.join(project.base_path, "ANALYSIS", project.dirname)
    file_handle = None
    if log_file_path:
        try:
            file_handle = open(log_file_path, 'w')
        except Exception as e:
            LOG.error('Could not open log file "{}"; reverting to standard '
                      'logger (error: {})'.format(log_file_path, e))
            log_file_path = None
    popen_object = execute_command_line(command_line, cwd=working_dir, shell=True,
                                        stdout=(file_handle or subprocess.PIPE),
                                        stderr=(file_handle or subprocess.PIPE))
    if not log_file_path:
        log_process_non_blocking(popen_object.stdout, LOG.info)
        log_process_non_blocking(popen_object.stderr, LOG.warn)
    return popen_object


def _find_existing_analysis_results(project_obj, sample_obj=None, analysis_type=None):
    """
    Business logic for locating pre-existing analysis results on disk.

    :param NGIProject project_obj: The NGIProject object with relevant NGISamples
    :param NGISample sample_obj: The relevant NGISample object
    :param str analysis_type: Indicating the analysis type to search for. Possible values: all, genotype, variantcall
        default is "all"
    :return: A list of paths to located pre-existing analysis results
    :rtype list:
    :raises TypeError: if an invalid analysis_type is specified
    """
    analysis_type = analysis_type or "all"
    if analysis_type not in ["all", "genotype", "variantcall"]:
        raise TypeError("{} is not a valid analysis type to locate results for".format(analysis_type))

    project_dir_path = os.path.join(project_obj.base_path, "ANALYSIS",
                                    project_obj.project_id, "piper_ngi")
    project_dir_pattern = os.path.join(project_dir_path, "[0123456789][0123456789]_*")
    LOG.debug("Searching for previous {} analysis output files in "
              "{}".format(analysis_type, project_dir_path))

    analysis_result_files = []
    for sample in project_obj:
        # skip if a specific sample is requested
        if sample_obj and sample.name != sample_obj.name:
            continue
        analysis_result_files.extend(
            glob.glob(os.path.join(project_dir_pattern,
                                   "{}.*".format(sample.name))))
        # include the hidden files as well
        analysis_result_files.extend(
            glob.glob(os.path.join(project_dir_pattern,
                                   ".{}.*".format(sample.name))))

    def _piper_output_dir(output_file):
        return os.path.relpath(output_file, project_dir_path).split(os.sep)[0]

    # filter the analysis result files based on analysis_type and return the list
    return filter(lambda x: analysis_type == "all" or (
        _piper_output_dir(x).endswith("genotype_concordance") and analysis_type == "genotype") or (
        not _piper_output_dir(x).endswith("genotype_concordance") and analysis_type == "variantcall"),
                  analysis_result_files)


def _remove_existing_analysis_results(project_obj, sample_obj=None, analysis_type=None):
    analysis_result_files = _find_existing_analysis_results(
        project_obj,
        sample_obj=sample_obj,
        analysis_type=analysis_type)
    LOG.info("Deleting {} files for samples {}".format(
        analysis_type or "all",
        ", ".join(project_obj.samples)))
    for analysis_result_file in analysis_result_files:
        LOG.debug("Deleting file {}".format(analysis_result_file))
        try:
            if os.path.isdir(analysis_result_file):
                shutil.rmtree(analysis_result_file)
            else:
                os.remove(analysis_result_file)
        except OSError as e:
            LOG.warn("Error when removing {}: {}".format(analysis_result_file, e))
    if not analysis_result_files:
        LOG.debug("No {} analysis files found to delete for project {} / samples {}".format(
            project_obj,
            ", ".join(project_obj.samples)))


def find_previous_genotype_analyses(project_obj, sample_obj):
    analysis_result_files = _find_existing_analysis_results(
        project_obj,
        analysis_type="genotype",
        sample_obj=sample_obj)
    return any((
        os.path.exists(
            os.path.join(
                os.path.dirname(
                    ar_file),
                ".{}.done".format(
                    os.path.basename(
                        ar_file)))) for ar_file in analysis_result_files))


def find_previous_sample_analyses(project_obj, sample_obj=None, include_genotype_files=False):
    """Find analysis results for a sample, including .failed and .done files.
    Doesn't throw an error if it can't read a directory.

    :param NGIProject project_obj: The NGIProject object with relevant NGISamples
    :param bool include_genotype_files: Include genotyping files (default False)

    :returns: A list of files
    :rtype: list
    """
    return _find_existing_analysis_results(
        project_obj,
        analysis_type="all" if include_genotype_files else "variantcall",
        sample_obj=sample_obj)


def remove_previous_genotype_analyses(project_obj):
    """Remove genotype concordance analysis results for a sample, including
    .failed and .done files.
    Doesn't throw an error if it can't read a directory, but does if it can't
    delete a file it knows about.

    :param NGIProject project_obj: The NGIProject object with relevant NGISamples

    :returns: Nothing
    :rtype: None
    """
    _remove_existing_analysis_results(project_obj, analysis_type="genotype")


def remove_previous_sample_analyses(project_obj, sample_obj=None):
    """Remove analysis results for a sample, including .failed and .done files.
    Doesn't throw an error if it can't read a directory, but does if it can't
    delete a file it knows about.

    :param NGIProject project_obj: The NGIProject object with relevant NGISamples
    :param NGISample sample_obj: The relevant NGISample object 

    :returns: Nothing
    :rtype: None
    """
    _remove_existing_analysis_results(project_obj, sample_obj=sample_obj, analysis_type="variantcall")


def get_finished_seqruns_for_sample(project_id, sample_id,
                                    include_failed_libpreps=False):
    """Find all the finished seqruns for a particular sample.

    :param str project_id: The id of the project
    :param str sample_id: The id of the sample

    :returns: A dict of {libprep_01: [seqrun_01, ..., seqrun_nn], ...}
    :rtype: dict
    """
    charon_session = CharonSession()
    sample_libpreps = charon_session.sample_get_libpreps(projectid=project_id,
                                                         sampleid=sample_id)
    libpreps = collections.defaultdict(list)
    for libprep in sample_libpreps['libpreps']:
        if libprep.get('qc') != "FAILED" or include_failed_libpreps:
            libprep_id = libprep['libprepid']
            for seqrun in charon_session.libprep_get_seqruns(projectid=project_id,
                                                             sampleid=sample_id,
                                                             libprepid=libprep_id)['seqruns']:
                seqrun_id = seqrun['seqrunid']
                aln_status = charon_session.seqrun_get(projectid=project_id,
                                                       sampleid=sample_id,
                                                       libprepid=libprep_id,
                                                       seqrunid=seqrun_id).get('alignment_status')
                if aln_status == "DONE":
                    libpreps[libprep_id].append(seqrun_id)
                else:
                    LOG.debug('Skipping seqrun "{}" due to alignment_status '
                              '"{}"'.format(seqrun_id, aln_status))
        else:
            LOG.info('Skipping libprep "{}" due to qc status '
                     '"{}"'.format(libprep, libprep.get("qc")))
    return dict(libpreps)

def get_valid_seqruns_for_sample(project_id, sample_id,
                                 include_failed_libpreps=False,
                                 include_done_seqruns=False,
                                 status_field="alignment_status"):
    """Find all the valid seqruns for a particular sample.

    :param str project_id: The id of the project
    :param str sample_id: The id of the sample
    :param bool include_failed_libpreps: Include seqruns for libreps that have failed QC
    :param bool include_done_seqruns: Include seqruns that are already marked DONE

    :returns: A dict of {libprep_01: [seqrun_01, ..., seqrun_nn], ...}
    :rtype: dict

    :raises ValueError: If status_field is not a valid value
    """
    valid_status_values = ("alignment_status", "genotype_status",)
    if status_field not in valid_status_values:
        raise ValueError('"status_field" argument must be one of {} '
                         '(value passed was "{}")'.format(", ".join(valid_status_values),
                                                          status_field))
    charon_session = CharonSession()
    sample_libpreps = charon_session.sample_get_libpreps(projectid=project_id,
                                                         sampleid=sample_id)
    libpreps = collections.defaultdict(list)
    for libprep in sample_libpreps['libpreps']:
        if libprep.get('qc') != "FAILED" or include_failed_libpreps:
            libprep_id = libprep['libprepid']
            for seqrun in charon_session.libprep_get_seqruns(projectid=project_id,
                                                             sampleid=sample_id,
                                                             libprepid=libprep_id)['seqruns']:
                seqrun_id = seqrun['seqrunid']
                try:
                    aln_status = charon_session.seqrun_get(projectid=project_id,
                                                           sampleid=sample_id,
                                                           libprepid=libprep_id,
                                                           seqrunid=seqrun_id)[status_field]
                except KeyError:
                    LOG.error('Field "{}" not available for seqrun "{}" in Charon '
                              'for project "{}" / sample "{}". Including as '
                              'valid.'.format(status_field, seqrun_id,
                                              project_id, sample_id))
                    aln_status = None
                if aln_status != "DONE" or include_done_seqruns:
                    libpreps[libprep_id].append(seqrun_id)
                else:
                    LOG.info('Skipping seqrun "{}" due to {}'
                             '"{}"'.format(seqrun_id,status_field, aln_status))
        else:
            LOG.info('Skipping libprep "{}" due to qc status '
                     '"{}"'.format(libprep, libprep.get("qc")))
    return dict(libpreps)


def record_analysis_details(project, job_identifier):
    """Write a yaml file enumerating exactly which fastq files we've started
    analyzing.
    """
    output_file_path = os.path.join(project.base_path, "ANALYSIS",
                                    project.dirname, "piper_ngi","logs",
                                    "{}.files".format(job_identifier))
    analysis_dict = {}
    proj_dict = analysis_dict[project.dirname] = {}
    for sample in project:
        samp_dict = proj_dict[sample.name] = {}
        for libprep in sample:
            lib_dict = samp_dict[libprep.name] = {}
            for seqrun in libprep:
                lib_dict[seqrun.name] = seqrun.fastq_files
    rotate_file(output_file_path)
    safe_makedir(os.path.dirname(output_file_path))
    with open(output_file_path, 'w') as f:
        f.write(yaml.dump(analysis_dict))


def create_project_obj_from_analysis_log(project_name, project_id,
                                         project_base_path, sample_id, workflow):
    """Using the log of seqruns used for a sample analysis, recreate a project
    object with relevant sample, libprep, and seqrun objects.
    """
    analysis_log_filename = "{}-{}-{}.files".format(project_id, sample_id, workflow)
    analysis_log_path = os.path.join(project_base_path, "ANALYSIS",
                                     project_id, "piper_ngi", "logs", analysis_log_filename)
    with open(analysis_log_path, 'r') as f:
        analysis_dict = yaml.load(f)
    project_obj = NGIProject(name=project_name, dirname=project_id,
                             project_id=project_id, base_path=project_base_path)
    sample_obj = project_obj.add_sample(sample_id, sample_id)
    for libprep_name, seqrun_dict in analysis_dict[project_id][sample_id].items():
        libprep_obj = sample_obj.add_libprep(libprep_name, libprep_name)
        for seqrun_name in seqrun_dict.keys():
            libprep_obj.add_seqrun(seqrun_name, seqrun_name)
    return project_obj


def check_for_preexisting_sample_runs(project_obj, sample_obj,
                                      restart_running_jobs, restart_finished_jobs,
                                      status_field="alignment_status"):
    """If any analysis is undergoing or has completed for this sample's
    seqruns, raise a RuntimeError.

    :param NGIProject project_obj: The project object
    :param NGISample sample_obj: The sample object
    :param boolean restart_running_jobs: command line parameter
    :param boolean restart_finished_jobs: command line parameter
    :param str status_field: The field to check in Charon (seqrun level)

    :raise RuntimeError if the status is RUNNING or DONE and the flags do not allow to continue
    """
    project_id = project_obj.project_id
    sample_id = sample_obj.name
    charon_session = CharonSession()
    sample_libpreps = charon_session.sample_get_libpreps(projectid=project_id,
                                                         sampleid=sample_id)
    for libprep in sample_libpreps['libpreps']:
        libprep_id = libprep['libprepid']
        for seqrun in charon_session.libprep_get_seqruns(projectid=project_id,
                                                         sampleid=sample_id,
                                                         libprepid=libprep_id)['seqruns']:
            seqrun_id = seqrun['seqrunid']
            aln_status = charon_session.seqrun_get(projectid=project_id,
                                                   sampleid=sample_id,
                                                   libprepid=libprep_id,
                                                   seqrunid=seqrun_id).get(status_field)
            if (aln_status == "RUNNING" or aln_status == "UNDER_ANALYSIS" and \
                not restart_running_jobs) or \
                (aln_status == "DONE" and not restart_finished_jobs):
                raise RuntimeError('Project/Sample "{}/{}" has a preexisting '
                                   'seqrun "{}" with status "{}"'.format(project_obj,
                                                                         sample_obj,
                                                                         seqrun_id,
                                                                         aln_status))


SBATCH_HEADER = """#!/bin/bash -l

#SBATCH -A {slurm_project_id}
#SBATCH -p {slurm_queue}
#SBATCH -n {num_cores}
#SBATCH -t {slurm_time}
#SBATCH -J {job_name}
#SBATCH -o {slurm_out_log}
#SBATCH -e {slurm_err_log}
"""

def create_sbatch_header(slurm_project_id, slurm_queue, num_cores, slurm_time,
                         job_name, slurm_out_log, slurm_err_log):
    """
    :param str slurm_project_id: The project ID to use for accounting (e.g. "b2013064")
    :param str slurm_queue: "node" or "core"
    :param int num_cores: How many cores to use (max 16 at the moment)
    :param str slurm_time: The time for to schedule (e.g. "0-12:34:56")
    :param str job_name: The name to use for the job (e.g. "Your Mom")
    :param str slurm_out_log: The path to use for the slurm stdout log
    :param str slurm_err_log: The path to use for the slurm stderr log
    """
    ## TODO check how many cores are available for a given slurm queue
    if num_cores > 16: num_cores = 16
    return SBATCH_HEADER.format(slurm_project_id=slurm_project_id,
                                slurm_queue=slurm_queue,
                                num_cores=num_cores,
                                slurm_time=slurm_time,
                                job_name=job_name,
                                slurm_out_log=slurm_out_log,
                                slurm_err_log=slurm_err_log)


def add_exit_code_recording(cl, exit_code_path):
    """Takes a command line and returns it with increased pizzaz"""
    record_exit_code = "; echo $? > {}".format(exit_code_path)
    if type(cl) is list:
        # This should work, right? Right
        cl = " ".join(cl)
    return cl + record_exit_code


def create_log_file_path(workflow_subtask, project_base_path, project_name,
                         project_id=None,sample_id=None, libprep_id=None, seqrun_id=None):
    file_base_pathname = _create_generic_output_file_path(workflow_subtask,
                                                          project_base_path,
                                                          project_name,
                                                          project_id,
                                                          sample_id,
                                                          libprep_id,
                                                          seqrun_id)
    return file_base_pathname + ".log"


def create_exit_code_file_path(workflow_subtask, project_base_path, project_name, project_id,
                               sample_id=None, libprep_id=None, seqrun_id=None):
    file_base_pathname = _create_generic_output_file_path(workflow_subtask,
                                                          project_base_path,
                                                          project_name,
                                                          project_id,
                                                          sample_id,
                                                          libprep_id,
                                                          seqrun_id)
    return file_base_pathname + ".exit"


def _create_generic_output_file_path(workflow_subtask, project_base_path, project_name, project_id,
                                     sample_id=None, libprep_id=None, seqrun_id=None):
    base_path = os.path.join(project_base_path, "ANALYSIS", project_id, "piper_ngi","logs")
    file_name = project_id
    if sample_id:
        file_name += "-{}".format(sample_id)
        if libprep_id:
            file_name += "-{}".format(libprep_id)
            if seqrun_id:
                file_name += "-{}".format(seqrun_id)
    file_name += "-{}".format(workflow_subtask)
    return os.path.join(base_path, file_name)
