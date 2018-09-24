"""Keeps track of running workflow processes"""
import contextlib
import glob
import json
import os
import re
import shelve

from ngi_pipeline.database.classes import CharonSession, CharonError
from ngi_pipeline.database.communicate import get_project_id_from_name
from ngi_pipeline.log.loggers import minimal_logger
from ngi_pipeline.utils.classes import with_ngi_config
from ngi_pipeline.utils.parsers import STHLM_UUSNP_SEQRUN_RE, \
                                       STHLM_UUSNP_SAMPLE_RE

LOG = minimal_logger(__name__)

def get_all_tracked_processes(config=None):
    """Returns all the processes that are being tracked locally,
    which is to say all the processes that have a record in our local
    process_tracking database.

    :param dict config: The parsed configuration file (optional)

    :returns: The dict of the entire database
    :rtype: dict
    """
    with get_shelve_database(config) as db:
        db_dict = {job_name: job_obj for job_name, job_obj in db.iteritems()}
    return db_dict



# This can be run intermittently to track the status of jobs and update the database accordingly,
# as well as to remove entries from the local database if the job has completed (but ONLY ONLY ONLY
# once the status has been successfully written to Charon!!)
## FIXME
@with_ngi_config
def check_update_jobs_status(projects_to_check=None, config=None, config_file_path=None):
    """Check and update the status of jobs associated with workflows/projects;
    this goes through every record kept locally, and if the job has completed
    (either successfully or not) AND it is able to update Charon to reflect this
    status, it deletes the local record.

    :param list projects_to_check: A list of project names to check (exclusive, optional)
    :param dict config: The parsed NGI configuration file; optional.
    :param list config_file_path: The path to the NGI configuration file; optional.
    """
    db_dict = get_all_tracked_processes()
    for job_name, project_dict in db_dict.iteritems():
        LOG.info("Checking workflow {} for project {}...".format(project_dict["workflow"],
                                                                 job_name))
        return_code = project_dict["p_handle"].poll()
        if return_code is not None:
            # Job finished somehow or another; try to update database.
            LOG.info('Workflow "{}" for project "{}" completed '
                     'with return code "{}". Attempting to update '
                     'Charon database.'.format(project_dict['workflow'],
                                               job_name, return_code))
            # Only if we succesfully write to Charon will we remove the record
            # from the local db; otherwise, leave it and try again next cycle.
            try:
                project_id = project_dict['project_id']
                if projects_to_check and project_id not in projects_to_check:
                    continue
                ## MARIO FIXME
                ### THIS IS NOT REALLY CORRECT, HERE TRIGGER KNOWS DETAILS ABOUT ENGINE!!!!
                if project_dict["workflow"] == "dna_alignonly":
                    #in this case I need to update the run level infomration
                    #I know that:
                    # I am running Piper at flowcell level, I need to know the folder where results are stored!!!
                    write_to_charon_alignment_results(job_name, return_code, project_dict["run_dir"])
                elif project_dict["workflow"] == "NGI":
                    write_to_charon_NGI_results(job_name, return_code, project_dict["run_dir"])
                else:
                    write_status_to_charon(project_id, return_code)
                LOG.info("Successfully updated Charon database.")
                try:
                    # This only hits if we succesfully update Charon
                    remove_record_from_local_tracking(job_name)
                except RuntimeError:
                    LOG.error(e)
                    continue
            except RuntimeError as e:
                LOG.warning(e)
                continue
        else:
            ## FIXME
            #this code duplication can be avoided
            if project_dict["workflow"] == "dna_alignonly":
                #in this case I need to update the run level infomration
                write_to_charon_alignment_results(job_name, return_code, project_dict["run_dir"])
            elif project_dict["workflow"] == "NGI":
                    write_to_charon_NGI_results(job_name, return_code, project_dict["run_dir"])

            LOG.info('Workflow "{}" for project "{}" (pid {}) '
                     'still running.'.format(project_dict['workflow'],
                                             job_name,
                                             project_dict['p_handle'].pid))


def remove_record_from_local_tracking(project, config=None):
    """Remove a record from the local tracking database.

    :param NGIProject project: The NGIProject object
    :param dict config: The parsed configuration file (optional)

    :raises RuntimeError: If the record could not be deleted
    """
    LOG.info('Attempting to remove local process record for '
             'project "{}"...'.format(project))
    with get_shelve_database(config) as db:
        try:
            del db[project]
            LOG.info("...successfully removed.")
        except KeyError:
            error_msg = ('Project "{}" not found in local process '
                         'tracking database.'.format(project))
            LOG.error(error_msg)
            raise RuntimeError(error_msg)


def write_status_to_charon(project_id, return_code):
    """Update the status of a workflow for a project in the Charon database.

    :param NGIProject project_id: The name of the project
    :param int return_code: The return code of the workflow process

    :raises RuntimeError: If the Charon database could not be updated
    """
    ## Consider keeping on CharonSession open. What's the time savings?
    charon_session = CharonSession()
    ## Is "CLOSED" correct here?
    status = "CLOSED" if return_code is 0 else "FAILED"
    try:
        charon_session.project_update(project_id, status=status)
    except CharonError as e:
        error_msg = ('Failed to update project status to "{}" for "{}" '
                     'in Charon database: {}'.format(status, project_id, e))
        raise RuntimeError(error_msg)


## FIXME This needs to be moved to the engine_ngi or else some generic format needs to be created
##       and a dict passed back to this function. Something like that.
## BUG is this just for seqrun results? Doesn't work for sample??
def write_to_charon_NGI_results(job_id, return_code, run_dir):
    """Update the status of a sequencing run after alignment.

    :param NGIProject project_id: The name of the project, sample, lib prep, flowcell id
    :param int return_code: The return code of the workflow process
    :param string run_dir: the directory where results are stored (I know that I am running piper)

    :raises RuntimeError: If the Charon database could not be updated
    """
    charon_session = CharonSession()
    # Consider moving this mapping to the CharonSession object or something
    if return_code is None:
        status = "RUNNING"
    elif return_code == 0:
        status = "DONE"
    else:
        ## TODO we need to differentiate between COMPUTATION_FAILED and DATA_FAILED
        ##      also there is IGNORE?
        status = "COMPUTATION_FAILED"
    try:
        m_dict = STHLM_UUSNP_SAMPLE_RE.match(job_id).groupdict()
        #m_dict = re.match(r'?P<project_name>\w\.\w+_\d+_\d+|\w{2}-\d+)_(?P<sample_id>[\w-]+)_(?P<libprep_id>\w|\w{2}\d{3}_\2)_(?P<seqrun_id>\d{6}_\w+_\d{4}_.{10})', job_id).groupdict()
        project_id = get_project_id_from_name(m_dict['project_name'])
        sample_id = m_dict['sample_id']
    except (TypeError, AttributeError):
        error_msg = "Could not parse project/sample ids from job id \"{}\"; cannot update Charon with results!".format(job_id)
        raise RuntimeError(error_msg)
    try:
        charon_session.sample_update(project_id, sample_id, status=status)
    except CharonError as e:
        error_msg = ('Failed to update sample status to "{}" for sample "{}" '
                     'in Charon database: {}'.format(status, project_id, sample_id, e))
        raise RuntimeError(error_msg)


## moved to piper_ngi.local_process_tracking
def record_process_sample(p_handle, workflow, project, sample, analysis_module, analysis_dir, config=None):
    LOG.info('Recording process id "{}" for project "{}", sample "{}", '
             'workflow "{}"'.format(p_handle.pid, project, sample, workflow))
    project_dict = { "workflow": workflow,
                     "p_handle": p_handle,
                     "analysis_module": analysis_module.__name__,
                     "project_id": project.project_id,
                     "run_dir": analysis_dir
                   }
    with get_shelve_database(config) as db:
        db_key = "{}_{}".format(project, sample)
        if db_key in db:
            error_msg = ('Project "{}" / sample "{}" has an entry in the '
                         'local db. '.format(project, sample))
            raise RuntimeError(error_msg)
        else:
            db[db_key] = project_dict
            LOG.info('Successfully recorded process id "{}" for project "{}" / '
                     'sample "{}" / workflow "{}"'.format(p_handle.pid,
                                                          project,
                                                          sample,
                                                          workflow))


# Don't switch the order of these or you'll break everything
@contextlib.contextmanager
@with_ngi_config
def get_shelve_database(config=None, config_file_path=None):
    """Context manager for opening the local process tracking database.
    Closes the db automatically on exit.
    """
    try:
        database_path = config["database"]["record_tracking_db_path"]
    except KeyError as e:
        error_msg = ("Could not get path to process tracking database "
                     "from provided configuration: key missing: {}".format(e))
        raise KeyError(error_msg)
    db = shelve.open(database_path)
    try:
        yield db
    finally:
        db.close()
