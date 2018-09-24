#!/bin/env python
"""
Build a bcbio-nextgen run config using info from the SampleSheet.csv and
project info from the StatusDB.
"""
from __future__ import print_function

import argparse
import collections
import couchdb
import datetime
import glob
import functools
import os
import re
import shlex
import shutil
import subprocess
import sys
import yaml

from bcbio.workflow import template
from scilifelab.utils.config import load_yaml_config_expand_vars
from textwrap import dedent

from scilifelab.log import minimal_logger
# Set up logging for this script
LOG = minimal_logger(__name__)

## TODO better error messages for log

## TODO add these to the cl
## Maybe this should just be assumed to have been done via e.g. .bashrc?
PRERUN_ENV_SETUP_CMDS = [ "source activate bcbio-nextgen",
                          "load_modules"]
LAUNCH_METHODS = (
                    'sbatch',
                    'local',
                 )

def launch_pipeline(run_config, launch_method):
    if launch_method.lower() not in LAUNCH_METHODS:
        LOG.warning("Launch method \"{}\" not available for run config {}. "\
                 "Valid choices are: {}".format(launch_method, run_config,
                                                ", ".join(LAUNCH_METHODS)))
        return None
    cmd = build_bcbio_cmd(run_config)
    if launch_method.lower() == "sbatch":
        create_sbatch_file()
        #create_sbatch_file(sbatch_output_dir=os.path.dirname(run_config))
        queue_sbatch(sbatch_file)
    elif launc_method.lower() == "local":
        execute_locally_nonblocking(cmd)



def execute_locally_nonblocking(command="", work_dir=None):
    """Executes the command line locally via subprocess.Popen to allow the code
    to continue processing multiple commands.

    :param str command: The command to be executed via the shell.
    :param str work_dir: The directory in which to execute the command.`

    :returns: The subprocess object.
    :rtype: subprocess.Popen
    """
    return subprocess.Popen(shlex.split(command),
                            bufsize=-1,
                            ## TODO These should redirect to the logger
                            ## TODO stdout=None,
                            ## TODO stderr=None,
                            # This is redundant but doesn't hurt I suppose
                            #env=dict(os.environ, my_additional_env="some_value"),
                            cwd=work_dir,)


def execute_remote_ssh_nonblocking(command, work_dir=None, user=None, host=None):
    """Exectutes the command line as user@host. Uses subprocess.Popen
    to allow the code to continue processing multiple commands.

    :param str command: The command to be executed remotely via the shell.
    :param str work_dir: The directory in which to execute the command.
    :param str user: The user on the remote machine.
    :param str host: The remote host (e.g. milou.uppmax.uu.se).

    :returns: The subprocess object.
    :rtype: subprocess.Popen
    """
    # Need to use bash -l so that we get the full login environment, including
    # all the variables set in .bashrc, .bash_profile, etc.
    remote_execution_cmd = "ssh {user}@{host} 'bash -l -c \"{}\"'".format(command)
    return execute_locally_nonblocking(remote_execution_cmd, work_dir)


def copy_sbatch_and_queue(sbatch_file_local, user, host, sbatch_file_remote_dir=None):
    """Copies the passed sbatch file to the host as the user and queues it.

    :param str sbatch_file_local: The path to the sbatch file to copy.
    :param str sbatch_file_remote_dir: The remote directory to copy to.
    :param str user: The user on the remote host.
    :param str host: The remote host.

    :returns: The slurm job id of the submitted batch job.
    :rtype: str

    :raises RuntimeError: If we could not connect to the remote host or
                          the file could not be queued.
    :raises subprocess.CalledProcessError: If return value is nonzero.
    :raises OSError: If sbatch or rsync is not available on the system.
    """
    if not sbatch_file_remote_dir:
        sbatch_file_remote_dir = os.path.dirname(sbatch_file_local)
    sbatch_file_remote_abspath = sbatch_file_remote_dir + \
                                 os.path.basename(sbatch_file_local)
    copy_cmd = shlex.split(
            "rsync -aPvz -e ssh {sbatch_file_local} {user}@{host}:{sbatch_file_remote_dir}/".format(
            sbatch_file_local=sbatch_file_local,
            user=user, host=host,
            sbatch_file_remote_dir=sbatch_file_remote_dir
            ))
    execution_cmd = shlex.split("ssh {user}@{host} \"sbatch {}\"".format(sbatch_file_remote_abspath))
    try:
        subprocess.check_call(copy_cmd)
    except subprocess.CalledProcessError:
        raise RuntimeError("Could not copy sbatch file {} to remote host {}".format(
            sbatch_file_local, host))
    try:
        s = subprocess.check_output(execution_cmd)
        return re.match('Submitted batch job (\d+)', s).groups()[0]
    except KeyError:
        raise RuntimeError("Could not determine slurm job ID; job may not have queued.")


def queue_sbatch(sbatch_file_path):
    """Submit the sbatch file to the slurm queue.
    
    :param sbatch_file_path: The path to the sbatch file.
    
    :return: Return code 0 if successful
    :rtype: int
 
    :raises subprocess.CalledProcessError: If return value is nonzero
    :raises OSError: If sbatch is not available on the system
    """
    return subprocess.check_return(shlex.split("sbatch {}".format(sbatch_file_path)))


## TODO hmm
def create_sbatch_file(sbatch_dir, sample_name, template_name, command_line, 
                       uppmax_project_group="a2010002", timelimit="7-00", numcores=16):
    """Creates an sbatch file that will start the bcbio run.

    :param str sbatch_dir: The directory that will store the sbatch launcher file.
    :param str sample_name: The sample name from the input file (stripped of read
                            number and extensions)
    :param str template_name: The name of the template/pipeline (e.g. gatk-variant)
    :param str command_line: The command that will be executed via slurm.
    :param str uppmax_project_group: The UPPMAX project group name (e.g. a2010002).
                                     Default a2010002.
    :param str timelimit: The timelimit to use when scheduling the job.
                          Default 7 days.
    :param str numcores: The number of cores to use when processing.
                         Default 16.

    :returns: The path to the sbatch file
    :rtype: str
    """
    sbatch_text = dedent("""    #!/bin/sh
    #SBATCH -p core
    #SBATCH -J bcbio-launcher
    #SBATCH -o bcbio-launcher.out.%j
    #SBATCH -e bcbio-launcher.err.%j
    #SBATCH -A {project}
    #SBATCH -t {timelimit}
    #SBATCH -n {numcores}

    {cl}
    """).format(project=uppmax_project_group, timelimit=timelimit,
                cl=command_line, numcores=numcores)
    ## TODO need function to make sure name doesn't contain weird chars (*/\.?)
    template_name = template_name.replace(".yaml", "")
    output_sbatch_file = os.path.join(sbatch_dir,
        "{sample_name}{template_name}.sbatch".format(sample_name=sample_name,
                                                     template_name=template_name))
    if os.path.exists(output_sbatch_file):
        shutil.move(output_sbatch_file,
                    output_sbatch_file + ".bak{}".format(
                        datetime.datetime.now().strftime("%Y-%m-%d-%H-%M-%S")))
    with open(output_sbatch_file, 'w') as f:
        f.write(sbatch_text)
    return output_sbatch_file


def build_run_configs(samples_dir=None, config_path=None, output_dir=None, upload_dir=None):
    """With a finished directory, generate a configuration file for each
    (pair of) sample file(s).

    :param str samples_dir: The path to the directory containing the fastq files.
    :param str config_path: That path to the config file with statusdb auth info,
                            sbatch parameters, and template file locations
    :param str output_dir: Where to store the config & work dirs. Default is a
                           directory named after the sample (or set of samples).
    :param str upload_dir: Where to store the processing output, if different
                           than output_dir.

    :returns: A list of dicts representing runs to execute

    :raises RuntimeError: If unable to connect to StatusDB
    """
    # Load info from config file
    try:
        config_yaml = load_yaml_config_expand_vars(config_path)
    except (IOError, TypeError) as e:
        LOG.error("Could not open config file for reading; cannot proceed ({})".format(e))
        return []
    try:
        pipeline_mappings = config_yaml["method_to_pipeline_mappings"]
    except KeyError as e:
        LOG.error("Config file must provide library construction method " \
                  "-> pipelines mapping. Cannot proceed.")
        return []
    files = glob.glob("{}/*.fastq*".format(samples_dir))
    supported_genomes = set(config_yaml.get("supported_genomes", []))

    # Connect to database
    status_db_config = config_yaml.get("statusdb")
    try:
        templates = config_yaml["templates"]
    except KeyError as e:
        ## TODO Logging
        LOG.error("Config file must provide template file locations. Cannot proceed.")
        return []
    import ipdb; ipdb.set_trace()
    LOG.info("Trying to connect to StatusDB... ")
    ## TODO do we need the user /password info?
    couch = couchdb.Server("http://{user}:{password}@{url}:{port}".format(
                    user=status_db_config.get("username"),
                    password=status_db_config.get("password"),
                    url=status_db_config.get("url"),
                    port=status_db_config.get("port")))
    #if not couch:
    #    raise RuntimeError("Couldn't connect to StatusDB or "\
    #                       "config file lacked authentication information.")

    ## TODO This is where the connection actually can time out
    proj_db = couch['projects']
    supported_genomes = set(supported_genomes)
    file_pairs = find_fastq_read_pairs(input_files)

    samples_to_process = []
    for sample_basename, sample_files in file_pairs.items():
        try:
            project_id = get_project_id_from_filename(sample_basename)
        except ValueError as e:
            # Could not determine project id
            LOG.warning(e)
            continue
        if not output_dir:
            output_dir = os.path.join(os.path.dirname(os.path.abspath(sample_files[0])), "project_{}".format(sample_basename))
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
        if not upload_dir:
            upload_dir = os.path.join(output_dir, "final")
        project_data_json = get_project_data_for_id(project_id, proj_db)
        adapter_seqs      = project_data_json.get("adapter_sequences")
        # Determine the library construction method from the Project
        lib_const_method  = project_data_json.get("details", {}).get("library_construction_method")
        # If we have no matching template, just put it through qc
        pipeline_list     = pipeline_mappings.get(lib_const_method) or ("qc_pipeline_template")
        reference_genome  = project_data_json.get("reference_genome")
        if reference_genome.lower() not in supported_genomes:
            # Unsupported organisms cannot be aligned
            reference_genome = None
            pipeline_list = ("qc_pipeline_template")
        ## TODO At the moment we just have multiple run config files -- each pipeline runs separately.
        ##      I need to think about when this might not work, conflicts and so on.
        ##      May eventually need some way to merge the two template files
        ##      or have template.py modify both pipelines when it adds attributes.
        # This could instead construct and then yield the list of pipelines,
        # if in the future we have pipelines that must run serially
        for template_name in pipeline_list:
            # Get the path to the template file; default is always just qc pipeline
            template_path = templates.get(template_name) or templates.get("qc_pipeline_template")
            namespace_args = argparse.Namespace(template=template_path,
                                                input_files=sample_files,
                                                out_dir=output_dir,
                                                upload_dir=upload_dir)
            if adapter_seqs or reference_genome:
                # Create the csv file that will be passed to the template creator
                project_csv = create_project_csv_from_dbinfo(sample_basename, output_dir,
                                                             adapter_seqs, reference_genome)
                namespace_args.__dict__["metadata"] = project_csv
            ## TODO bcbio.workflow.template doesn't expand $ENV_VARS -- fix
            config_file_path, work_dir = template.setup(namespace_args)
            samples_to_process.append({ 'sample_basename': sample_basename,
                                        'template_name': os.path.basename(template_name),
                                        'run_config': config_file_path,
                                        'work_dir': work_dir})
    return samples_to_process



## TODO load run params, dist params from configuration file
#def build_bcbio_cmd(bcbio_run_params=None, bcbio_dist_params=None):
def build_bcbio_cmd(run_config_path, sys_config_path=None, work_dir_path=None):
    """Builds a bcbio_nextgen.py command line to launch processing.

    :param str run_config_path: The path to the run configuration.
    :param str sys_config_path: The path to the config file containing
                                parallelization parameters (optional).
    :param str work_dir_path: The path to the work directory (optional).

    :returns: A str that can be executed via the shell.
    :rtype: str

    :raises KeyError: if a required bcbio_nextgen parameter is missing
                      (system_config, run_config)
    """
    if not sys_config:
        bcbio_run_params = {
               "system_config_path": "$SYSCONFIG",
               "run_config_path": run_config,
               "numcores": 16,}
        # See bcbio_nextgen.py --help for available options
        bcbio_dist_params = {"resources": ["account:b2013064","timelimit=4-00:00:00"]}
        if work_dir_path:
               bcbio_run_params.update({"work_dir": work_dir_path})

    # Note that this assumes that the environment executing this commmand will have
    # all the requisite modules loaded and the environment activated
    cl = ["bcbio_nextgen.py"]
    cl.append("-n {} ".format(bcbio_run_params.get("numcores") or "1"))
    if bcbio_run_params.get("work_dir"):
        cl.append("--workdir {} ".format(bcbio_run_params["work_dir"]))
    if bcbio_dist_params:
        cl.append("--scheduler slurm")
        cl.append("--paralleltype ipython")
        cl.append("--queue {}".format(bcbio_dist_params.get("queue") or "core"))
        cl.append("--timeout {}".format(bcbio_dist_params.get("timeout") or "600"))
        cl.append("--retries {}".format(bcbio_dist_params.get("retries") or "1"))
        for res_param in bcbio_dist_params.get("resources") or []:
            cl.append("--resources {}".format(res_param))
    cl.append(bcbio_run_params["system_config_path"])
    cl.append(bcbio_run_params["run_config_path"])
    cl_text = " ".join(cl)
    return cl_text
