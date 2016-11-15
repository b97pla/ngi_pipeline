"""The Piper automated launcher script."""
from __future__ import print_function

import os
import re
import subprocess
import sys
import time

from ngi_pipeline.conductor.classes import NGIProject, NGIChipGenotypes
from ngi_pipeline.engines.piper_ngi import command_creation_config, local_process_tracking, workflows
from ngi_pipeline.engines.piper_ngi import utils as sample_runs_handler
from ngi_pipeline.log.loggers import minimal_logger
from ngi_pipeline.utils import classes, filesystem, slurm

LOG = minimal_logger(__name__)


class PiperLauncher(object):

    def __init__(self,
                 project_obj,
                 sample_obj,
                 restart_finished_jobs=False,
                 restart_running_jobs=False,
                 keep_existing_data=False,
                 generate_bqsr_bam=False,
                 level="sample",
                 config=None,
                 **kwargs):

        # set fields passed as arguments
        self.project_obj = project_obj
        self.sample_obj = sample_obj
        self.restart_finished_jobs = restart_finished_jobs
        self.restart_running_jobs = restart_running_jobs
        self.keep_existing_data = keep_existing_data
        self.generate_bqsr_bam = generate_bqsr_bam
        self.level = level
        self.config = config

        # set default fields if not passed
        self.exec_mode = kwargs.get("exec_mode", "sbatch")
        # dependency injection stuff
        self.workflow_handler = kwargs.get("workflow_handler", workflows)
        self.local_process_tracking_handler = kwargs.get("local_process_tracking_handler", local_process_tracking)
        self.preexisting_sample_runs_handler = kwargs.get("preexisting_sample_runs_handler", sample_runs_handler)
        self.sample_runs_handler = kwargs.get("sample_runs_handler", sample_runs_handler)
        self.filesystem_handler = kwargs.get("filesystem_handler", filesystem)
        self.command_handler = kwargs.get("command_handler", command_creation_config)
        self.launch_handler = kwargs.get("launch_handler", sys.modules[__name__])
        self.slurm_handler = kwargs.get("slurm_handler", slurm)

        # set derived fields
        self.status_field = "alignment_status" if self.level == "sample" else "genotype_status"
        self.project_data_directory = os.path.join(self.project_obj.base_path, "DATA", self.project_obj.dirname)
        self.sample_data_directory = os.path.join(self.project_data_directory, self.sample_obj.dirname)

    def analyze(self):
        # get the list of workflow subtasks to perform
        workflow_subtasks = self.determine_workflow_subtasks()

        # assert that we support the exec_mode (e.g. sbatch or local)
        try:
            self.assert_exec_mode()
        except NotImplementedError as ne:
            LOG.error('Processing project "{}" / sample "{}" / workflow "{}" failed: {}'.format(
                self.project_obj,
                self.sample_obj,
                ",".join(workflow_subtasks),
                ne))
            return False

        # assert that the sample is ok to be started (this will raise an exception if that's not the case)
        try:
            self.assert_sample_should_be_started()
        except RuntimeError:
            raise

        # iterate over the workflow subtasks and start them, cleaning up as we go
        for workflow_subtask in workflow_subtasks:
            # assert that a workflow subtask is not running on the sample at this point,
            # and if it is, that it is killed if that's desired. otherwise skip this sample
            if not self.assert_sample_is_not_running(workflow_subtask):
                continue

            LOG.info('Launching "{}" analysis for sample "{}" in project "{}"'.format(
                workflow_subtask, self.sample_obj, self.project_obj))

            # clean up previous analysis results for this sample
            self.clean_previous_analysis_results()

            # determine the libpreps and seqruns that should be included in the analysis
            try:
                libpreps_seqruns = self.determine_seqruns_to_be_analyzed()
            except RuntimeError as err:
                LOG.error('Processing project "{}" / sample "{}" / workflow "{}" failed: {}'.format(
                    self.project_obj, self.sample_obj, workflow_subtask, err))
                continue

            # locate the files needed for the analysis
            try:
                chip_genotypes_for_analysis = self.locate_chip_genotype_files_for_project()
                fastq_files_for_analysis = self.locate_fastq_files_for_project_sample(valid_seqruns=libpreps_seqruns)
                old_files_for_analysis = self.locate_preexisting_data_to_include_in_analysis()
            except ValueError as ve:
                LOG.error('Processing project "{}" / sample "{}" / workflow "{}" failed: {}'.format(
                    self.project_obj, self.sample_obj, workflow_subtask, ve))
                continue

            # create the file paths for log and exit code
            exit_code_path = self.create_exit_code_path(workflow_subtask)
            log_file_path = self.create_log_file_path(workflow_subtask)

            # rotate the log file
            self.rotate_log_file(log_file_path)

            # set up a local project object based on the located files to use for command creation
            local_project_obj = self.create_local_copy_of_project_obj(
                fastq_files_for_analysis, chip_genotypes_for_analysis)

            # create the setup command and get the setup xml output path
            setup_command, setup_xml_path = self.create_setup_command(local_project_obj, workflow_subtask)

            # create the piper command
            piper_command = self.create_piper_command(
                local_project_obj, workflow_subtask, setup_xml_path, exit_code_path)

            # submit the commands for execution
            submitted_job_id = self.submit_commands(
                [setup_command, piper_command], workflow_subtask, old_files_for_analysis)

            # verify the job status of the submitted commands
            try:
                self.submitted_job_status(submitted_job_id)
            except RuntimeError as err:
                LOG.error('command file for sample {}/{} did not submit properly: {}'.format(
                    self.project_obj, self.sample_obj, err))

            # record the job details in the local tracking database
            try:
                self.record_job_status_to_databases(workflow_subtask, submitted_job_id)
            except (RuntimeError, ValueError) as err:
                LOG.error("error recording job status to database: {}".format(err))

    def assert_exec_mode(self):
        implemented_modes = ["sbatch"]
        if self.exec_mode not in implemented_modes:
            raise NotImplementedError("{} execution mode not implemented. Only {} available.".format(
                self.exec_mode, ",".join(implemented_modes)
            ))

    def assert_sample_should_be_started(self):
        try:
            self.preexisting_sample_runs_handler.check_for_preexisting_sample_runs(
                self.project_obj,
                self.sample_obj,
                self.restart_running_jobs,
                self.restart_finished_jobs,
                self.status_field)
        except RuntimeError:
            raise

    def assert_sample_is_not_running(self, workflow_subtask):
        if self.restart_running_jobs:
            self.local_process_tracking_handler.kill_running_sample_analysis(
                workflow_subtask, self.project_obj.project_id, self.sample_obj.name)
        return not self.local_process_tracking_handler.is_sample_analysis_running_local(
            workflow_subtask, self.project_obj.project_id, self.sample_obj.name)

    def clean_previous_analysis_results(self):
        if not self.keep_existing_data:
            if self.level == "sample":
                self.preexisting_sample_runs_handler.remove_previous_sample_analyses(self.project_obj, self.sample_obj)
            else:
                self.preexisting_sample_runs_handler.remove_previous_genotype_analyses(self.project_obj)

    def create_exit_code_path(self, workflow_subtask):
        return self.sample_runs_handler.create_exit_code_file_path(
            workflow_subtask,
            self.project_obj.base_path,
            self.project_obj.dirname,
            self.project_obj.project_id,
            sample_id=self.sample_obj.name)

    def create_local_copy_of_project_obj(self, fastq_files, chip_genotype_files):
        # initialize the local NGIProject
        local_project_obj = NGIProject(
            self.project_obj.name,
            self.project_obj.dirname,
            self.project_obj.project_id,
            self.project_obj.base_path,
            chip_genotypes=map(
                lambda f: NGIChipGenotypes(name=os.path.basename(f)),
                chip_genotype_files) if chip_genotype_files else None)

        # attach a sample object
        local_sample_obj = local_project_obj.add_sample(
            self.sample_obj.name,
            self.sample_obj.dirname)

        # create libpreps and seqruns based on supplied fastq files
        for fastq_file in fastq_files:
            # split the path backwards to get the libprep, seqrun and fastq names
            # this reverses the full path, splits on '/', takes the first three elements and reverses them back
            (fastq_name, seqrun_name, libprep_name) = map(
                lambda x: x[::-1],
                fastq_file[::-1].split(os.path.sep)[0:3])
            local_libprep_obj = local_sample_obj.add_libprep(name=libprep_name, dirname=libprep_name)
            local_seqrun_obj = local_libprep_obj.add_seqrun(name=seqrun_name, dirname=seqrun_name)
            local_seqrun_obj.add_fastq_files(fastq_name)
        return local_project_obj

    def create_log_file_path(self, workflow_subtask):
        return self.sample_runs_handler.create_log_file_path(
            workflow_subtask,
            self.project_obj.base_path,
            self.project_obj.dirname,
            self.project_obj.project_id,
            sample_id=self.sample_obj.name)

    def create_piper_command(self, local_project_obj, workflow_subtask, setup_xml_path, exit_code_path):
        return self.command_handler.build_piper_cl(
            project=local_project_obj,
            workflow_name=workflow_subtask,
            setup_xml_path=setup_xml_path,
            exit_code_path=exit_code_path,
            config=self.config,
            exec_mode=self.exec_mode,
            generate_bqsr_bam=self.generate_bqsr_bam)

    def create_setup_command(self, local_project_obj, workflow_subtask):
        return self.command_handler.build_setup_xml(
            local_project_obj,
            next(local_project_obj),
            workflow_subtask,
            self.exec_mode == "sbatch",
            self.config
        )

    def determine_seqruns_to_be_analyzed(self, include_failed_libpreps=False):
        valid_seqruns = self.preexisting_sample_runs_handler.get_valid_seqruns_for_sample(
            project_id=self.project_obj.project_id,
            sample_id=self.sample_obj.name,
            include_failed_libpreps=include_failed_libpreps,
            include_done_seqruns=self.restart_finished_jobs,
            status_field=self.status_field)
        if not valid_seqruns:
            raise ValueError('No valid libpreps/seqruns found for project/sample "{}/{}"'.format(
                self.project_obj, self.sample_obj))
        return valid_seqruns

    def determine_workflow_subtasks(self):
        return self.workflow_handler.get_subtasks_for_level(level=self.level)

    def locate_chip_genotype_files_for_project(self):
        chip_genotype_files = self.filesystem_handler.locate_chip_genotypes_in_dir(self.project_data_directory)
        return chip_genotype_files if chip_genotype_files else None

    def locate_fastq_files_for_project_sample(self, realpath=False, valid_seqruns=None):
        fastq_files = self.filesystem_handler.fastq_files_under_dir(
            self.sample_data_directory, realpath=realpath)
        if not fastq_files:
            raise ValueError('No valid fastq files found for project/sample {}/{}'.format(
                self.project_obj, self.sample_obj))
        # flatten the valid_seqruns dict to a list of key1/val1, key1/val2 combinations
        valid_seqrun_paths = [
                os.path.join(libprep, seqrun)
                for libprep, seqruns in valid_seqruns.items()
                for seqrun in seqruns] if valid_seqruns else []
        # filter fastq files against valid seqruns if specified,
        # this assumes a hierarchy like .../libprep/seqrun/file.fastq.gz
        return filter(
            lambda fq: any(
                map(
                    lambda pth: os.path.split(fq)[0].endswith(pth),
                    valid_seqrun_paths)),
            fastq_files)

    def locate_preexisting_data_to_include_in_analysis(self, include_genotype_files=True):
        return self.preexisting_sample_runs_handler.find_previous_sample_analyses(
            self.project_obj, self.sample_obj,  include_genotype_files=include_genotype_files)

    def record_job_status_to_databases(self, workflow_subtask, submitted_job_id):
        try:
            self.local_process_tracking_handler.record_process_sample(
                project=self.project_obj,
                sample=self.sample_obj,
                analysis_module_name="piper_ngi",
                workflow_subtask=workflow_subtask,
                config=self.config,
                **submitted_job_id)
        except (RuntimeError, ValueError):
            raise

    def rotate_log_file(self, log_file_path):
        try:
            self.filesystem_handler.rotate_file(log_file_path)
        except OSError:
            raise

    def submit_commands(self, *args):
        if self.exec_mode == "sbatch":
            return self.submit_sbatch_commands(*args)

    def submit_sbatch_commands(self, command_list, workflow_subtask, files_needed_for_analysis):
        slurm_job_id = self.launch_handler.sbatch_piper_sample(
            command_list,
            workflow_subtask,
            self.project_obj,
            self.sample_obj,
            restart_finished_jobs=self.restart_finished_jobs,
            files_to_copy=files_needed_for_analysis,
            config=self.config)
        return {"slurm_job_id": slurm_job_id}

    def submitted_job_status(self, *args):
        if self.exec_mode == "sbatch":
            return self.submitted_sbatch_job_status(*args)

    def submitted_sbatch_job_status(self, job_id):
        slurm_job_id = job_id.get("slurm_job_id")
        for i in xrange(10):
            try:
                job_status = self.slurm_handler.get_slurm_job_status(slurm_job_id)
                return job_status
            except ValueError:
                time.sleep(2)
        raise RuntimeError("slurm job id {} cannot be found".format(slurm_job_id))


@classes.with_ngi_config
def analyze(project, sample,
            exec_mode="sbatch",
            restart_finished_jobs=False,
            restart_running_jobs=False,
            keep_existing_data=False,
            level="sample",
            genotype_file=None,
            config=None, config_file_path=None,
            generate_bqsr_bam=False):
    """
    This is the "old" analyze method with preserved signature. Now it wraps around
    the PiperLauncher object and calls its analyze method
    """
    launcher = PiperLauncher(
        project,
        sample,
        restart_finished_jobs=restart_finished_jobs,
        restart_running_jobs=restart_running_jobs,
        keep_existing_data=keep_existing_data,
        generate_bqsr_bam=generate_bqsr_bam,
        level=level,
        config=config,
        exec_mode=exec_mode
    )
    try:
        launcher.analyze()
    except RuntimeError:
        raise
    except Exception:
        raise

@classes.with_ngi_config
def sbatch_piper_sample(command_line_list, workflow_name, project, sample,
                        libprep=None, restart_finished_jobs=False, files_to_copy=None,
                        config=None, config_file_path=None):
    """sbatch a piper sample-level workflow.

    :param list command_line_list: The list of command lines to execute (in order)
    :param str workflow_name: The name of the workflow to execute
    :param NGIProject project: The NGIProject
    :param NGISample sample: The NGISample
    :param dict config: The parsed configuration file (optional)
    :param str config_file_path: The path to the configuration file (optional)
    """
    job_identifier = "{}-{}-{}".format(project.project_id, sample, workflow_name)
    # Paths to the various data directories
    project_dirname = project.dirname
    perm_analysis_dir = os.path.join(project.base_path, "ANALYSIS", project_dirname, "piper_ngi", "")
    scratch_analysis_dir = os.path.join("$SNIC_TMP/ANALYSIS/", project_dirname, "piper_ngi", "")
    #ensure that the analysis dir exists
    filesystem.safe_makedir(perm_analysis_dir)
    try:
        slurm_project_id = config["environment"]["project_id"]
    except KeyError:
        raise RuntimeError('No SLURM project id specified in configuration file '
                           'for job "{}"'.format(job_identifier))
    slurm_queue = config.get("slurm", {}).get("queue") or "core"
    num_cores = config.get("slurm", {}).get("cores") or 16
    slurm_time = config.get("piper", {}).get("job_walltime", {}).get(workflow_name) or "4-00:00:00"
    slurm_out_log = os.path.join(perm_analysis_dir, "logs", "{}_sbatch.out".format(job_identifier))
    slurm_err_log = os.path.join(perm_analysis_dir, "logs", "{}_sbatch.err".format(job_identifier))
    for log_file in slurm_out_log, slurm_err_log:
        filesystem.rotate_file(log_file)
    sbatch_text = sample_runs_handler.create_sbatch_header(slurm_project_id=slurm_project_id,
                                       slurm_queue=slurm_queue,
                                       num_cores=num_cores,
                                       slurm_time=slurm_time,
                                       job_name="piper_{}".format(job_identifier),
                                       slurm_out_log=slurm_out_log,
                                       slurm_err_log=slurm_err_log)
    sbatch_text_list = sbatch_text.split("\n")
    sbatch_extra_params = config.get("slurm", {}).get("extra_params", {})
    for param, value in sbatch_extra_params.iteritems():
        sbatch_text_list.append("#SBATCH {} {}\n\n".format(param, value))
    modules_to_load = config.get("piper", {}).get("load_modules", [])
    if modules_to_load:
        sbatch_text_list.append("\n# Load required modules for Piper")
        for module_name in modules_to_load:
            sbatch_text_list.append("module load {}".format(module_name))

    # Fastq files to copy
    fastq_src_dst_list = []
    directories_to_create = set()
    for libprep in sample:
        for seqrun in libprep:
            project_specific_path = os.path.join(project.dirname,
                                                     sample.dirname,
                                                     libprep.dirname,
                                                     seqrun.dirname)
            directories_to_create.add(os.path.join("$SNIC_TMP/DATA/",
                                                       project_specific_path))
            for fastq in seqrun.fastq_files:
                src_file = os.path.join(project.base_path, "DATA",
                                            project_specific_path, fastq)
                dst_file = os.path.join("$SNIC_TMP/DATA/",
                                            project_specific_path,
                                            fastq)
                fastq_src_dst_list.append([src_file, dst_file])

    sbatch_text_list.append("echo -ne '\\n\\nCopying fastq files at '")
    sbatch_text_list.append("date")
    if fastq_src_dst_list:
        for directory in directories_to_create:
            sbatch_text_list.append("mkdir -p {}".format(directory))
        for src_file, dst_file in fastq_src_dst_list:
            sbatch_text_list.append("rsync -rptoDLv {} {}".format(src_file, dst_file))
    else:
        raise ValueError(('No valid fastq files available to process for '
                          'project/sample {}/{}'.format(project, sample)))

    # Pre-existing analysis files
    if files_to_copy:
        sbatch_text_list.append("echo -ne '\\n\\nCopying pre-existing analysis files at '")
        sbatch_text_list.append("date")

        sbatch_text_list.append("if [ ! -d {output directory} ]; then")
        sbatch_text_list.append("mkdir {output directory} ")
        sbatch_text_list.append("fi")
        sbatch_text_list.append(("rsync -rptoDLv {input_files} "
                                 "{output_directory}/").format(input_files=" ".join(files_to_copy),
                                                               output_directory=scratch_analysis_dir))
        # Delete pre-existing analysis files after copy
        sbatch_text_list.append("echo -ne '\\n\\nDeleting pre-existing analysis files at '")
        sbatch_text_list.append("date")
        sbatch_text_list.append("rm -rf {input_files}".format(input_files=" ".join(files_to_copy)))

    sbatch_text_list.append("echo -ne '\\n\\nExecuting command lines at '")
    sbatch_text_list.append("date")
    sbatch_text_list.append("# Run the actual commands")
    for command_line in command_line_list:
        sbatch_text_list.append(command_line)


    piper_status_file = sample_runs_handler.create_exit_code_file_path(workflow_subtask=workflow_name,
                                                   project_base_path=project.base_path,
                                                   project_name=project.dirname,
                                                   project_id=project.project_id,
                                                   sample_id=sample.name)
    sbatch_text_list.append("\nPIPER_RETURN_CODE=$?")

    #Precalcuate md5sums
    sbatch_text_list.append('MD5FILES="$SNIC_TMP/ANALYSIS/{}/piper_ngi/05_processed_alignments/*.bam'.format(project.project_id))
    sbatch_text_list.append('$SNIC_TMP/ANALYSIS/{}/piper_ngi/05_processed_alignments/*.table'.format(project.project_id))
    sbatch_text_list.append('$SNIC_TMP/ANALYSIS/{}/piper_ngi/07_variant_calls/*.genomic.vcf.gz'.format(project.project_id))
    sbatch_text_list.append('$SNIC_TMP/ANALYSIS/{}/piper_ngi/07_variant_calls/*.annotated.vcf.gz"'.format(project.project_id))
    sbatch_text_list.append('for f in $MD5FILES')
    sbatch_text_list.append('do')
    sbatch_text_list.append("    md5sum $f | awk '{printf $1}' > $f.md5 &")
    sbatch_text_list.append('done')
    sbatch_text_list.append('wait')
    
    #Copying back files
    sbatch_text_list.append("echo -ne '\\n\\nCopying back the resulting analysis files at '")
    sbatch_text_list.append("date")
    sbatch_text_list.append("mkdir -p {}".format(perm_analysis_dir))
    sbatch_text_list.append("rsync -rptoDLv {}/ {}/".format(scratch_analysis_dir, perm_analysis_dir))
    sbatch_text_list.append("\nRSYNC_RETURN_CODE=$?")

    # Record job completion status
    sbatch_text_list.append("if [[ $RSYNC_RETURN_CODE == 0 ]]")
    sbatch_text_list.append("then")
    sbatch_text_list.append("  if [[ $PIPER_RETURN_CODE == 0 ]]")
    sbatch_text_list.append("  then")
    sbatch_text_list.append("    echo '0'> {}".format(piper_status_file))
    sbatch_text_list.append("  else")
    sbatch_text_list.append("    echo '1'> {}".format(piper_status_file))
    sbatch_text_list.append("  fi")
    sbatch_text_list.append("else")
    sbatch_text_list.append("  echo '2'> {}".format(piper_status_file))
    sbatch_text_list.append("fi")

    # Write the sbatch file
    sbatch_dir = os.path.join(perm_analysis_dir, "sbatch")
    filesystem.safe_makedir(sbatch_dir)
    sbatch_outfile = os.path.join(sbatch_dir, "{}.sbatch".format(job_identifier))
    filesystem.rotate_file(sbatch_outfile)
    with open(sbatch_outfile, 'w') as f:
        f.write("\n".join(sbatch_text_list))
    LOG.info("Queueing sbatch file {} for job {}".format(sbatch_outfile, job_identifier))
    # Queue the sbatch file
    p_handle = filesystem.execute_command_line("sbatch {}".format(sbatch_outfile),
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE)
    p_out, p_err = p_handle.communicate()
    try:
        slurm_job_id = re.match(r'Submitted batch job (\d+)', p_out).groups()[0]
    except AttributeError:
        raise RuntimeError('Could not submit sbatch job for workflow "{}": '
                           '{}'.format(job_identifier, p_err))
    # Detail which seqruns we've started analyzing so we can update statuses later
    sample_runs_handler.record_analysis_details(project, job_identifier)
    return int(slurm_job_id)
