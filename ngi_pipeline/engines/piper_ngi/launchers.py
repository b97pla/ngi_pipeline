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
        self.project_analysis_directory = os.path.join(self.project_obj.base_path, "ANALYSIS", self.project_obj.dirname)
        self.sample_data_directory = os.path.join(self.project_data_directory, self.sample_obj.dirname)
        self.piper_analysis_directory = os.path.join(self.project_analysis_directory, "piper_ngi")

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
            return None

        # assert that the sample is ok to be started (this will raise an exception if that's not the case)
        try:
            self.assert_sample_should_be_started()
        except RuntimeError:
            raise

        submitted_job_ids = []

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

            # get an updated list of files that will go into the analysis from the local project object
            files_for_analysis = self.get_files_from_project_obj(local_project_obj)

            # create the setup command and get the setup xml output path
            setup_command, setup_xml_path = self.create_setup_command(local_project_obj, workflow_subtask)

            # create the piper command
            piper_command = self.create_piper_command(
                local_project_obj, workflow_subtask, setup_xml_path, exit_code_path)

            # submit the commands for execution
            submitted_job_id = self.submit_commands(
                [setup_command, piper_command],
                workflow_subtask,
                files_for_analysis,
                old_files_for_analysis,
                log_file_path,
                exit_code_path)

            submitted_job_ids.append(submitted_job_id)

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

        return submitted_job_ids

    def assert_exec_mode(self):
        implemented_modes = ["sbatch"]
        if self.exec_mode not in implemented_modes:
            raise NotImplementedError("{} execution mode not implemented. Only {} available.".format(
                self.exec_mode, ",".join(implemented_modes)))

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
                lambda f: NGIChipGenotypes(name=os.path.basename(f), dirname=os.path.basename(os.path.dirname(f))),
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
            self.config)

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

    def generate_sbatch_script(
            self,
            command_list,
            workflow_subtask,
            files_needed_for_analysis,
            old_files_needed_for_analysis,
            log_file_path,
            exit_code_path):

        # assert that we have a project id for slurm
        try:
            slurm_project_id = self.config["environment"]["project_id"]
        except KeyError:
            raise RuntimeError(
                "No SLURM project id specified in configuration file for project/sample {}/{}".format(
                    self.project_obj, self.sample_obj))

        # Paths to the various data directories
        scratch_analysis_dir = self.piper_analysis_directory.replace(self.project_obj.base_path, "$SNIC_TMP")
        scratch_data_dir = self.project_data_directory.replace(self.project_obj.base_path, "$SNIC_TMP")

        # get a list of source and destination file pairs to rsync
        files_to_rsync = map(
            lambda x: [os.path.join(self.project_data_directory, x), os.path.join(scratch_data_dir, x)],
            files_needed_for_analysis)

        # transform the list of files to statements
        rsync_fastq_file_statements = list(set(
            map(
                lambda x: "mkdir -p {}".format(os.path.dirname(x[1])),
                files_to_rsync)))
        rsync_fastq_file_statements.extend(
            map(
                lambda x: "rsync -rptoDLv {} {}".format(x[0], x[1]),
                files_to_rsync))

        # set up syncing and removal of previous analysis results, if applicable
        rsync_preexisting_results_statement = """
echo -ne "

Copying pre-existing analysis files at $(date)"
mkdir -p {scratch_analysis_dir}
rsync -rptoDLv {files_on_source} {scratch_analysis_dir}
echo -ne "

Deleting pre-existing analysis files at $(date)"
rm -rf {files_on_source}
        """.format(**{
            "scratch_analysis_dir": scratch_analysis_dir,
            "files_on_source": " ".join(old_files_needed_for_analysis)}) if old_files_needed_for_analysis else ""

        # calculate checksums of the large files
        checksum_calculation_statements = """
for f in $(find {scratch_analysis_dir} -type f -size +100M)
do
  md5sum $f |awk '{{printf $1}}' > $f.md5 &
done
wait
        """.format(scratch_analysis_dir=scratch_analysis_dir)

        # set up syncing the analysis folder back
        rsync_analysis_results_statements = """
rsync -rptoDLv {scratch_analysis_dir}{sep} {local_analysis_dir}{sep}
        """.format(**{
            "local_analysis_dir": self.piper_analysis_directory,
            "scratch_analysis_dir": scratch_analysis_dir,
            "sep": os.path.sep})

        # prepare the dictionary to use for replacing placeholders in the script template
        sbatch_parameters = {
            "slurm_project_id": slurm_project_id,
            "slurm_queue": self.config.get("slurm", {}).get("queue", "core"),
            "num_cores": self.config.get("slurm", {}).get("cores", "16"),
            "num_nodes": self.config.get("slurm", {}).get("nodes", "1"),
            "slurm_time": self.config.get("piper", {}).get("job_walltime", {}).get(workflow_subtask, "10-00:00:00"),
            "job_name": "piper_{}".format(self.job_identifier(workflow_subtask)),
            "sbatch_extra_params": "\n".join(
                map(
                    lambda (k, v): "#SBATCH {} {}".format(k, v),
                    self.config.get("slurm", {}).get("extra_params", {}).items())),
            "load_module_statements": "\n".join(
                map(
                    lambda x: "module load {}".format(x),
                    self.config.get("piper", {}).get("load_modules", []))),
            "rsync_fastq_file_statements": "\n".join(rsync_fastq_file_statements),
            "rsync_preexisting_results_statement": rsync_preexisting_results_statement,
            "command_line_statements": "\n".join(command_list),
            "checksum_calculation_statements": checksum_calculation_statements,
            "rsync_analysis_results_statements": rsync_analysis_results_statements,
            "piper_status_file": exit_code_path
        }

        # rotate the log files
        for log_stream in ["out", "err"]:
            log_file = log_file_path.replace(".log", "_sbatch.{}".format(log_stream))
            self.rotate_log_file(log_file)
            sbatch_parameters["slurm_{}_log".format(log_stream)] = log_file

        return sbatch_script_template().format(**sbatch_parameters)

    @staticmethod
    def get_files_from_project_obj(local_project_obj):
        project_files = [chip_genotype.name for chip_genotype in local_project_obj.chip_genotypes or []]
        for local_sample_obj in local_project_obj:
            project_files.extend(PiperLauncher.get_files_from_sample_obj(local_sample_obj))
        return project_files

    @staticmethod
    def get_files_from_sample_obj(local_sample_obj):
        sample_files = []
        for libprep in local_sample_obj:
            for seqrun in libprep:
                subdir = os.path.join(
                    local_sample_obj.dirname,
                    libprep.dirname,
                    seqrun.dirname)
                for fqfile in seqrun.fastq_files:
                    sample_files.append(
                        os.path.join(subdir, fqfile))
        return sample_files

    def job_identifier(self, workflow_subtask):
        return "-".join([self.project_obj.name, self.sample_obj.name, workflow_subtask])

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
            lambda fq: not valid_seqrun_paths or any(
                map(
                    lambda pth: os.path.split(fq)[0].endswith(pth),
                    valid_seqrun_paths)),
            fastq_files)

    def locate_preexisting_data_to_include_in_analysis(self, include_genotype_files=True):
        return self.preexisting_sample_runs_handler.find_previous_sample_analyses(
            self.project_obj, self.sample_obj,  include_genotype_files=include_genotype_files)

    def record_analysis_details(self, workflow_subtask):
        self.sample_runs_handler.record_analysis_details(self.project_obj, self.job_identifier(workflow_subtask))

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

    def submit_sbatch_commands(
            self,
            command_list,
            workflow_subtask,
            files_needed_for_analysis,
            old_files_needed_for_analysis,
            log_file_path,
            exit_code_path):

        # generate the sbatch script for the analysis
        sbatch_script = self.generate_sbatch_script(
            command_list,
            workflow_subtask,
            files_needed_for_analysis,
            old_files_needed_for_analysis,
            log_file_path,
            exit_code_path)

        # write the sbatch script to disk
        sbatch_script_file = self.write_sbatch_script(workflow_subtask, sbatch_script)

        job_identifier = self.job_identifier(workflow_subtask)
        LOG.info("Queueing sbatch file {} for job {}".format(sbatch_script_file, job_identifier))

        # queue the sbatch file
        p_handle = self.filesystem_handler.execute_command_line(
            "sbatch {}".format(sbatch_script_file),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
        p_out, p_err = p_handle.communicate()

        # parse the slurm job_id from stdout
        try:
            slurm_job_id = re.match(r'Submitted batch job (\d+)', p_out).groups()[0]
        except AttributeError:
            raise RuntimeError(
                "Could not submit sbatch job for workflow \"{}\": {}".format(job_identifier, p_err))

        # write details which seqruns we've started analyzing so we can update statuses later
        self.record_analysis_details(workflow_subtask)

        return {"slurm_job_id": int(slurm_job_id)}

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

    def write_sbatch_script(self, workflow_subtask, sbatch_script):
        sbatch_outfile = os.path.join(
            self.piper_analysis_directory, "sbatch", "{}.sbatch".format(
                self.job_identifier(workflow_subtask)))
        self.filesystem_handler.safe_makedir(os.path.dirname(sbatch_outfile))
        self.rotate_log_file(sbatch_outfile)
        with open(sbatch_outfile, 'w') as f:
            f.write(sbatch_script)
        return sbatch_outfile


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


def sbatch_script_template():
    return """
#!/bin/bash -l

#SBATCH -A {slurm_project_id}
#SBATCH -p {slurm_queue}
#SBATCH -n {num_cores}
#SBATCH -N {num_nodes}
#SBATCH -t {slurm_time}
#SBATCH -J {job_name}
#SBATCH -o {slurm_out_log}
#SBATCH -e {slurm_err_log}
{sbatch_extra_params}

# Load required modules for Piper
{load_module_statements}

echo -ne "

Copying fastq files at $(date)"

{rsync_fastq_file_statements}

{rsync_preexisting_results_statement}

echo -ne "

Executing command lines at $(date)"

# Run the actual commands
{command_line_statements}

PIPER_RETURN_CODE=$?

echo -ne "

Calculating checksums for selected files at $(date)"
{checksum_calculation_statements}

echo -ne "

Copying back the resulting analysis files at $(date)"

{rsync_analysis_results_statements}

RSYNC_RETURN_CODE=$?

# Record job completion status
if [[ $RSYNC_RETURN_CODE == 0 ]]
then
    if [[ $PIPER_RETURN_CODE == 0 ]]
    then
        echo 0 > {piper_status_file}
    else
        echo 1 > {piper_status_file}
    fi
else
    echo 2 > {piper_status_file}
fi
    """