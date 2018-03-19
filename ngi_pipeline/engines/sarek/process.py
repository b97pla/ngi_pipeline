import datetime
import errno
import os
import re
import subprocess

from ngi_pipeline.engines.sarek.exceptions import SlurmStatusNotRecognizedError
from ngi_pipeline.utils.filesystem import execute_command_line, safe_makedir, chdir
from ngi_pipeline.utils.slurm import get_slurm_job_status as core_get_slurm_job_status


class ProcessConnector(object):
    # connector for running processes, subclasses can run in other environments e.g. SLURM

    def __init__(self, cwd=None):
        self.cwd = cwd or os.curdir

    def execute_process(self, command_line, working_dir=None, **extra_args):
        working_dir = working_dir or self.cwd
        safe_makedir(working_dir)
        with chdir(working_dir):
            try:
                proc = execute_command_line(command_line, shell=False, cwd=working_dir)
                return proc.pid
            except RuntimeError as re:
                raise


class ProcessStatus(object):

    @classmethod
    def get_type_from_processid_and_exit_code_path(cls, processid, exit_code_path):
        return ProcessRunning if cls.is_process_running(processid) else \
            ProcessExitStatus.get_type_from_exit_code_path(exit_code_path)

    @staticmethod
    def is_process_running(processid):
        try:
            os.kill(processid, 0)
            return True
        except OSError as ose:
            if ose.errno != errno.ESRCH:
                raise
        return False


class ProcessRunning(ProcessStatus):
    pass


class ProcessStopped(ProcessStatus):
    pass


class JobStatus(ProcessStatus):

    @staticmethod
    def is_process_running(jobid):
        return SlurmConnector.get_slurm_job_status(jobid) == ProcessRunning


class ProcessExitStatus(ProcessStatus):

    @staticmethod
    def get_type_from_exit_code_path(exit_code_path):
        exit_code = ProcessExitStatus.get_exit_code(exit_code_path)
        if exit_code is None:
            return ProcessExitStatusUnknown
        if exit_code != 0:
            return ProcessExitStatusFailed
        return ProcessExitStatusSuccessful

    @staticmethod
    def get_exit_code(exit_code_path):
        try:
            with open(exit_code_path, "r") as fh:
                return int(fh.readline().strip())
        except IOError as e:
            if e.errno != errno.ENOENT:
                raise
        except (ValueError, TypeError) as e:
            pass
        return None


class ProcessExitStatusSuccessful(ProcessExitStatus):
    pass


class ProcessExitStatusFailed(ProcessExitStatus):
    pass


class ProcessExitStatusUnknown(ProcessExitStatus):
    pass


class SlurmConnector(ProcessConnector):
    # a connector for running SLURM jobs

    JOB_HEADER_TEMPLATE = """#! /bin/bash -l

#SBATCH -A {slurm_project}
#SBATCH -J "{slurm_job_name}"
#SBATCH -p {slurm_partition}
#SBATCH -N {slurm_nodes}
#SBATCH -n {slurm_cores}
#SBATCH -t {slurm_job_time}
#SBATCH -o "{slurm_stdout}"
#SBATCH -e "{slurm_stderr}"
#SBATCH -D "{slurm_working_directory}"
#SBATCH --mail-user {slurm_mail_user}
#SBATCH --mail-type {slurm_mail_events}
"""

    SLURM_DEFAULTS = {
        "slurm_partition": "core",
        "slurm_nodes": 1,
        "slurm_cores": 1,
        "slurm_job_time": "48:00:00",
        "slurm_mail_events": "NONE"
    }

    def __init__(self, slurm_project, slurm_mail_user, cwd=None, **slurm_args):
        super(SlurmConnector, self).__init__(cwd=cwd)
        slurm_args["slurm_project"] = slurm_project
        slurm_args["slurm_mail_user"] = slurm_mail_user
        self.slurm_parameters = self.SLURM_DEFAULTS.copy()
        self.slurm_parameters.update(slurm_args)

    def _slurm_script_from_command_line(
            self,
            command_line,
            working_dir,
            exit_code_path,
            job_name):
        # create the script in the passed working directory
        slurm_script = os.path.join(
            working_dir, "{}.{}.sbatch".format(job_name, datetime.datetime.now().strftime("%s")))
        slurm_stdout = "{}.out".format(slurm_script)
        slurm_stderr = "{}.err".format(slurm_script)

        with open(slurm_script, "w") as fh:
            fh.write(
                SlurmConnector.JOB_HEADER_TEMPLATE.format(
                    slurm_job_name=job_name,
                    slurm_stdout=slurm_stdout,
                    slurm_stderr=slurm_stderr,
                    slurm_working_directory=self.slurm_parameters.get("slurm_working_directory", working_dir),
                    **self.slurm_parameters))
            # append any extra SLURM arguments passed
            for slurm_extra_arg in self.slurm_parameters.get("slurm_extra_args", []):
                fh.write("#SBATCH {}\n".format(slurm_extra_arg))

            fh.write("\necho \"\" > \"{}\"\n".format(exit_code_path))
            fh.write("{}\n".format(command_line))
            fh.write("echo \"$?\" > \"{}\"\n".format(exit_code_path))

        return slurm_script

    def execute_process(self, command_line, working_dir=None, exit_code_path=None, job_name=None):
        # create the working dir if it does not exist already
        working_dir = working_dir or self.cwd
        exit_code_path = exit_code_path or os.devnull
        job_name = job_name or command_line.replace(" ", "_")[0:20]
        safe_makedir(working_dir)
        with chdir(working_dir):
            slurm_script = self._slurm_script_from_command_line(
                command_line,
                working_dir,
                exit_code_path,
                job_name)
            # submit the sbatch file
            sbatch_command_line = "sbatch {}".format(slurm_script)
            proc = execute_command_line(
                sbatch_command_line,
                shell=False,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE)
            stdout, stderr = proc.communicate()
            try:
                slurm_job_id = re.match(r'Submitted batch job (\d+)', stdout).groups()[0]
                return slurm_job_id
            except AttributeError:
                raise RuntimeError(
                    'Could not submit sbatch job for workflow "{}": {}'.format(job_name, stderr))

    @staticmethod
    def get_slurm_job_status(slurm_job_id):
        try:
            return ProcessRunning if core_get_slurm_job_status(slurm_job_id) is None else ProcessStopped
        except RuntimeError as e:
            raise SlurmStatusNotRecognizedError(slurm_job_id, e.message)
