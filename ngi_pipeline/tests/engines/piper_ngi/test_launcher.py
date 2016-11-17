import mock
import unittest
import os

from ngi_pipeline.engines.piper_ngi.launchers import PiperLauncher
from ngi_pipeline.tests.engines.piper_ngi.test_utils import ProjectData
import ngi_pipeline.engines.piper_ngi.workflows
import ngi_pipeline.engines.piper_ngi.local_process_tracking
import ngi_pipeline.engines.piper_ngi.utils
import ngi_pipeline.utils.filesystem
import ngi_pipeline.utils.slurm
import ngi_pipeline.engines.piper_ngi.command_creation_config


class TestPiperLauncher(unittest.TestCase):

    def setUp(self):

        # get an example project object
        self.project_data = ProjectData()
        project_obj = next(self.project_data.next())

        # inject mocks to mask the external calls
        self.mocks = {
            "workflow_handler": mock.create_autospec(
                ngi_pipeline.engines.piper_ngi.workflows, spec_set=True),
            "local_process_tracking_handler": mock.create_autospec(
                ngi_pipeline.engines.piper_ngi.local_process_tracking, spec_set=True),
            "preexisting_sample_runs_handler": mock.create_autospec(
                ngi_pipeline.engines.piper_ngi.utils, spec_set=True),
            "sample_runs_handler": mock.create_autospec(
                ngi_pipeline.engines.piper_ngi.utils, spec_set=True),
            "filesystem_handler": mock.create_autospec(
                ngi_pipeline.utils.filesystem, spec_set=True),
            "command_handler": mock.create_autospec(
                ngi_pipeline.engines.piper_ngi.command_creation_config, spec_set=True),
            "launch_handler": mock.create_autospec(
                ngi_pipeline.engines.piper_ngi.launchers, spec_set=True),
            "slurm_handler": mock.create_autospec(
                ngi_pipeline.utils.slurm, spec_set=True)
        }

        # create a test instance and inject mocks
        self.launcher = PiperLauncher(
            project_obj,
            project_obj.samples.values()[0],
            **self.mocks
        )

    def test_analyze(self):
        # create a mock from the test instance of PiperLauncher
        mock_launcher = mock.create_autospec(self.launcher, spec_set=True, instance=True)
        # set the attributes on the mock according to the test instance
        mock_launcher.configure_mock(**self.launcher.__dict__)
        # set up the return values of the instance methods according to the test scenario
        workflow_subtasks = ["test_subtask_1", "test_subtask_2"]
        libprep_seqruns = {
            "libprep1": ["seqrun_11", "seqrun_12"],
            "libprep2": ["seqrun_21", "seqrun_22"]}
        chip_genotypes_for_analysis = [
            os.path.join("path", "to", "genotype", "file1"),
            os.path.join("path", "to", "genotype", "file2")]
        fastq_files_for_analysis = [
            os.path.join("path", "to", "libprep1", "seqrun11", "fastq_file1"),
            os.path.join("path", "to", "libprep1", "seqrun11", "fastq_file2"),
            os.path.join("path", "to", "libprep1", "seqrun12", "fastq_file1"),
            os.path.join("path", "to", "libprep1", "seqrun12", "fastq_file2"),
            os.path.join("path", "to", "libprep2", "seqrun21", "fastq_file1"),
            os.path.join("path", "to", "libprep2", "seqrun21", "fastq_file2"),
            os.path.join("path", "to", "libprep2", "seqrun22", "fastq_file1"),
            os.path.join("path", "to", "libprep2", "seqrun22", "fastq_file2")]
        old_files_for_analysis = [
            os.path.join("path", "to", "old", "file1"),
            os.path.join("path", "to", "old", "file2")]
        exit_code_path = os.path.join("path", "to", "exit", "code", "file")
        log_path = os.path.join("path", "to", "log", "file")
        local_project_obj = self.launcher.project_obj
        piper_setup_cmd = "this is a test piper setup command"
        piper_setup_xml_path = os.path.join("path", "to", "piper", "setup", "xml")
        piper_cmd = "this is a test piper command"
        submitted_job_ids = ({"slurm_job_id": 1}, {"slurm_job_id": 2})
        mock_launcher.determine_workflow_subtasks.return_value = workflow_subtasks
        mock_launcher.determine_seqruns_to_be_analyzed.return_value = libprep_seqruns
        mock_launcher.create_exit_code_path.return_value = exit_code_path
        mock_launcher.create_log_file_path.return_value = log_path
        mock_launcher.locate_chip_genotype_files_for_project.return_value = chip_genotypes_for_analysis
        mock_launcher.locate_fastq_files_for_project_sample.return_value = fastq_files_for_analysis
        mock_launcher.locate_preexisting_data_to_include_in_analysis.return_value = old_files_for_analysis
        mock_launcher.create_local_copy_of_project_obj.return_value = local_project_obj
        mock_launcher.create_setup_command.return_value = (piper_setup_cmd, piper_setup_xml_path)
        mock_launcher.create_piper_command.return_value = piper_cmd
        mock_launcher.submit_commands.side_effect = submitted_job_ids

        # have the analyze method in the mocked instance call a wrapper with the actual code as a side effect
        def _analyze_wrapper():
            return PiperLauncher.analyze(mock_launcher)

        mock_launcher.analyze.side_effect = _analyze_wrapper

        # if exec_mode assertion fails, analyze method should return None
        mock_launcher.assert_exec_mode.side_effect = NotImplementedError
        self.assertIsNone(mock_launcher.analyze())
        mock_launcher.assert_exec_mode.side_effect = None

        # if a sample should not be started, a RuntimeError should be raised
        mock_launcher.assert_sample_should_be_started.side_effect = RuntimeError
        with self.assertRaises(RuntimeError):
            mock_launcher.analyze()
        mock_launcher.assert_sample_should_be_started.side_effect = None

        # if a subtask is already running, nothing should be done and an empty list returned
        mock_launcher.assert_sample_is_not_running.return_value = False
        self.assertListEqual([], mock_launcher.analyze())
        self.assertListEqual(
            map(mock.call, workflow_subtasks),
            mock_launcher.assert_sample_is_not_running.call_args_list)
        mock_launcher.clean_previous_analysis_results.assert_not_called()
        mock_launcher.assert_sample_is_not_running.reset_mock()

        # if no valid seqruns are found for a subtask, it should be skipped
        mock_launcher.assert_sample_is_not_running.return_value = True
        mock_launcher.determine_seqruns_to_be_analyzed.side_effect = RuntimeError
        self.assertListEqual([], mock_launcher.analyze())
        mock_launcher.locate_chip_genotype_files_for_project.assert_not_called()
        mock_launcher.determine_seqruns_to_be_analyzed.side_effect = None

        # a raised exception when locating files should result in a subtask being skipped
        for mocked_method in [
            mock_launcher.locate_chip_genotype_files_for_project,
            mock_launcher.locate_fastq_files_for_project_sample,
            mock_launcher.locate_preexisting_data_to_include_in_analysis
        ]:
            mocked_method.side_effect = ValueError
            self.assertListEqual([], mock_launcher.analyze())
            mock_launcher.create_exit_code_path.assert_not_called()
            mocked_method.side_effect = None

        # now, executing the analyze method should return job ids for the workflows
        mock_launcher.locate_fastq_files_for_project_sample.reset_mock()
        self.assertListEqual(list(submitted_job_ids), mock_launcher.analyze())
        # do a bunch of assertions on the method calls
        self.assertListEqual(
            [mock.call(valid_seqruns=libprep_seqruns) for i in xrange(2)],
            mock_launcher.locate_fastq_files_for_project_sample.call_args_list)
        for mocked_method in [
            mock_launcher.create_exit_code_path,
            mock_launcher.create_log_file_path]:
            self.assertListEqual(
                map(mock.call, workflow_subtasks),
                mocked_method.call_args_list)
        mock_launcher.rotate_log_file.assert_called_with(log_path)
        mock_launcher.create_local_copy_of_project_obj.assert_called_with(
            fastq_files_for_analysis, chip_genotypes_for_analysis)
        self.assertListEqual(
            map(lambda x: mock.call(local_project_obj, x),
                workflow_subtasks),
            mock_launcher.create_setup_command.call_args_list)
        self.assertListEqual(
            map(lambda x: mock.call(local_project_obj, x, piper_setup_xml_path, exit_code_path),
                workflow_subtasks),
            mock_launcher.create_piper_command.call_args_list)
        self.assertListEqual(
            map(lambda x: mock.call([piper_setup_cmd, piper_cmd], x, old_files_for_analysis),
                workflow_subtasks),
            mock_launcher.submit_commands.call_args_list)
        self.assertListEqual(
            map(mock.call, submitted_job_ids),
            mock_launcher.submitted_job_status.call_args_list)
        self.assertListEqual([
            mock.call(
                workflow_subtasks[i],
                list(submitted_job_ids)[i]) for i in xrange(2)],
            mock_launcher.record_job_status_to_databases.call_args_list)
        mock_launcher.submit_commands.side_effect = submitted_job_ids

        # failure to verify the status of the submitted job should not affect the execution
        mock_launcher.submitted_job_status.side_effect = RuntimeError
        self.assertListEqual(list(submitted_job_ids), mock_launcher.analyze())
        mock_launcher.submitted_job_status.side_effect = None
        mock_launcher.submit_commands.side_effect = submitted_job_ids

        # failure to record the status of the submitted job should not affect the execution
        for eff in (RuntimeError, ValueError):
            mock_launcher.record_job_status_to_databases.side_effect = eff
            self.assertListEqual(list(submitted_job_ids), mock_launcher.analyze())
            mock_launcher.submit_commands.side_effect = submitted_job_ids
        mock_launcher.record_job_status_to_databases.side_effect = None

    def test_assert_exec_mode(self):
        # assert with valid exec mode
        try:
            self.launcher.assert_exec_mode()
        except Exception as e:
            self.fail("should not have raised {}".format(e))

        self.launcher.exec_mode = "not-sbatch"
        with self.assertRaises(NotImplementedError):
            self.launcher.assert_exec_mode()

    def test_assert_sample_should_be_started(self):
        # not raising an exception should mean that the sample can be started
        try:
            self.launcher.assert_sample_should_be_started()
        except Exception as e:
            self.fail("should not have raised {}".format(e))

        # a raised exception should be passed on
        self.launcher.preexisting_sample_runs_handler.check_for_preexisting_sample_runs.side_effect = RuntimeError
        with self.assertRaises(RuntimeError):
            self.launcher.assert_sample_should_be_started()

    def test_assert_sample_is_not_running(self):
        # assert correct result depending on whether sample is running
        for is_running in (True, False):
            self.launcher.local_process_tracking_handler.is_sample_analysis_running_local.return_value = is_running
            # the result should be independent of whether running analysis should be restarted
            for restart_running in (True, False):
                self.launcher.restart_running_jobs = restart_running
                self.assertEqual(not is_running, self.launcher.assert_sample_is_not_running("subtask"))

    def test_clean_previous_analysis_results(self):
        # if previous results should be kept, nothing should happen
        self.launcher.keep_existing_data = True
        self.launcher.clean_previous_analysis_results()
        self.launcher.preexisting_sample_runs_handler.remove_previous_sample_analyses.assert_not_called()
        self.launcher.preexisting_sample_runs_handler.remove_previous_genotype_analyses.assert_not_called()
        # if data should not be kept, the corresponding cleanup routines should be called
        expected_calls = {
            "sample": [
                self.launcher.preexisting_sample_runs_handler.remove_previous_sample_analyses,
                self.launcher.preexisting_sample_runs_handler.remove_previous_genotype_analyses],
            "genotype": [
                self.launcher.preexisting_sample_runs_handler.remove_previous_genotype_analyses,
                self.launcher.preexisting_sample_runs_handler.remove_previous_sample_analyses]}
        self.launcher.keep_existing_data = False
        for level, mocks in expected_calls.items():
            map(lambda x: x.reset_mock(), mocks)
            self.launcher.level = level
            self.launcher.clean_previous_analysis_results()
            self.assertEqual(1, mocks[0].call_count)
            mocks[1].assert_not_called()

    def test_create_local_copy_of_project_obj(self):
        # spoof a list of fastq files and genotype files
        fastq_files = []
        genotype_files = []
        for i in xrange(2):
            genotype_files.append(os.path.join("this", "is", "a", "path", "to", "genotype_file_{}.vcf".format(i)))
            genotype_files.append(os.path.join("this", "is", "a", "path", "to", "genotype_file_{}.vcf.gz".format(i)))
            libprep = "libprep{}".format(i)
            for j in xrange(3):
                seqrun = "seqrun{}".format(j)
                for k in xrange(4):
                    fastq_files.append(
                        os.path.join(
                            "this", "is", "path", "to", libprep, seqrun, "{}_{}_fastq{}.fastq.gz".format(
                                libprep, seqrun, k)))

        # create a copy of the original object with the spoofed files attached
        local_project_obj = self.launcher.create_local_copy_of_project_obj(fastq_files, genotype_files)

        # assert genotype files were attached
        self.assertItemsEqual(
            map(os.path.basename, genotype_files),
            map(lambda x: x.name, local_project_obj.chip_genotypes))
        # assert attributes were copied
        for attr in ("name", "dirname", "project_id", "base_path"):
            self.assertEqual(getattr(self.launcher.project_obj, attr), getattr(local_project_obj, attr))

        # assert attributes on sample were copied
        local_sample_obj = next(iter(local_project_obj))
        for attr in ("name", "dirname"):
            self.assertEqual(getattr(self.launcher.sample_obj, attr), getattr(local_sample_obj, attr))
        fastq_names = []
        for libprep_obj in local_sample_obj:
            for seqrun_obj in libprep_obj:
                for fastq_name in seqrun_obj:
                    # assert that libprep and seqrun folder names were parsed correctly and
                    # the corresponding fastq file attached
                    self.assertTrue(fastq_name.startswith("{}_{}_fastq".format(libprep_obj.name, seqrun_obj.name)))
                    fastq_names.append(fastq_name)
        # assert that all fastq files were attached
        self.assertItemsEqual(map(os.path.basename, fastq_files), fastq_names)

    def test_determine_seqruns_to_be_analyzed(self):
        expected_seqruns = {"this": "is", "a": "dict", "of": "\"seqruns\""}
        self.launcher.preexisting_sample_runs_handler.get_valid_seqruns_for_sample.return_value = expected_seqruns
        self.assertDictEqual(expected_seqruns, self.launcher.determine_seqruns_to_be_analyzed())
        self.launcher.preexisting_sample_runs_handler.get_valid_seqruns_for_sample.return_value = {}
        with self.assertRaises(ValueError):
            self.launcher.determine_seqruns_to_be_analyzed()

    def test_locate_chip_genotype_files_for_project(self):
        expected_genotype_files = ["this", "is", "a", "bunch", "of", "genotype", "files"]
        self.launcher.filesystem_handler.locate_chip_genotypes_in_dir.return_value = expected_genotype_files
        self.assertListEqual(expected_genotype_files, self.launcher.locate_chip_genotype_files_for_project())
        self.launcher.filesystem_handler.locate_chip_genotypes_in_dir.return_value = []
        self.assertIsNone(self.launcher.locate_chip_genotype_files_for_project())

    def test_locate_fastq_files_for_project_sample(self):
        self.launcher.filesystem_handler.fastq_files_under_dir.return_value = []
        with self.assertRaises(ValueError):
            self.launcher.locate_fastq_files_for_project_sample()

        # spoof some valid and invalid seqruns with corresponding fastq files
        valid_seqruns = {}
        valid_fastq_files = []
        invalid_fastq_files = []
        for i in xrange(4):
            libprep = "libprep{}".format(i)
            if i < 3:
                valid_seqruns[libprep] = []
            for j in xrange(5):
                seqrun = "seqrun{}".format(j)
                fastq_files = []
                for k in xrange(3):
                    fastq_files.append(
                        os.path.join("path", "to", libprep, seqrun, "{}_{}_fastq{}.fastq.gz".format(
                            libprep, seqrun, k)))
                if i in [1, 2] and j < 3:
                    valid_seqruns[libprep].append(seqrun)
                    valid_fastq_files.extend(fastq_files)
                else:
                    invalid_fastq_files.extend(fastq_files)

        # assert that fastq files are not filtered if no valid seqruns are passed
        self.launcher.filesystem_handler.fastq_files_under_dir.return_value = valid_fastq_files
        self.assertListEqual(valid_fastq_files, self.launcher.locate_fastq_files_for_project_sample())

        # assert that fastq files are filtered against valid seqruns
        self.launcher.filesystem_handler.fastq_files_under_dir.return_value = valid_fastq_files + invalid_fastq_files
        self.assertListEqual(
            valid_fastq_files,
            self.launcher.locate_fastq_files_for_project_sample(
                valid_seqruns=valid_seqruns))

    def test_locate_preexisting_data_to_include_in_analysis(self):
        expected_analysis_results = ["this", "is", "some", "data", "files"]
        self.launcher.preexisting_sample_runs_handler.find_previous_sample_analyses.return_value = \
            expected_analysis_results
        self.assertListEqual(
            expected_analysis_results,
            self.launcher.locate_preexisting_data_to_include_in_analysis())

    def test_record_job_status_to_databases(self):
        expected_kwargs = {
            "project": self.launcher.project_obj,
            "sample": self.launcher.sample_obj,
            "analysis_module_name": "piper_ngi",
            "workflow_subtask": "workflow_subtask",
            "config": self.launcher.config,
            "slurm_job_id": 123456
        }
        for effect in (ValueError, RuntimeError):
            self.launcher.local_process_tracking_handler.record_process_sample.reset_mock()
            self.launcher.local_process_tracking_handler.record_process_sample.side_effect = effect
            with self.assertRaises(effect):
                self.launcher.record_job_status_to_databases(
                    expected_kwargs["workflow_subtask"],
                    {"slurm_job_id": expected_kwargs["slurm_job_id"]})
            self.launcher.local_process_tracking_handler.record_process_sample.assert_called_once_with(
                **expected_kwargs)

    def test_rotate_log_file(self):
        expected_log_file = "some expected log file path"
        self.launcher.rotate_log_file(expected_log_file)
        self.launcher.filesystem_handler.rotate_file.assert_called_once_with(expected_log_file)
        self.launcher.filesystem_handler.rotate_file.side_effect = OSError
        with self.assertRaises(OSError):
            self.launcher.rotate_log_file(expected_log_file)

    def _submit_sbatch_commands_helper(self, test_fn):
        expected_args = [
            ["expected command 1",
             "expected command 2"],
            "workflow_subtask",
            self.launcher.project_obj,
            self.launcher.sample_obj]
        expected_kwargs = {
            "restart_finished_jobs": self.launcher.restart_finished_jobs,
            "files_to_copy": [
                "file 1 needed for analysis",
                "file 2 needed for analysis",
                "file 3 needed for analysis"],
            "config": self.launcher.config}
        expected_slurm_job_id = {"slurm_job_id": 123456}
        self.launcher.launch_handler.sbatch_piper_sample.return_value = expected_slurm_job_id["slurm_job_id"]
        self.assertDictEqual(
            expected_slurm_job_id,
            test_fn(
                expected_args[0],
                expected_args[1],
                expected_kwargs["files_to_copy"]))
        self.launcher.launch_handler.sbatch_piper_sample.assert_called_once_with(*expected_args, **expected_kwargs)

    def test_submit_commands(self):
        # when exec_mode is "sbatch", this method will only be a wrapper around submit_sbatch_commands
        self.launcher.exec_mode = "sbatch"
        self._submit_sbatch_commands_helper(self.launcher.submit_commands)

    def test_submit_sbatch_commands(self):
        self._submit_sbatch_commands_helper(self.launcher.submit_sbatch_commands)

    def _submitted_sbatch_job_status_helper(self, test_fn):
        job_id = {"slurm_job_id": 123456}
        expected_job_status = "test status"

        # first just a straightforward call
        self.launcher.slurm_handler.get_slurm_job_status.return_value = expected_job_status
        self.assertEqual(expected_job_status, test_fn(job_id))

        # next test the case where no result is found
        # for speed, patch the time.sleep method
        with mock.patch("ngi_pipeline.engines.piper_ngi.launchers.time.sleep"):
            self.launcher.slurm_handler.get_slurm_job_status.side_effect = ValueError
            with self.assertRaises(RuntimeError):
                test_fn(job_id)

    def test_submitted_job_status(self):
        self.launcher.exec_mode = "sbatch"
        self._submitted_sbatch_job_status_helper(self.launcher.submitted_job_status)

    def test_submitted_sbatch_job_status(self):
        self._submitted_sbatch_job_status_helper(self.launcher.submitted_sbatch_job_status)
