import os

from ngi_pipeline.engines.sarek.models.resources import Runfolder, SampleFastq
from ngi_pipeline.utils.filesystem import is_index_file


class SarekAnalysisSample:

    def __init__(self, project_object, sample_object, analysis_object, restart_options=None):
        self.analysis_object = analysis_object
        self.sampleid = sample_object.name
        self.projectid = project_object.project_id
        self.project_base_path = project_object.base_path
        self.restart_options = {
            key: False for key in ["restart_failed_jobs", "restart_finished_jobs", "restart_running_jobs"]}
        self.restart_options.update(restart_options or {})
        self.sample_ngi_object = sample_object

    def sample_data_path(self):
        return self.analysis_object.sample_data_path(self.project_base_path, self.projectid, self.sampleid)

    def sample_analysis_path(self):
        return self.analysis_object.sample_analysis_path(self.project_base_path, self.projectid, self.sampleid)

    def sample_seqrun_path(self, libprepid, seqrunid):
        return os.path.join(self.sample_data_path(), libprepid, seqrunid)

    def sample_analysis_exit_code_path(self):
        return self.analysis_object.sample_analysis_exit_code_path(
            self.project_base_path, self.projectid, self.sampleid)

    def sample_analysis_tsv_file(self):
        return self.analysis_object.sample_analysis_tsv_file(self.project_base_path, self.projectid, self.sampleid)

    def _get_sample_librep(self, libprepid):
        return list(filter(lambda lp: lp.name == libprepid, self.sample_ngi_object)).pop()

    def sample_libprep_ids(self):
        return list(map(lambda x: x.name, self.sample_ngi_object))

    def libprep_seqrun_ids(self, libprepid):
        return list(map(lambda x: x.name, self._get_sample_librep(libprepid)))

    def libpreps_to_analyze(self):
        return filter(
            lambda lp: self.analysis_object.libprep_should_be_started(
                self.projectid, self.sampleid, lp.name),
            self.sample_ngi_object)

    def seqruns_to_analyze(self, libprep):
        return filter(
            lambda sr: self.analysis_object.seqrun_should_be_started(
                self.projectid, self.sampleid, libprep.name, sr.name, self.restart_options),
            libprep)

    def runid_and_fastq_files_for_sample(self):
        for libprep in self.libpreps_to_analyze():
            for libprep_fastq in self._get_runid_and_fastq_files_for_libprep(libprep):
                yield libprep_fastq

    def _get_runid_and_fastq_files_for_libprep(self, libprep):
        for seqrun in self.seqruns_to_analyze(libprep):
            for seqrun_fastq in self._get_runid_and_fastq_files_for_seqrun(libprep.name, seqrun):
                yield seqrun_fastq

    def _get_runid_and_fastq_files_for_seqrun(self, libprepid, seqrun):
        from ngi_pipeline.engines.sarek.models.sarek import SarekGermlineAnalysis
        # the runfolder is represented with a Runfolder object
        runfolder = Runfolder(self.sample_seqrun_path(libprepid, seqrun.name))
        # iterate over the fastq file pairs belonging to the seqrun, filter out files with index reads
        sample_fastq_objects = map(
            lambda f: SampleFastq(os.path.join(runfolder.path, f)),
            filter(lambda fq: not is_index_file(fq), seqrun.fastq_files))
        for sample_fastq_file_pair in SarekGermlineAnalysis.sample_fastq_file_pair(sample_fastq_objects):
            runid = "{}.{}.{}".format(
                runfolder.flowcell_id, sample_fastq_file_pair[0].lane_number, sample_fastq_file_pair[0].sample_number)
            yield [runid] + map(lambda x: x.path, sample_fastq_file_pair)
