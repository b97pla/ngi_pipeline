import mock
import unittest
import tempfile
import os
import shutil

from ngi_pipeline.conductor.classes import NGIProject
import ngi_pipeline.engines.piper_ngi.utils

class ProjectData():

    def __init__(self):
        self.proj_data = (
            {
                "proj_name": "Y.Mom_15_01",
                "proj_id": "P1155",
                "sequencing_facility": "NGI-S",
                "basepath": tempfile.mkdtemp(),
                "samples": [{
                    "libprep_id": "P1155_prepA",
                    "seqrun_id": "P1155_seqrunA",
                    "sample_name": "P1155_101"},
                {
                    "libprep_id": "P1155_prepAA",
                    "seqrun_id": "P1155_seqrunAA",
                    "sample_name": "P1155_1011"}]},
            {
                "proj_name": "Y.Mom_15_011",
                "proj_id": "P11551",
                "sequencing_facility": "NGI-S",
                "basepath": tempfile.mkdtemp(),
                "samples": [{
                    "libprep_id": "P11551_prepA",
                    "seqrun_id": "P11551_seqrunA",
                    "sample_name": "P11551_101"},
                {
                    "libprep_id": "P11551_prepAA",
                    "seqrun_id": "P11551_seqrunAA",
                    "sample_name": "P11551_1011"}]},
            {
                "proj_name": "OK-9999",
                "proj_id": "OK-9999",
                "sequencing_facility": "NGI-U",
                "basepath": tempfile.mkdtemp(),
                "samples": [{
                    "libprep_id": "OK-9999.SampleA.v1",
                    "seqrun_id": "161019_ST-1234_0123_ATEST123CXX",
                    "sample_name": "SampleA"},
                    {
                    "libprep_id": "OK-9999.SampleA.v11",
                    "seqrun_id": "161019_ST-1234_0123_ATEST123CXX",
                    "sample_name": "SampleAA"}]},
            {
                "proj_name": "OK-999",
                "proj_id": "OK-999",
                "sequencing_facility": "NGI-U",
                "basepath": tempfile.mkdtemp(),
                "samples": [{
                    "libprep_id": "OK-999.SampleA.v1",
                    "seqrun_id": "161019_ST-1234_0123_ATEST123CXX",
                    "sample_name": "SampleA"},
                    {
                    "libprep_id": "OK-999.SampleA.v11",
                    "seqrun_id": "161019_ST-1234_0123_ATEST123CXX",
                    "sample_name": "SampleAA"}]})
        self.piper_output_folders = {
            "01_this_is_raw_alignments": [
                ".raw.bam",
                ".raw.bam.bai"
            ],
            "02_this_is_preliminary_qc": [
                ".raw.qualimap.qc"
            ],
            "03_genotype_concordance": [
                ".genotype_concordance"
            ],
            "04_this_is_merged_bams": [
                ".bam"
                ".bam.bai"
            ],
            "05_this_is_final_bams": [
                ".clean.dedup.bam",
                ".clean.dedup.bam.bai",
                ".metrics"
            ],
            "06_this_is_final_qc": [
                ".clean.dedup.qualimap.qc"
            ],
            "07_this_is_genotype_calls": [
                ".clean.dedup.recal.raw.genomic.gvcf.gz",
                ".clean.dedup.recal.raw.genomic.raw.vcf.gz",
                ".clean.dedup.recal.raw.genomic.raw.vcf.gz.tbi"
            ],
            "08_this_is_miscallaneous_files": [],
            "logs": [
                "-merge.process.variantcall.sbatch.out",
                "-merge.process.variantcall.sbatch.err"
            ],
            "setup_xml_files": [
                "-merge.process.variantcall.setup.xml"
            ],
            "sbatch_files": [
                "-merge.process.variantcall.sbatch"
            ],
        }

    def next(self):
        for data in self.proj_data:
            proj_obj = NGIProject(
                name=data["proj_name"],
                dirname=data["proj_name"],
                project_id=data["proj_id"],
                base_path=data["basepath"])
            for sample_data in data["samples"]:
                sample = proj_obj.add_sample(
                    name=sample_data["sample_name"],
                    dirname=sample_data["sample_name"])
                libprep = sample.add_libprep(
                    name=sample_data["libprep_id"],
                    dirname=sample_data["libprep_id"]
                )
                seqrun = libprep.add_seqrun(
                    name=sample_data["seqrun_id"],
                    dirname=sample_data["seqrun_id"]
                )
            yield proj_obj

    def generate_sample_analysis_results(self, proj_obj, sample_obj):
        # generate file names corresponding to analysis results for a project object and its samples
        analysis_result_files = []
        for output_path, suffixes in self.piper_output_folders.items():
            folder = os.path.join(proj_obj.base_path, "ANALYSIS", proj_obj.project_id, "piper_ngi", output_path)
            for suffix in suffixes:
                for (pfx, sfx) in [("", ""), ("", ".out"), (".", ".done"), (".", ".fail")]:
                    analysis_result_files.append(
                        os.path.join(folder,
                                     "{}{}{}{}".format(pfx, sample_obj.name, suffix, sfx)))
        return analysis_result_files

    def create_sample_analysis_results(self, proj_obj, sample_obj, files_to_create=None):
        # create files corresponding to analysis results for a project object and its samples
        if files_to_create is None:
            files_to_create = self.generate_sample_analysis_results(proj_obj, sample_obj)
        for analysis_result_file in files_to_create:
            try:
                os.makedirs(os.path.dirname(analysis_result_file))
            except OSError as ose:
                pass
            open(analysis_result_file, "w").close()
        return files_to_create


class TestPreviousSampleAnalyses(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # set up test data
        cls.engine_name = "piper_ngi"
        cls.workflow_name = "merge_process_variantcall"

    def setUp(self):
        # create new test project objects
        self.project_data_objects = ProjectData()
        self.created_files = {}
        for proj_obj in self.project_data_objects.next():
            self.created_files[proj_obj.project_id] = {}
            for sample_obj in proj_obj.samples.values():
                self.created_files[proj_obj.project_id][sample_obj.name] = self.project_data_objects.create_sample_analysis_results(
                    proj_obj, sample_obj)
        self.maxDiff = None

    def tearDown(self):
        # remove test data files
        for proj_obj in self.project_data_objects.next():
            shutil.rmtree(proj_obj.base_path, ignore_errors=True)

    def test_find_previous_genotype_analyses(self):
        for proj_obj in self.project_data_objects.next():
            for sample_obj in proj_obj.samples.values():
                existing_files = filter(
                    lambda x: os.path.basename(os.path.dirname(x)) == "03_genotype_concordance",
                    self.created_files[proj_obj.project_id][sample_obj.name])
                self.assertEqual(
                    len(existing_files) > 0,
                    ngi_pipeline.engines.piper_ngi.utils.find_previous_genotype_analyses(proj_obj, sample_obj),
                    "find_previous_genotype_analyses did not evaluate to {} for {}:{}".format(
                        len(existing_files) > 0, proj_obj.project_id, sample_obj.name))

    def test_find_previous_sample_analyses(self):
        for proj_obj in self.project_data_objects.next():
            all_existing_files = []
            for sample_obj in proj_obj.samples.values():
                existing_files = filter(
                    lambda x: os.path.basename(os.path.dirname(x))[0:3] in ["0{}_".format(i+1) for i in xrange(8)],
                    self.created_files[proj_obj.project_id][sample_obj.name])
                self.assertItemsEqual(
                    existing_files,
                    ngi_pipeline.engines.piper_ngi.utils.find_previous_sample_analyses(proj_obj, sample_obj, include_genotype_files=True))#,
                   # "find_previous_sample_analyses did not return expected files for {}:{}".format(
                   #     proj_obj.project_id, sample_obj.name))
                # keep the existing_files for the final test
                all_existing_files.extend(existing_files)
            # filter out the genotyping files for the final test
            all_existing_files = filter(
                lambda x: os.path.basename(os.path.dirname(x))[0:3] != "03_",
                all_existing_files)
            self.assertItemsEqual(
                all_existing_files,
                ngi_pipeline.engines.piper_ngi.utils.find_previous_sample_analyses(proj_obj, sample_obj=None, include_genotype_files=False),
                "find_previous_sample_analyses did not return expected files for {}".format(proj_obj.project_id))

    def _remove_previous_type_analyses(self, result_type):

        def _wrap_test_function(expected_files, *args, **kwargs):
            observed_files = []

            def _in_existing_files(f):
                observed_files.append(f)
                if f not in expected_files:
                    raise AssertionError(f)

            with mock.patch('ngi_pipeline.engines.piper_ngi.utils.os.remove',
                            new=_in_existing_files), mock.patch(
                'ngi_pipeline.engines.piper_ngi.utils.shutil.rmtree', new=_in_existing_files):
                getattr(ngi_pipeline.engines.piper_ngi.utils, "remove_previous_{}_analyses".format(result_type))(*args, **kwargs)
                return observed_files

        for proj_obj in self.project_data_objects.next():
            for sample_obj in proj_obj.samples.values():
                if result_type == "sample":
                    existing_files = filter(
                        lambda x: os.path.basename(os.path.dirname(x))[0:3] in [
                            "0{}_".format(i + 1) for i in xrange(8) if i != 2],
                        self.created_files[proj_obj.project_id][sample_obj.name])
                    try:
                        removed_files = _wrap_test_function(existing_files, proj_obj, sample_obj=sample_obj)
                    except AssertionError as ae:
                        self.fail(
                            "remove_previous_{}_analyses unexpectedly tried to remove path {} for {}:{}".format(
                                result_type, ae.message, proj_obj.project_id, sample_obj.name))
                    self.assertItemsEqual(
                        existing_files,
                        removed_files,
                        "remove_previous_{}_analyses did not remove everything as expected for {}:{}".format(
                            result_type, proj_obj.project_id, sample_obj.name))
            existing_project_files = [
                existing_file for sample_files in self.created_files[proj_obj.project_id].values()
                for existing_file in sample_files]
            if result_type == "sample":
                existing_files = filter(
                    lambda x: os.path.basename(os.path.dirname(x))[0:3] in [
                        "0{}_".format(i + 1) for i in xrange(8) if i != 2],
                    existing_project_files)
                kwargs = {"sample_obj": None}
            else:
                existing_files = filter(
                    lambda x: os.path.basename(os.path.dirname(x))[0:3] == "03_",
                    existing_project_files)
                kwargs = {}
            try:
                removed_files = _wrap_test_function(existing_files, proj_obj, **kwargs)
            except AssertionError as ae:
                self.fail(
                    "remove_previous_{}_analyses unexpectedly tried to remove path {} for {}".format(
                        result_type, ae.message, proj_obj.project_id))
            self.assertItemsEqual(
                existing_files,
                removed_files,
                "remove_previous_{}_analyses did not remove everything as expected for {}".format(
                    result_type, proj_obj.project_id))

    def test_remove_previous_sample_analyses(self):
        self._remove_previous_type_analyses("sample")

    def test_remove_previous_genotype_analyses(self):
        self._remove_previous_type_analyses("genotype")
