import json
import os

from ngi_pipeline.conductor.classes import NGIProject
from ngi_pipeline.database.classes import CharonSession, CharonError
from ngi_pipeline.log.loggers import minimal_logger
from ngi_pipeline.utils.filesystem import locate_chip_genotypes_for_project

LOG = minimal_logger(__name__)


def create_charon_entries_from_project(project, best_practice_analysis="whole_genome_reseq",
                                       sequencing_facility="NGI-S",
                                       force_overwrite=False, delete_existing=False,
                                       retry_on_fail=True):
    """Given a project object, creates the relevant entries in Charon.
    This code is remarkably shoddy as I created it in a hurry and then later
    it became a part of the pipeline. Use at your own risk! Ha ha.

    :param NGIProject project: The NGIProject object
    :param str best_practice_analysis: The workflow to assign for this project (default "variant_calling")
    :param str sequencing_facility: The facility that did the sequencing
    :param bool force_overwrite: If this is set to true, overwrite existing entries in Charon (default false)
    :param bool delete_existing: Don't just update existing entries, delete them and create new ones (default false)
    """
    charon_session = CharonSession()
    update_failed=False
    try:
        status = "OPEN"
        LOG.info('Creating project "{}" with status "{}", best practice analysis "{}", '
                 'and sequencing_facility {}'.format(project, status,
                                                     best_practice_analysis,
                                                     sequencing_facility))
        charon_session.project_create(projectid=project.project_id,
                                      name=project.name,
                                      status=status,
                                      best_practice_analysis=best_practice_analysis,
                                      sequencing_facility=sequencing_facility)
        LOG.info('Project "{}" created in Charon.'.format(project))
    except CharonError as e:
        if e.status_code == 400:
            if force_overwrite:
                LOG.warn('Overwriting data for project "{}"'.format(project))
                charon_session.project_update(projectid=project.project_id,
                                              name=project.name,
                                              status=status,
                                              best_practice_analysis=best_practice_analysis,
                                              sequencing_facility=sequencing_facility)
                LOG.info('Project "{}" updated in Charon.'.format(project))
            else:
                LOG.info('Project "{}" already exists; moving to samples...'.format(project))
        else:
            raise
    for sample in project:
        # TODO: parse genotype files and set status only on genotyped samples instead of on all samples in project
        genotype_status = "AVAILABLE" if project.chip_genotypes is not None else None
        if delete_existing:
            LOG.warn('Deleting existing sample "{}"'.format(sample))
            try:
                charon_session.sample_delete(projectid=project.project_id,
                                             sampleid=sample.name)
            except CharonError as e:
                update_failed=True
                LOG.error('Could not delete sample "{}": {}'.format(sample, e))
        try:
            analysis_status = "TO_ANALYZE"
            LOG.info('Creating sample "{}" with analysis_status "{}"'.format(sample, analysis_status))
            charon_session.sample_create(projectid=project.project_id,
                                         sampleid=sample.name,
                                         analysis_status=analysis_status,
                                         genotype_status=genotype_status)
            LOG.info('Project/sample "{}/{}" created in Charon.'.format(project, sample))
        except CharonError as e:
            if e.status_code == 400:
                if force_overwrite:
                    LOG.warn('Overwriting data for project "{}" / '
                             'sample "{}"'.format(project, sample))
                    charon_session.sample_update(projectid=project.project_id,
                                                 sampleid=sample.name,
                                                 analysis_status=analysis_status,
                                                 genotype_status=genotype_status)
                    LOG.info('Project/sample "{}/{}" updated in Charon.'.format(project, sample))
                else:
                    LOG.info('Project "{}" / sample "{}" already exists; moving '
                             'to libpreps'.format(project, sample))
            else:
                update_failed=True
                LOG.error(e)
                continue
        for libprep in sample:
            if delete_existing:
                LOG.warn('Deleting existing libprep "{}"'.format(libprep))
                try:
                    charon_session.libprep_delete(projectid=project.project_id,
                                                  sampleid=sample.name,
                                                  libprepid=libprep.name)
                except CharonError as e:
                    LOG.warn('Could not delete libprep "{}": {}'.format(libprep, e))
            try:
                qc = "PASSED"
                LOG.info('Creating libprep "{}" with qc status "{}"'.format(libprep, qc))
                charon_session.libprep_create(projectid=project.project_id,
                                              sampleid=sample.name,
                                              libprepid=libprep.name,
                                              qc=qc)
                LOG.info(('Project/sample/libprep "{}/{}/{}" created in '
                          'Charon').format(project, sample, libprep))
            except CharonError as e:
                if e.status_code == 400:
                    if force_overwrite:
                        LOG.warn('Overwriting data for project "{}" / '
                                 'sample "{}" / libprep "{}"'.format(project, sample,
                                                                     libprep))
                        charon_session.libprep_update(projectid=project.project_id,
                                                      sampleid=sample.name,
                                                      libprepid=libprep.name,
                                                      qc=qc)
                        LOG.info(('Project/sample/libprep "{}/{}/{}" updated in '
                                  'Charon').format(project, sample, libprep))
                    else:
                        LOG.debug(e)
                        LOG.info('Project "{}" / sample "{}" / libprep "{}" already '
                                 'exists; moving to libpreps'.format(project, sample, libprep))
                else:
                    update_failed=True
                    LOG.error(e)
                    continue
            for seqrun in libprep:
                if delete_existing:
                    LOG.warn('Deleting existing seqrun "{}"'.format(seqrun))
                    try:
                        charon_session.seqrun_delete(projectid=project.project_id,
                                                     sampleid=sample.name,
                                                     libprepid=libprep.name,
                                                     seqrunid=seqrun.name)
                    except CharonError as e:
                        update_failed=True
                        LOG.error('Could not delete seqrun "{}": {}'.format(seqrun, e))
                try:
                    alignment_status="NOT_RUNNING"
                    LOG.info('Creating seqrun "{}" with alignment_status "{}"'.format(seqrun, alignment_status))
                    charon_session.seqrun_create(projectid=project.project_id,
                                                 sampleid=sample.name,
                                                 libprepid=libprep.name,
                                                 seqrunid=seqrun.name,
                                                 alignment_status=alignment_status,
                                                 total_reads=0,
                                                 mean_autosomal_coverage=0)
                    LOG.info(('Project/sample/libprep/seqrun "{}/{}/{}/{}" '
                              'created in Charon').format(project, sample,
                                                          libprep, seqrun))
                except CharonError as e:
                    if e.status_code == 400:
                        if force_overwrite:
                            LOG.warn('Overwriting data for project "{}" / '
                                     'sample "{}" / libprep "{}" / '
                                     'seqrun "{}"'.format(project, sample,
                                                          libprep, seqrun))
                            charon_session.seqrun_update(projectid=project.project_id,
                                                         sampleid=sample.name,
                                                         libprepid=libprep.name,
                                                         seqrunid=seqrun.name,
                                                         alignment_status=alignment_status,
                                                         total_reads=0,
                                                         mean_autosomal_coverage=0)
                            LOG.info(('Project/sample/libprep/seqrun "{}/{}/{}/{}" '
                                      'updated in Charon').format(project, sample,
                                                                  libprep, seqrun))
                        else:
                            LOG.info('Project "{}" / sample "{}" / libprep "{}" / '
                                     'seqrun "{}" already exists; next...'.format(project, sample,
                                                                                  libprep, seqrun))
                    else:
                        update_failed=True
                        LOG.error(e)
                        continue

    if update_failed:
        if retry_on_fail:
            create_charon_entries_from_project(project, best_practice_analysis=best_practice_analysis,
                                       sequencing_facility=sequencing_facility,
                                       force_overwrite=force_overwrite, delete_existing=delete_existing,
                                       retry_on_fail=False)
        else:
            raise CharonError("A network error blocks Charon updating.")

def recreate_project_from_db(analysis_top_dir, project_name, project_id):
    project_dir = os.path.join(analysis_top_dir, "DATA", project_name)
    # locate chip genotype files and recreate the chip_genotypes based on those
    chip_genotype_files = locate_chip_genotypes_for_project(
        project_name,
        config={'qc': {'chip_genotypes_root': [os.path.join(analysis_top_dir, "DATA")]}})
    project_obj = NGIProject(name=project_name,
                             dirname=project_name,
                             project_id=project_id,
                             base_path=analysis_top_dir,
                             chip_genotypes=map(
                                 lambda f: NGIChipGenotypes(name=os.path.basename(f)),
                                 chip_genotype_files) if chip_genotype_files is not None else None)
    charon_session = CharonSession()
    try:
        samples_dict = charon_session.project_get_samples(project_id)["samples"]
    except CharonError as e:
        raise RuntimeError("Could not access samples for project {}: {}".format(project_id, e))
    for sample in samples_dict:
        sample_id = sample.get("sampleid")
        sample_dir = os.path.join(project_dir, sample_id)
        sample_obj = project_obj.add_sample(name=sample_id, dirname=sample_id)
        sample_obj.status = sample.get("status", "unknown")
        try:
            libpreps_dict = charon_session.sample_get_libpreps(project_id, sample_id)["libpreps"]
        except CharonError as e:
            raise RuntimeError("Could not access libpreps for project {} / sample {}: {}".format(project_id,sample_id, e))
        for libprep in libpreps_dict:
            libprep_id = libprep.get("libprepid")
            libprep_obj = sample_obj.add_libprep(name=libprep_id,  dirname=libprep_id)
            libprep_obj.status = libprep.get("status", "unknown")
            try:
                seqruns_dict = charon_session.libprep_get_seqruns(project_id, sample_id, libprep_id)["seqruns"]
            except CharonError as e:
                raise RuntimeError("Could not access seqruns for project {} / sample {} / "
                                   "libprep {}: {}".format(project_id, sample_id, libprep_id, e))
            for seqrun in seqruns_dict:
                # e.g. 140528_D00415_0049_BC423WACXX
                seqrun_id = seqrun.get("seqrunid")
                seqrun_obj = libprep_obj.add_seqrun(name=seqrun_id, dirname=seqrun_id)
                seqrun_obj.status = seqrun.get("status", "unknown")
    return project_obj
