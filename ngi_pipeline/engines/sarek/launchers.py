from ngi_pipeline.engines.sarek.database import CharonConnector, TrackingConnector
from ngi_pipeline.engines.sarek.exceptions import SampleNotValidForAnalysisError
from ngi_pipeline.engines.sarek.models import SarekAnalysis
from ngi_pipeline.engines.sarek.process import SlurmConnector


def analyze(analysis_object):
    """
    This is the main entry point for launching the Sarek analysis pipeline. This gets called by NGI pipeline for
    projects having the corresponding best_practice_analysis in Charon. It's called per project and the passed analysis
    object contains (some?) parameters for the analysis.

    Refer to the ngi_pipeline.conductor.launchers.launch_analysis method to see how the analysis object is created.

    :param analysis_object: an ngi_pipeline.conductor.classes.NGIAnalysis object holding parameters for the analysis
    :return: None
    """
    # get a SlurmConnector that will take care of submitting analysis jobs
    slurm_project_id = analysis_object.config["environment"]["project_id"]
    slurm_mail_user = analysis_object.config["mail"]["recipient"]
    slurm_conector = SlurmConnector(
        slurm_project_id, slurm_mail_user, cwd="/scratch", slurm_mail_events="TIME_LIMIT_80")

    # get a CharonConnector that will interface with the Charon database
    charon_connector = CharonConnector(analysis_object.config, analysis_object.log)

    # get a TrackingConnector that will interface with the local SQLite database
    tracking_connector = TrackingConnector(analysis_object.config, analysis_object.log)

    analysis_object.log.info("Launching SAREK analysis for {}".format(analysis_object.project.project_id))
    analysis_engine = SarekAnalysis.get_analysis_instance_for_project(
        analysis_object.project.project_id,
        analysis_object.config,
        analysis_object.log,
        charon_connector=charon_connector,
        tracking_connector=tracking_connector,
        process_connector=slurm_conector)

    # iterate over the samples and launch analysis for each
    for sample_object in analysis_object.project:
        try:
            analysis_engine.analyze_sample(sample_object, analysis_object)
        except SampleNotValidForAnalysisError as e:
            analysis_object.log.error(e)
