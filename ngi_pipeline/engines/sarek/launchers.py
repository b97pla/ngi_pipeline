from ngi_pipeline.engines.sarek.models import SarekAnalysis


def analyze(analysis_object):
    """
    This is the main entry point for launching the Sarek analysis pipeline. This gets called by NGI pipeline for
    projects having the corresponding best_practice_analysis in Charon. It's called per project and the passed analysis
    object contains (some?) parameters for the analysis.

    Refer to the ngi_pipeline.conductor.launchers.launch_analysis method to see how the analysis object is created.

    :param analysis_object: an ngi_pipeline.conductor.classes.NGIAnalysis object holding parameters for the analysis
    :return: None
    """
    analysis_object.log.info("Launching SAREK analysis for {}".format(analysis_object.project.project_id))
    analysis_engine = SarekAnalysis.get_analysis_instance_for_project(
        analysis_object.project.project_id,
        analysis_object.config,
        analysis_object.log,
        charon_connector=None)

    # filter out samples that do not fulfil the conditions to be launched
    samples_to_launch = filter(
        lambda sample: analysis_engine.sample_should_be_started(analysis_object.project.project_id, sample.name),
        analysis_object.project)

    # launch analysis for each sample
    map(lambda sample: analysis_engine.analyze_sample(sample, analysis_object), samples_to_launch)


