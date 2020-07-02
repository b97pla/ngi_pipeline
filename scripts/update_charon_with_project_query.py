#!/bin/env python

import argparse
import os

from ngi_pipeline.engines import sarek
from ngi_pipeline.engines.sarek.database import TrackingConnector
from ngi_pipeline.utils.classes import with_ngi_config
from ngi_pipeline.log.loggers import minimal_logger

LOG = minimal_logger(__name__, debug=True)


class DiskTrackingSession(object):
    """
    This is an object to replace the SQLAlchemy object injected into the TrackingConnector, in order to replace the
    database connections
    """
    def __init__(self, analyses=None):
        self.analyses = analyses or list()

    def add(self, db_obj):
        self.analyses.append(db_obj)

    def all(self, *args, **kwargs):
        for analysis in self.analyses:
            yield analysis

    def commit(self, *args, **kwargs):
        pass

    def delete(self, *args, **kwargs):
        pass

    def filter(self, *args, **kwargs):
        return self

    def query(self, *args, **kwargs):
        return self


@with_ngi_config
def update_charon_with_project(project, sample=None, config=None):
    project_base_path = os.path.join(
        config["analysis"]["base_root"],
        config["analysis"]["upps_root"],
        config["analysis"]["top_dir"])

    project_analysis_dir = os.path.join(
        project_base_path,
        "ANALYSIS",
        project)

    db_session = DiskTrackingSession()

    for sample in os.listdir(project_analysis_dir):
        db_session.add(
            TrackingConnector._SampleAnalysis(
                project_id=project,
                project_name=project,
                sample_id=os.path.basename(sample),
                project_base_path=project_base_path,
                workflow="SarekGermlineAnalysis",
                engine="sarek",
                process_id=999999999999)
            )

    tracking_connector = TrackingConnector(
        config,
        LOG,
        tracking_session=db_session)
    sarek.local_process_tracking.update_charon_with_local_jobs_status(
        config=config,
        log=LOG,
        tracking_connector=tracking_connector)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Update Charon with Sarek analysis status and analysis results for a project independent of the "
                    "local processing db")

    parser.add_argument("-p", "--project", required=True)
    parser.add_argument("-s", "--sample", required=False)
    parser.add_argument("-c", "--config", required=False)

    args = parser.parse_args()
    update_charon_with_project(args.project, sample=args.get("sample"), config=args.get("config"))
