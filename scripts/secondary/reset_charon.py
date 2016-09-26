"""Keeps track of running workflow processes"""
import json
import shelve
import os
import glob
import re


import argparse
import time
import os

from ngi_pipeline.database.classes import CharonSession, CharonError
from ngi_pipeline.log.loggers import minimal_logger


def main(project):

    charon_session = CharonSession()
    samples = charon_session.project_get_samples(project)
    for sample in samples["samples"]:
        
        charon_session.sample_update(project, sample["sampleid"],
                                     status="FRESH",
                                     analysis_status= "TO_ANALYZE",
                                     genotype_status=None,
                                     total_autosomal_coverage="0",
                                     total_sequenced_reads="0")
        for sample_prep in charon_session.sample_get_libpreps(project, sample["sampleid"])['libpreps']:
            seqruns = charon_session.libprep_get_seqruns(project, sample["sampleid"], sample_prep["libprepid"])['seqruns']
            for seqrun in seqruns:
                charon_session.seqrun_update(project, sample["sampleid"], sample_prep["libprepid"], seqrun["seqrunid"],
                                             mean_autosomal_coverage = "0",
                                             alignment_status  = "NOT_RUNNING")




if __name__ == '__main__':
    parser = argparse.ArgumentParser("Clean the specify project: restore it to default --- no sample will be deleted")
    parser.add_argument("-p", "--project", dest="project", help=("Project to restore "))
    args_dict = vars(parser.parse_args())
    main(**args_dict)


