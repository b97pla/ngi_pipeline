#!/bin/env python
"""Print the analysis status for samples, libpreps, and seqruns within a project."""

from __future__ import print_function

import argparse
import collections
import functools
import os
import sys
import tempfile
import time

from ngi_pipeline.conductor.flowcell import organize_projects_from_flowcell
from ngi_pipeline.database.classes import CharonSession, CharonError
from ngi_pipeline.engines.piper_ngi.local_process_tracking import update_charon_with_local_jobs_status
from ngi_pipeline.utils.classes import with_ngi_config
from ngi_pipeline.utils.filesystem import locate_project

print_stderr = functools.partial(print, file=sys.stderr)

def project_summarize(projects, verbosity=0):
    if type(verbosity) is not int or verbosity < 0:
        print_stderr('Invalid verbosity level ("{}"); must be a positive '
                     'integer; falling back to 0')
        verbosity = 0
    update_charon_with_local_jobs_status(quiet=True) # Don't send mails
    charon_session = CharonSession()
    projects_list = []
    for project in projects:
        try:
            project = os.path.basename(locate_project(project))
        except ValueError as e:
            print_stderr("Skipping project: {}".format(e))
            continue
        print_stderr('Gathering information for project "{}"...'.format(project))
        project_dict = {}
        try:
            project = charon_session.project_get(project)
        except CharonError as e:
            print_stderr('Project "{}" not found in Charon; skipping ({})'.format(project, e), file=sys.stderr)
            continue
        project_dict['name'] = project['name']
        project_dict['id'] = project['projectid']
        project_dict['status'] = project['status']
        samples_list = project_dict['samples'] = []
        for sample in charon_session.project_get_samples(project['projectid']).get('samples', []):
            sample_dict = {}
            sample_dict['id'] = sample['sampleid']
            sample_dict['analysis_status'] = sample['analysis_status']
            sample_dict['coverage'] = sample['total_autosomal_coverage']
            libpreps_list = sample_dict['libpreps'] = []
            samples_list.append(sample_dict)
            for libprep in charon_session.sample_get_libpreps(project['projectid'],
                                                              sample['sampleid']).get('libpreps', []):
                libprep_dict = {}
                libprep_dict['id'] = libprep['libprepid']
                libprep_dict['qc'] = libprep['qc']
                seqruns_list = libprep_dict['seqruns'] = []
                libpreps_list.append(libprep_dict)
                for seqrun in charon_session.libprep_get_seqruns(project['projectid'],
                                                                 sample['sampleid'],
                                                                 libprep['libprepid']).get('seqruns', []):
                    seqrun_dict = {}
                    seqrun_dict['id'] = seqrun['seqrunid']
                    seqrun_dict['alignment_status'] = seqrun['alignment_status']
                    seqrun_dict['coverage'] = seqrun['mean_autosomal_coverage']
                    if seqrun.get('total_reads'):
                        seqrun_dict['total_reads'] = seqrun['total_reads']
                    seqruns_list.append(seqrun_dict)
        projects_list.append(project_dict)


    if verbosity in (0, 1):
        projects_status_list = []
        #projects_by_status = collections.defaultdict(dict)
        #samples_by_status = collections.defaultdict(set)
        #libpreps_by_status = collections.defaultdict(set)
        #seqruns_by_status = collections.defaultdict(set)
        for project_dict in projects_list:
            project_status_dict = {}
            project_status_dict['name'] = "{} ({})".format(project_dict['name'], project_dict['id'])
            project_status_dict['status'] = project_dict['status']
            samples_by_status = project_status_dict['samples_by_status'] = collections.defaultdict(set)
            libpreps_by_status = project_status_dict['libpreps_by_status'] = collections.defaultdict(set)
            seqruns_by_status = project_status_dict['seqruns_by_status'] = collections.defaultdict(set)
            for sample_dict in project_dict.get('samples', []):
                #samples_by_status[sample_dict['analysis_status']].add(sample_dict['id'])
                sample_status = sample_dict['analysis_status']
                libpreps = sample_dict.get('libpreps')
                if libpreps:
                    if not any([libprep["seqruns"] for libprep in libpreps]):
                        sample_status = "NO_SEQRUNS"
                    else:
                        for libprep_dict in libpreps:
                            libpreps_by_status[libprep_dict['qc']].add(libprep_dict['id'])
                            for seqrun_dict in libprep_dict.get('seqruns', []):
                                seqruns_by_status[seqrun_dict['alignment_status']].add(seqrun_dict['id'])
                else:
                    sample_status = "NO_LIBPREPS"
                samples_by_status[sample_status].add(sample_dict['id'])
            projects_status_list.append(project_status_dict)

        print_items = (("Samples", "samples_by_status"),
                       ("Libpreps", "libpreps_by_status"),
                       ("Seqruns", "seqruns_by_status"),)

        for project_dict in projects_status_list:
            print_stderr("\nProject\n-------")
            print_stderr("    Name:   {:>40}".format(project_dict['name']))
            print_stderr("    Status: {:>40}".format(project_dict['status']))
            for name, dict_key in print_items:
                status_dict = project_dict[dict_key]
                print_stderr("{}\n{}".format(name, "-"*len(name)))
                total_items = sum(map(len, status_dict.values()))
                # Sort by analysis value
                for status, item_set in sorted(status_dict.iteritems(), key=lambda key_value: key_value[0]):
                    num_items = len(item_set)
                    percent = (100.00 * num_items) / total_items
                    print_stderr("    Status: {:<20} ({:>3}/{:<3}) ({:>6.2f}%)".format(status,
                                                                                       num_items,
                                                                                       total_items,
                                                                                       percent))
                    if verbosity == 1:
                        for item in sorted(item_set):
                            print_stderr("        {}".format(item))
            print_stderr("")

    else: # Verbosity is 2+, maximum verbosity
        output_template = "{}{:<30}{:>{rspace}}"
        for project_dict in projects_list:
            offset = 0
            indent = " " * offset
            rspace = 80 - offset
            print_stderr(output_template.format(indent, "Project name:", project_dict['name'], rspace=rspace))
            print_stderr(output_template.format(indent, "Project ID:", project_dict['id'], rspace=rspace))
            print_stderr(output_template.format(indent, "Project status:", project_dict['status'], rspace=rspace))
            for sample_dict in project_dict['samples']:
                print_stderr("")
                offset = 4
                indent = " " * offset
                rspace = 80 - offset
                print_stderr(output_template.format(indent, "Sample ID:", sample_dict['id'], rspace=rspace))
                print_stderr(output_template.format(indent, "Sample analysis status:", sample_dict['analysis_status'], rspace=rspace))
                print_stderr(output_template.format(indent, "Sample coverage:", sample_dict['coverage'], rspace=rspace))
                for libprep_dict in sample_dict['libpreps']:
                    print_stderr("")
                    offset = 8
                    indent = " " * offset
                    rspace = 80 - offset
                    print_stderr(output_template.format(indent, "Libprep ID:", libprep_dict['id'], rspace=rspace))
                    print_stderr(output_template.format(indent, "Libprep qc status:", libprep_dict['qc'], rspace=rspace))
                    for seqrun_dict in libprep_dict['seqruns']:
                        print_stderr("")
                        offset = 12
                        indent = " " * offset
                        rspace = 80 - offset
                        print_stderr(output_template.format(indent, "Seqrun ID:", seqrun_dict['id'], rspace=rspace))
                        print_stderr(output_template.format(indent, "Seqrun alignment status:", seqrun_dict['alignment_status'], rspace=rspace))
                        print_stderr(output_template.format(indent, "Seqrun mean auto. coverage:", seqrun_dict['coverage'], rspace=rspace))
                        if "total_reads" in seqrun_dict:
                            print_stderr(output_template.format(indent, "Seqrun total reads:", seqrun_dict['total_reads'], rspace=rspace))
            print_stderr("\n")


def flowcell_summarize(flowcells, brief=False, verbose=False):
    projects_to_analyze = \
            organize_projects_from_flowcell(demux_fcid_dirs=flowcells,
                                            quiet=True, create_files=False)
    project_summarize(projects_to_analyze, brief=brief, verbose=verbose)


if __name__=="__main__":
    parser = argparse.ArgumentParser()
    verbosity = parser.add_mutually_exclusive_group()
    verbosity.add_argument("-v", "--verbosity", action="count", default=0,
            help="Increase output verbosity (try -v, -vv)")
    subparsers = parser.add_subparsers(help="Summarize project or flowcell.")
    project_parser = subparsers.add_parser('project')
    project_parser.add_argument('project_dirs', nargs='+',
            help=('The name or ID (in Charon) of or path to one or more projects to be summarized.'))

    flowcell_parser = subparsers.add_parser('flowcell')
    flowcell_parser.add_argument('flowcell_dirs', nargs='+',
            help=('The name of or path to one or more flowcell directories to be summarized.'))


    args = parser.parse_args()

    if "project_dirs" in args:
        project_summarize(args.project_dirs, args.verbosity)
    elif "flowcell_dirs" in args:
        flowcell_summarize(args.flowcell_dirs, args.verbosity)
    else:
        parser.print_usage()
