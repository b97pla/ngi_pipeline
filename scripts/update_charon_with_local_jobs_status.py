#!/bin/env python

import argparse
import importlib


if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-e", "--engine", required=True)

    # E.g. piper
    engine = parser.parse_args().engine.lower()
    #if not engine.endswith("_ngi"):
        # half-hearted attempt to make this more flexible
    #    engine = "{}_ngi".format(engine)

    module = "ngi_pipeline.engines.{}".format(engine)
    ## This should at some point be refactored. Sigh.
    update_function = importlib.import_module(module).local_process_tracking.update_charon_with_local_jobs_status
    update_function()
