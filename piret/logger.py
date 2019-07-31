#! /usr/bin/env python

import os
"""
Create log instance and files.

log files
"""
import logging


def create_logger(workdir):
    """
    Create logger instances.

    also creates WORKDIR/process.log

    """
    log_file = os.path.join(workdir, "piret.log")
    if os.path.exists(workdir) is False:
        os.makedirs(workdir)
    open(log_file, 'a').close()
    logger = logging.basicConfig(format='%(asctime)s %(message)s',
                                 filename=log_file, level=logging.INFO)
    return logger
