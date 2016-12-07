#! /usr/bin/env python


"""
Create log instance and files.

log files
"""
import logging


def create_logger(workdir, type, level):
    """
    Create logger instances.

    also creates WORKDIR/process.log
    """
    logger = logging.getLogger(type)
    formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
    fn = '/'.join([workdir, '.'.join([type, 'log'])])
    hdlr = logging.FileHandler(fn)
    hdlr.setFormatter(formatter)
    logger.addHandler(hdlr)
    logger.setLevel(level)
    return logger
