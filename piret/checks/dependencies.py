#! /usr/bin/env python

"""Check design."""
from __future__ import print_function
import sys
from plumbum import local


class CheckDependencies():
    """Check if third party dependencies are in the path."""

    def __init__(self, package):
        """Initialize."""
        self.package = package
        # self.logger = logger

    def check_thirdparty(self):
        """Check if thirdparty tool is in the path."""
        try:
            local[self.package]
            print('%s package exist' % self.package)
            # self.logger.info('%s package exist' % self.package)
            return True
        except:
            exit_msg = """%s is not installed or not in your PATH!
                          Please install. 
                          See README for instructions on how to install it."""
            sys.exit(exit_msg % self.package)


def check_depen(tool_list):
    """A function that checks the list of tool are in path or not."""
    for tool in tool_list:
        check = CheckDependencies(package=tool)
        check.check_thirdparty()

def check_pydepen():
    """A function that checks the list of tool are in path or not."""
    try:
         import luigi
         import plumbum
         import pandas
         import Bio
         import gffutils
         import os
         import sys
         import argparse
         import logging
         print("All required python package exist")
    except:
        exit_msg = """python package luigi is not installed! \n Please install
        it using python setup.py install and try again"""
        sys.exit(exit_msg)