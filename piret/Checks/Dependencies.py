#! /usr/bin/env python

"""Check design."""
from __future__ import print_function
import sys
from plumbum import local


class CheckDependencies():
    """Check if third party dependencies are in the path."""

    def __init__(self, package, logger):
        """Initialize."""
        self.package = package
        self.logger = logger

    def check_thirdparty(self):
        """Check if thirdparty tool is in the path."""
        try:
            local[self.package]
            print('%s package exist' % self.package)
            # self.logger.info('%s package exist' % self.package)
            return True
        except:
            exit_msg = """%s is not installed or not in your PATH! Please install. See README for instructions on how to install it."""
            sys.exit(exit_msg % self.package)

    def check_pypackage(self):
        """Check if python package is installed."""
        try:
            import self.package
            return True
        except:
            exit_msg = """python package %s is not installed! \n Please install
            it and try again. Run INSTALL.sh or read README for instruction on
            installations"""
            sys.exit(exit_msg % self.package)

    def check_all_thirdparty():
        """Check if given list of third party tools are installed."""


thirdparty_tools = ["bamtools", "samtools", ]

# CheckDependencies("bamtools")

def check_depen(tool_list, logger):
    """A function that checks the list of tool are in path or not."""
    for tool in tool_list:
        check = CheckDependencies(package=tool, logger=logger)
        check.check_thirdparty()