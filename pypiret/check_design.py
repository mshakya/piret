#! /usr/bin/env python

import sys


class CheckDesign():
    """Check experimental design are tab delimited."""

    def __init__(self):
        """Initialize."""
        # self.design_file = design_file

    def check_tab(self, design_file):
        """Check if the file is tab delimited and has more than 2 columns."""
        with open(design_file, 'r') as dfile:
            header = dfile.readline()
            num_tab = len(header.split("\t"))
            if num_tab > 2:
                line_status = []
                for line in dfile:
                    if len(line.split("\t")) == num_tab:
                        line_status.append(True)
                    else:
                        line_status.append(False)
                if all(line_status) is True:
                    return True
                else:
                    return False
                    sys.exit('Experimental design file is not correctly formatted!'
                             )
            else:
                return False
                sys.exit('Experimental design file is not correctly formatted!'
                         )

    def check_header(self, design_file):
        """Check if headers are present."""
        with open(design_file, 'r') as dfile:
            header = dfile.readline()

            if header.split("\t")[0] == "ID":
                pass
            else:
                return False
                sys.exit('The first column header must be labeled as ID!')

            if header.split("\t")[1] == "Files":
                pass
            else:
                return False
                sys.exit('The second column header must be labeled as Files!')

            if header.split("\t")[2] == "group":
                pass
            else:
                return False
                sys.exit('The third column header must be labeled as Files!')

    def check_sample_name(self, design_file):
        """Check if sample names are unique."""
        with open(design_file, 'r') as dfile:
            next(dfile)
            sample_names = []
            for line in dfile:
                sample_names.append(line.split("\t")[0])
            if len(sample_names) == len(set(sample_names)):
                return True
            else:
                return False
                sys.exit('Sample names must be unique!')

    def check_file_name(self, design_file):
        """Check if file names are unique."""
        with open(design_file, 'r') as dfile:
            next(dfile)
            file_names = []
            for line in dfile:
                file_names.append(line.split("\t")[1])
            if len(file_names) == len(set(file_names)):
                return True
            else:
                return False
                sys.exit('Input files must be unique!')

    def check_group_name(self, design_file):
        """Check if there are at least two groups."""
        with open(design_file, 'r') as dfile:
            next(dfile)
            groups = []
            for line in dfile:
                groups.append(line.split("\t")[2])
            if len(groups) < 2:
                return False
                sys.exit('There must be at least two groups for comparisons!')
            else:
                return True
