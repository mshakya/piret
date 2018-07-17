#! /usr/bin/env python

"""Check design."""
from __future__ import print_function
import sys
import pandas as pd
import os


class CheckDesign():
    """Check experimental design are tab delimited."""

    def __init__(self, exp):
        """Initialize."""
        self.design_file = exp

    def tab(self):
        """Check if the file is tab delimited and has more than 2 columns."""
        with open(self.design_file, 'r') as dfile:
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
                    sys.exit('Experimental design file is not correctly \
                              formatted!')
            else:
                sys.exit('Experimental design file is not correctly formatted!\n \
                    There is only one or no column!'
                         )

    def header(self):
        """Check if headers are present."""
        with open(self.design_file, 'r') as dfile:
            header = dfile.readline().rstrip()
            if header == "#SampleID\tFiles\tGroup":
                return True
            else:
                sys.exit('The header of experimental design file must be #SampleID\tFiles\tGroup\n')

    def sample_name(self):
        """Check if sample names are unique."""
        with open(self.design_file, 'r') as dfile:
            next(dfile)
            sample_names = []
            for line in dfile:
                sample_names.append(line.split("\t")[0])
            if len(sample_names) == len(set(sample_names)):
                return True
            else:
                return False
                sys.exit('Sample names must be unique!')

    def file_name(self):
        """Check if file names are unique."""
        with open(self.design_file, 'r') as dfile:
            next(dfile)
            file_names = []
            for line in dfile:
                file_names.append(line.split("\t")[1])
            if len(file_names) == len(set(file_names)):
                return True
            else:
                sys.exit('Input files must be unique!')

    def group_name(self):
        """Check if there are at least two groups."""
        with open(self.design_file, 'r') as dfile:
            next(dfile)
            groups = []
            for line in dfile:
                groups.append(line.split("\t")[2])
            if len(groups) < 2:
                sys.exit('There must be at least two groups for comparisons!')
            else:
                return True

    def fastq_exists(self):
        """Check if specified fastq files are present."""
        fastq_files = self.extract_sample_fastqs()
        for samp, fastqs in fastq_files.items():
            if isinstance(fastqs, str):
                files = fastqs.split(":")
                for f in files:
                    if os.path.exists(f):
                        pass
                    else:
                        err_message = "%s does not exist"
                        sys.exit(err_message % f)
            elif isinstance(fastqs, list):
                for pairs in fastqs:
                    files = pairs.split(";")
                    for f in files:
                        if os.path.exists(f):
                            pass
                        else:
                            err_message = "%s does not exist"
                            sys.exit(err_message % f)
        return True

    def extract_sample_fastqs(self):
        """Extract sample and fastqc from design file."""
        with open(self.design_file, 'r') as dfile:
            next(dfile)
            sample = []
            sample_fastqs = {}
            for line in dfile:
                sample = line.split("\t")[0]
                if ';' in line.split("\t")[1]:
                    fastq_pairs = line.split("\t")[1].split(";")
                else:
                    fastq_pairs = line.split("\t")[1]
                sample_fastqs[sample] = fastq_pairs
        return sample_fastqs

    def sample_suff(self, method):
        """Check if enough number of samples are present for Deseq analysis."""
        exp_design = pd.read_csv(self.design_file, sep='\t')
        if method in ['both', 'Deseq']:
            group_count = exp_design['Group'].value_counts()
            for e in [1, 2]:
                if e in group_count.values:
                    sys.exit("There must be at least three samples \
                              to run DeSeq")
                    break
                else:
                    pass
        return True
