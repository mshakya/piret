#! /usr/bin/env python
"""A suite to test file_exist."""
import sys
import os
import unittest
import pytest
dir_path = os.path.dirname(os.path.realpath(__file__))
lib_path = os.path.abspath(os.path.join(dir_path, '..'))
sys.path.append(lib_path)
from pypiret.Runs import FastQC  # noqa
import luigi
import shutil

class TestFaQC(unittest.TestCase):
    """Unittest testcase."""

    def file_len(self,fname):
        with open(fname) as f:
            for i, l in enumerate(f):
                pass
        return i + 1

    def test_faqc_function(self):
        """Test FaQC function."""
        dir_path = os.path.dirname(os.path.realpath(__file__))
        bindir = os.path.abspath(os.path.join(dir_path, '..', 'bin'))
        luigi.interface.build([FastQC.PairedRunQC(fastqs=["tests/data/fastqs/BTT_test15_R1.fastq",
                                                       "tests/data/fastqs/BTT_test15_R1.fastq"],
                                                sample="samp1",
                                                numCPUs=1,
                                                qc_outdir="test_faqc",
                                                bindir=bindir,
                                                faqc_min_L=50,
                                                n_cutoff=5)],
                                                local_scheduler=True,
                                                workers=1)
        

    def test_fq_out(self):
        """Test if fq has expected line number"""
        fq_lines = self.file_len("test_faqc/samp1.1.trimmed.fastq")
        assert(fq_lines) == 188420

    def test_output(self):
        """Test if 5 files are generated as output"""
        filenames =  [x[2] for x in os.walk("test_faqc")][0]
        assert len(filenames) == 5
    
    @classmethod
    def tearDownClass(cls):
        shutil.rmtree('test_faqc')



if __name__ == '__main__':
    unittest.main()
