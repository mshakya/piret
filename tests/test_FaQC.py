#! /usr/bin/env python

import sys
import os
import unittest
import shutil
from luigi.interface import build
DIR = os.path.dirname(os.path.realpath(__file__))
LIB = os.path.abspath(os.path.join(DIR, '..'))
sys.path.append(LIB)
from pypiret.Runs import FaQC




class TestFaQC(unittest.TestCase):
    """Unittest testcase."""

    def file_len(self, fname):
        """A function to test filename"""
        with open(fname) as f:
            for i, l in enumerate(f):
                pass
        return i + 1

    def test_faqc_function(self):
        """Test FaQC function."""
        dirpath = os.path.dirname(os.path.realpath(__file__))
        bindir = os.path.abspath(os.path.join(dirpath, '..', 'bin'))
        build([FaQC.PairedRunQC(fastqs=["tests/data/fastqs/BTT_test15_R1.fastq.gz",
                                        "tests/data/fastqs/BTT_test15_R1.fastq.gz"],
                                sample="samp1",
                                num_cpus=1,
                                qc_outdir="test_faqc",
                                faqc_min_L=50,
                                n_cutoff=5)],
              local_scheduler=True, workers=1)

    def test_fq_out(self):
        """Test if fq has expected line number"""
        fq_lines = self.file_len("test_faqc/samp1.1.trimmed.fastq")
        assert(fq_lines) == 189508

    def test_output(self):
        """Test if 5 files are generated as output"""
        filenames = [x[2] for x in os.walk("test_faqc")][0]
        assert len(filenames) == 5

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree('test_faqc')


if __name__ == '__main__':
    unittest.main()
