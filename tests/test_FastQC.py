#! /usr/bin/env python
"""A suite to test file_exist."""
import sys
import os
import unittest
dir_path = os.path.dirname(os.path.realpath(__file__))
lib_path = os.path.abspath(os.path.join(dir_path, '..'))
sys.path.append(lib_path)
from pypiret.Runs import FastQC  # noqa


class TestFaQC(unittest.TestCase):
    """Unittest testcase."""

    def test_faqc_function(self):
        """Test."""
        FastQC.PairedRunQC(fastqs=["tests/data/fastqs/BTT_test15_R1.fastq",
                                   "tests/data/fastqs/BTT_test15_R1.fastq"],
                           sample="samp1",
                           numCPUs=1,
                           outdir="test_faqc",
                           bindir="bin")

    def tearDown(self):
        """Tear down."""
        self.test_faqc_function()


if __name__ == '__main__':
    unittest.main()
