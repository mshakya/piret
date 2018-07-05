#! /usr/bin/env python
"""A suite to test confirm_fasta."""
import os
import sys
import unittest
dir_path = os.path.dirname(os.path.realpath(__file__))
lib_path = os.path.abspath(os.path.join(dir_path, '..'))
sys.path.append(lib_path)
from pypiret.Runs.srna import extract_pp 

class TestPP(unittest.TestCase):
    """Unittest testcase."""
    bam_file = "tests/data/samp5.bam"
    
    def test_split(self):
        """Test if true."""
        extract_pp.prop_paired(self, self.bam_file)
        assert os.path.exists(self.bam_file) == 1

    def test_split(self):
        """Test if true."""
        extract_pp.merge_prop_paired(self, self.bam_file)
        assert os.path.exists(self.bam_file) == 1

if __name__ == '__main__':
    unittest.main()
