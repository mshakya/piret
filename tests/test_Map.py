#! /usr/bin/env python
"""A suite to test file_exist."""
import sys
import os
import unittest
dir_path = os.path.dirname(os.path.realpath(__file__))
lib_path = os.path.abspath(os.path.join(dir_path, '..'))
sys.path.append(lib_path)
from pypiret.Runs import Map  # noqa


class TestMap(unittest.TestCase):
    """Unittest testcase."""

    def test_hisat_map(self):
        """Test if program SystemExit when tab delimited file is not given."""
        Map.HisatIndex(fasta="tests/data/test_prok.fa",
                       hi_index="tests/test_index",
                       bindir="bin",
                       numCPUs=1)

    def tearDown(self):
        """Tear down."""
        self.test_hisat_map()

# class Test
if __name__ == '__main__':
    unittest.main()
