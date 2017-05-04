#! /usr/bin/env python
"""A suite to test file_exist."""
import sys
import unittest
sys.path.append('../')
from pypiret.Runs import Map  # noqa


class TestMap(unittest.TestCase):
    """Unittest testcase."""

    def test_hisat_map(self):
        """Test if program SystemExit when tab delimited file is not given."""
        design = Map.HisatIndex(fasta="tests/data/test_prok.fa",
                                hi_index="tests/test_index",
                                bindir="bin")

    def tearDown(self):
        """Tear down."""
        self.widget.dispose()

# class Test
if __name__ == '__main__':
    unittest.main()
