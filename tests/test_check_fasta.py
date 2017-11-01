#! /usr/bin/env python
"""A suite to test confirm_fasta."""
import os
import sys
import unittest
dir_path = os.path.dirname(os.path.realpath(__file__))
lib_path = os.path.abspath(os.path.join(dir_path, '..'))
sys.path.append(lib_path)
from pypiret import CheckFasta  # noqa


class TestConfirmFasta(unittest.TestCase):
    """Unittest testcase."""

    def test_true(self):
        """Test if true."""
        fasta = CheckFasta()
        self.assertTrue(fasta.confirm_fasta("tests/data/test_fa.fa"))

    def test_raises(self):
        """Test if the TypeError is raised."""
        fasta = CheckFasta()
        with self.assertRaises(TypeError):
            fasta.confirm_fasta("tests/data/test_experimental_design.csv")


if __name__ == '__main__':
    unittest.main()
