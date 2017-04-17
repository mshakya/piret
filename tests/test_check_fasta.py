#! /usr/bin/env python
"""A suite to test confirm_fasta."""
import sys
import unittest
sys.path.append('../')
from pypiret import check_fasta  # noqa


class TestConfirmFasta(unittest.TestCase):
    """Unittest testcase."""

    def test_true(self):
        """Test if true."""
        self.assertTrue(check_fasta.confirm_fasta("tests/data/test_fa.fa"))

    def test_raises(self):
        """Test if true."""
        self.assertRaises(TypeError, check_fasta.confirm_fasta,
                          "tests/data/test_experimental_design.txt")


if __name__ == '__main__':
    unittest.main()
