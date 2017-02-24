#! /usr/bin/env python
"""A suite to test file_exist."""
import sys
import unittest
sys.path.append('../')
from pypiret import CheckDesign  # noqa


class TestCheckTab(unittest.TestCase):
    """Unittest testcase."""

    def test_false(self):
        """Test if false."""
        check_design_csv = CheckDesign("test_experimental_design.csv")
        self.assertFalse(check_design_csv.check_tab)

    def test_true(self):
        """Test if true."""
        check_design_tsv = CheckDesign("test_experimental_design.txt")
        self.assertTrue(check_design_tsv.check_tab)


if __name__ == '__main__':
    unittest.main()
