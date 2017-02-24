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
        check_design_csv = CheckDesign()
        self.assertFalse(check_design_csv.check_tab("__init__.py"))

    def test_true(self):
        """Test if true."""
        check_design_tsv = CheckDesign()
        self.assertTrue(check_design_tsv.check_tab("test_experimental_design.txt"))


class TestCheckHeader(unittest.TestCase):
    """Unittest testcase."""

    def test_false(self):
        """Test if False."""
        check_design = CheckDesign()
        self.assertFalse(check_design.check_header("test_experimental_design.csv"))

    def test_true(self):
        """Test if True."""
        check_design = CheckDesign()
        self.assertTrue(check_design.check_header("test_experimental_design.txt"))

if __name__ == '__main__':
    unittest.main()
