#! /usr/bin/env python

"""A suite to test file_exist."""
import unittest
import sys
sys.path.append('../')
from pypiret import file_exist # noqa


class TestFileExist(unittest.TestCase):
    """Unittest testcase."""

    def test_true(self):
        """Test if true."""
        self.assertTrue(file_exist("__init__.py"))

    def test_fail(self):
        """Test if false."""
        self.assertFalse(file_exist("not_present.py"))

if __name__ == '__main__':
    unittest.main()
