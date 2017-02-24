#! /usr/bin/env python

"""
    A suite to test file_exist.

    test
"""

import sys
sys.path.append('../')
import unittest
from pypiret import file_exist


class TestFileExist(unittest.TestCase):
    """
    Unittest testcase.

    a test case
    """

    def test_true(self):
        """Test if true."""
        self.assertTrue(file_exist("__init__.py"))

    def test_fail(self):
        """Test if false."""
        self.assertFalse(file_exist("not_present.py"))

if __name__ == '__main__':
    unittest.main()
