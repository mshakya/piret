#! /usr/bin/env python
"""A suite to test file_exist."""
import sys
import unittest
sys.path.append('../')
from pypiret import CheckDesign  # noqa


class TestCheckDesign(unittest.TestCase):
    """Unittest testcase."""

    def test_tab_false(self):
        """Test if false."""
        check_design = CheckDesign()
        self.assertFalse(check_design.check_tab("__init__.py"))

    def test_tab_true(self):
        """Test if true."""
        check_design = CheckDesign()
        self.assertTrue(check_design.check_tab("test_experimental_design.txt"))

    def test_header_false(self):
        """Test if False."""
        check_design = CheckDesign()
        self.assertFalse(check_design.check_header("test_experimental_design.csv"))

    def test_header_true(self):
        """Test if True."""
        check_design = CheckDesign()
        self.assertTrue(check_design.check_header("test_experimental_design.txt"))

    def test_sample_name_true(self):
        """Test if true."""
        check_sample_name = CheckDesign()
        self.assertTrue(check_sample_name.check_sample_name("test_experimental_design.txt"))

    def test_file_name_true(self):
        """Test if true."""
        check_sample_name = CheckDesign()
        self.assertTrue(check_sample_name.check_file_name("test_experimental_design.txt"))

    def test_group_name_true(self):
        """Test if true."""
        check_group_name = CheckDesign()
        self.assertTrue(check_group_name.check_group_name("test_experimental_design.txt"))

    def test_file_sep(self):
        """Test if correctly generated."""
        get_fastqs = CheckDesign()
        self.assertEqual(get_fastqs.extract_sample_fastqs("test_experimental_design.txt")['samp1'],
                         "fastqs/BTT_test15_R1.fastq:fastqs/BTT_test15_R2.fastq")

    def test_multi(self):
        """Test if correctly generated."""
        get_fastqs = CheckDesign()
        self.assertEqual(get_fastqs.extract_sample_fastqs("test_experimental_design.txt")['samp2'],
                         ["fastqs/BTT_test25_R1.fastq:fastqs/BTT_test25_R2.fastq",
                          "fastqs/BTT_test22_R1.fastq:fastqs/BTT_test22_R2.fastq"])


# class Test
if __name__ == '__main__':
    unittest.main()
