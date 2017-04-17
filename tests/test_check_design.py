#! /usr/bin/env python
"""A suite to test file_exist."""
import sys
import unittest
sys.path.append('../')
from pypiret.Checks import CheckDesign  # noqa

test_exp_file = "tests/data/test_experimental_design.txt"
test_exp_file_wp = "tests/data/test_experimental_design_wrong_path.txt"
negatest_exp_file = "tests/data/test_experimental_design.csv"


class TestCheckDesign(unittest.TestCase):
    """Unittest testcase."""

    def test_tab_false(self):
        """Test if program SystemExit when tab delimited file is not given."""
        design = CheckDesign("tests/__init__.py")
        with self.assertRaises(SystemExit):
            design.tab()

    def test_tab_true(self):
        """Test if it returns True when tab delimited file is given."""
        design = CheckDesign(test_exp_file)
        self.assertTrue(design.tab())

    def test_header_false(self):
        """Test if program SystemExit when header is not correct."""
        design = CheckDesign(negatest_exp_file)
        with self.assertRaises(SystemExit):
            design.header()

    def test_header_true(self):
        """Test if True header is returned."""
        design = CheckDesign(test_exp_file)
        self.assertTrue(design.header())

    def test_sample_name_true(self):
        """Test if true."""
        designs = CheckDesign(test_exp_file)
        self.assertTrue(designs.sample_name())

    def test_file_name_true(self):
        """Test if true."""
        design = CheckDesign(test_exp_file)
        self.assertTrue(design.file_name())

    def test_group_name_true(self):
        """Test if true."""
        group_name = CheckDesign(test_exp_file)
        self.assertTrue(group_name.group_name())

    def test_file_sep(self):
        """Test if correctly generated."""
        fastqs = CheckDesign(test_exp_file)
        self.assertEqual(fastqs.extract_sample_fastqs()['samp1'],
                         "tests/data/fastqs/BTT_test15_R1.fastq:"
                         "tests/data/fastqs/BTT_test15_R2.fastq")

    def test_multi(self):
        """Test if correctly generated."""
        get_fastqs = CheckDesign(test_exp_file)
        self.assertEqual(get_fastqs.extract_sample_fastqs()['samp2'],
                         ["tests/data/fastqs/BTT_test25_R1.fastq:"
                          "tests/data/fastqs/BTT_test25_R2.fastq",
                          "tests/data/fastqs/BTT_test22_R1.fastq:"
                          "tests/data/fastqs/BTT_test22_R2.fastq"])

    def test_fastq_notexists(self):
        """Test if Fastq files in experimental design are not present."""
        see_fastqs = CheckDesign(test_exp_file_wp)
        with self.assertRaises(SystemExit):
            see_fastqs.fastq_exists()

    def test_fastq_exists(self):
        """Test if Fastq files in experimental design are present."""
        see_fastqs = CheckDesign(test_exp_file)
        self.assertTrue(see_fastqs.fastq_exists())

    def test_sample_suff(self):
        """Test if there are enough samples, when both method is chosen."""
        design = CheckDesign(test_exp_file)
        self.assertTrue(design.sample_suff(method='both'))

# class Test
if __name__ == '__main__':
    unittest.main()
