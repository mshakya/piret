#! /usr/bin/env python
"""A suite to test file_exist."""
import sys
import os
import luigi
import unittest
import shutil
dir_path = os.path.dirname(os.path.realpath(__file__))
lib_path = os.path.abspath(os.path.join(dir_path, '..'))
sys.path.append(lib_path)
from pypiret.Runs import Map  # noqa


class TestGFF2GTF(unittest.TestCase):
    """Unittest testcase."""

    def setUp(self):
        """Setting tup the directory."""
        if os.path.exists("tests/test_gff2gtf") is False:
            os.makedirs("tests/test_gff2gtf")

    def test_gtf_creation(self):
        """Test if GFF2GTF works."""
        luigi.interface.build([Map.GFF2GTF(gff_file="tests/data/test_prok.gff",
                                           bindir="bin",
                                           workdir="tests/test_gff2gtf")],
                              local_scheduler=True)

        self.assertTrue(os.path.exists("tests/test_gff2gtf/test_prok.gtf"))

    def test_gtfs_creation(self):
        """Test if GFF2GTF works for multiple gffs."""
        luigi.interface.build([Map.GFF2GTF(gff_file="tests/data/test_prok.gff,tests/data/eukarya_test.gff3",
                                           bindir="bin",
                                           workdir="tests/test_gff2gtf")],
                              local_scheduler=True)

        self.assertTrue(os.path.exists("tests/test_gff2gtf/test_prok.gtf"))
        self.assertTrue(os.path.exists("tests/test_gff2gtf/eukarya_test.gtf"))
        with open("tests/test_gff2gtf/test_prok.gtf") as pg:
            first_line_fourth_column = pg.readline().split("\t")[3]
            self.assertEqual(first_line_fourth_column, "462")
        with open("tests/test_gff2gtf/eukarya_test.gtf") as pg:
            first_line_fifth_column = pg.readline().split("\t")[4]
            self.assertEqual(first_line_fifth_column, "2754")

    def tearDown(self):
        """Remove created files and directories."""
        shutil.rmtree("tests/test_gff2gtf/")


class TestCreateSplice(unittest.TestCase):
    """Unittest testcase."""

    def setUp(self):
        """Setting up the directory."""
        if os.path.exists("tests/test_createsplice") is False:
            os.makedirs("tests/test_createsplice")

    def test_splice_creation(self):
        """Test if CreateSplice works."""
        luigi.interface.build([Map.GFF2GTF(gff_file="tests/data/test_prok.gff",
                                           bindir="bin",
                                           workdir="tests/test_createsplice"),
                               Map.CreateSplice(gff_file="tests/data/test_prok.gff",
                                                bindir="bin",
                                                workdir="tests/test_createsplice")],
                              local_scheduler=True)
        self.assertTrue(os.path.exists("tests/test_createsplice/test_prok.splice"))

    def test_splices_creation(self):
        """Test if GFF2GTF works for multiple gffs."""
        luigi.interface.build([Map.GFF2GTF(gff_file="tests/data/test_prok.gff,tests/data/eukarya_test.gff3",
                                           bindir="bin",
                                           workdir="tests/test_createsplice"),
                               Map.CreateSplice(gff_file="tests/data/test_prok.gff,tests/data/eukarya_test.gff3",
                                                bindir="bin",
                                                workdir="tests/test_createsplice")],
                              local_scheduler=True)

        self.assertTrue(os.path.exists("tests/test_createsplice/test_prok.splice"))
        self.assertTrue(os.path.exists("tests/test_createsplice/eukarya_test.splice"))
        with open("tests/test_createsplice/eukarya_test.splice") as pg:
            first_line_second_column = pg.readline().split("\t")[1]
            self.assertEqual(first_line_second_column, "220")

    # def tearDown(self):
    #     """Remove created files and directories."""
    #     shutil.rmtree("tests/test_createsplice/")

# class Test
if __name__ == '__main__':
    unittest.main()
