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
from pypiret.Runs import Map, FastQC  # noqa


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
            lines=pg.readlines()
        self.assertEqual(lines[5].split("\t")[4], "410")
        with open("tests/test_gff2gtf/eukarya_test.gtf") as pg:
            lines=pg.readlines()
        self.assertEqual(lines[5].split("\t")[4], "3007")

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
        luigi.interface.build([Map.GFF2GTF(gff_file="tests/data/test_prok1.1.gff,tests/data/eukarya_test.gff3",
                                           bindir="bin",
                                           workdir="tests/test_createsplice"),
                               Map.CreateSplice(gff_file="tests/data/test_prok1.1.gff,tests/data/eukarya_test.gff3",
                                                bindir="bin",
                                                workdir="tests/test_createsplice")],
                              local_scheduler=True)

        self.assertTrue(os.path.exists("tests/test_createsplice/test_prok1.1.splice"))
        self.assertTrue(os.path.exists("tests/test_createsplice/eukarya_test.splice"))


    def tearDown(self):
        """Remove created files and directories."""
        shutil.rmtree("tests/test_createsplice/")


class TestHisatMapping(unittest.TestCase):
    """Unittest testcase."""

    def setUp(self):
        """Setting up the working directory."""
        if os.path.exists("tests/test_hisatmapping") is False:
            os.makedirs("tests/test_hisatmapping")

    def tearDown(self):
        # """Remove created files and directories."""
        shutil.rmtree("tests/test_hisatmapping/")

    def test_hisat(self):
        """Test hisat index creation, FastQC, and mapping."""
        dir_path = os.path.dirname(os.path.realpath(__file__))
        bindir = os.path.abspath(os.path.join(dir_path, '..', 'bin'))
        luigi.interface.build([
            Map.HisatIndex(fasta="tests/data/test_prok.fa",
                           hi_index="tests/test_hisatmapping/prok_index",
                           bindir="bin",
                           numCPUs=1),
            FastQC.PairedRunQC(fastqs=["tests/data/BTT_test15_R1.1000.fastq",
                                       "tests/data/BTT_test15_R2.1000.fastq"],
                               sample="samp1",
                               numCPUs=1,
                               outdir="tests/test_hisatmapping/trimming_results",
                               bindir=bindir,
                               faqc_min_L=50,
                               n_cutoff=5),
            Map.Hisat(fastq1="tests/data/BTT_test15_R1.1000.fastq",
                      fastq2="tests/data/BTT_test15_R2.1000.fastq",
                      numCPUs=1,
                      indexfile="tests/test_hisatmapping/prok_index",
                      spliceFile="",
                      mappingLogFile="tests/test_hisatmapping/mapping.log",
                      unalned="tests/test_hisatmapping//unligned.fastq",
                      outsam="tests/test_hisatmapping/samp1.sam",
                      ref_file="tests/data/test_prok.fa",
                      bindir=bindir)], local_scheduler=True)

        self.assertTrue(os.path.exists("tests/test_hisatmapping/prok_index.8.ht2l"))
        self.assertTrue(os.path.exists("tests/test_hisatmapping/samp1.sam"))

# class Test
if __name__ == '__main__':
    unittest.main()
