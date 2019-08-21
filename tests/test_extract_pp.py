#! /usr/bin/env python
"""A suite to test confirm_fasta."""
import os
import sys
import unittest
from plumbum.cmd import rm
DIR = os.path.dirname(os.path.realpath(__file__))
LIB = os.path.abspath(os.path.join(DIR, '..'))
sys.path.append(LIB)
from pypiret.Runs.srna import ExtractPP
from pypiret.Runs.srna import FindNovelRegions, CompileGFF


class TestPP(unittest.TestCase):
    """Unittest testcase."""
    bam_file = "tests/data/samp5.bam"
    workdir = "tests"
    pp = ExtractPP(kingdom="prokarya",
                       workdir="tests/",
                       map_dir="tests/data/",
                       sample="samp5",
                       num_cpus=1)

    def test_split(self):
        """Test if true."""

        self.pp.prop_paired(self.bam_file)

        assert os.path.exists(os.path.join(self.workdir, "tmp", "samp5_f163.bam")) == 1

    def test_merge(self):
        """Test if true."""
        self.pp.prop_paired(self.bam_file)
        self.pp.merge_prop_paired("samp5", "forward.samp5.pp.bam",
                                  "backward.samp5.pp.bam")
        assert os.path.exists("forward.samp5.pp.bam") == 1
        assert os.path.exists("backward.samp5.pp.bam") == 1

    def test_sort(self):
        """Test if true."""

        self.pp.sort_bam("forward.samp5.pp.bam")
        assert os.path.exists("forward.samp5.pp_srt.bam") == 1

    @classmethod
    def tearDownClass(cls):
        rm["-rf", os.path.join("tests", "tmp")]()
        rm["-rf", 'forward.samp5.pp.bam']()
        rm["-rf", 'backward.samp5.pp.bam']()
        rm["-rf", 'forward.samp5.pp_srt.bam']()


class TestCov(unittest.TestCase):
    """Unittest testcase."""
    fn = FindNovelRegions(kingdom="prokarya", workdir="tests/",
                          gff_file="tests/data/test_prok.gff",
                          map_dir="tests/data", sample="samp5")

    def test_ref(self):
        self.fn.get_genome_ref("tests/data/samp5.fw_srt.bam",
                                        "test_samp5_size.txt")
        assert os.path.exists("test_samp5_size.txt")

    @classmethod
    def tearDownClass(cls):
        rm["-rf", 'test_samp5_size.txt']()

    def test_cov(self):
        self.fn.genome_coverage("tests/data/samp5.fw_srt.bam",
                                "samp5.fw_srt.bam.genome_size")
        assert os.path.exists("samp5.fw_srt.bam.genome_size")

    # def test_novel(self):
    #     FindNovelRegions.NovelRegions(self, "tests/data/test_prok.gff",
    #                                    "samp5_fw_novel.bedfile",
    #                                    "test.txt.novel")
    #     assert os.path.exists("test.txt.novel") == 1


class TestNovelGFF(unittest.TestCase):
    """Unittest testcase."""
    cg = CompileGFF(kingdom="prokarya", workdir="tests/",
                    gff_file="tests/data/test_prok.gff",
                    fastq_dic={})

    def test_compile_bed(self):
        self.cg.compile_NovelRegions(["tests/data/bedfile1.txt",
                                       "tests/data/bedfile2.txt"],
                                       "all_novel.txt")
        assert os.path.exists("all_novel.txt") == 1

    def test_make_gff(self):
        self.cg.make_gff("all_novel.txt", "+", "test.gff")
        assert os.path.exists("test.gff") == 1

    @classmethod
    def tearDownClass(cls):
        rm["-rf", 'test.gff']()
        rm["-rf", 'all_novel.txt']()

if __name__ == '__main__':
    unittest.main()
