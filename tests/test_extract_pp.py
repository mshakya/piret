#! /usr/bin/env python
"""A suite to test confirm_fasta."""
import os
import sys
import unittest
from plumbum.cmd import rm
dir_path = os.path.dirname(os.path.realpath(__file__))
lib_path = os.path.abspath(os.path.join(dir_path, '..'))
sys.path.append(lib_path)
from pypiret.Runs.srna import ExtractPP
from pypiret.Runs.srna import FindNovelRegions, CompileGFF


class TestPP(unittest.TestCase):
    """Unittest testcase."""
    bam_file = "tests/data/samp5.bam"

    def test_split(self):
        """Test if true."""
        ExtractPP.prop_paired(self, "samp5", self.bam_file)
        assert os.path.exists("/tmp/samp5_f163.bam") == 1

    def test_merge(self):
        """Test if true."""

        ExtractPP.merge_prop_paired(self, "samp5", "forward.samp5.pp.bam",
                                    "backward.samp5.pp.bam")
        assert os.path.exists("forward.samp5.pp.bam") == 1
        assert os.path.exists("backward.samp5.pp.bam") == 1

    def test_sort(self):
        """Test if true."""
        ExtractPP.sort_bam(self, "forward.samp5.pp.bam")
        assert os.path.exists("forward.samp5.pp_srt.bam") == 1

    @classmethod
    def tearDownClass(cls):
        rm["-rf", 'forward.samp5.pp.bam']()
        rm["-rf", 'backward.samp5.pp.bam']()
        rm["-rf", 'forward.samp5.pp_srt.bam']()


class TestCov(unittest.TestCase):
    """Unittest testcase."""

    def test_ref(self):
        FindNovelRegions.get_genome_ref(self, "tests/data/samp5.fw_srt.bam",
                                        "test_samp5_size.txt")
        assert os.path.exists("test_samp5_size.txt")

    @classmethod
    def tearDownClass(cls):
        rm["-rf", 'test_samp5_size.txt']()

    # def test_cov(self):
    #     FindNovelRegions.genome_coverage(self, "forward.samp5.pp_srt.bam",
    #                                 "forward.samp5.pp_srt.bam.genome_size",
    #                                 "test.txt")
    #     assert os.path.exists("test.txt")

    # def test_novel(self):
    #     FindNovelRegions.novel_regions(self, "tests/data/test_prok.gff",
    #                                    "samp5_fw_novel.bedfile",
    #                                    "test.txt.novel")
    #     assert os.path.exists("test.txt.novel") == 1


class TestNovelGFF(unittest.TestCase):
    """Unittest testcase."""

    def test_compile_bed(self):
        CompileGFF.compile_novel_regions(self, ["tests/data/bedfile1.txt",
                                                "tests/data/bedfile2.txt"],
                                         "all_novel.txt")
        assert os.path.exists("all_novel.txt") == 1

    def test_make_gff(self):
        CompileGFF.make_gff(self, "all_novel.txt", "+", "test.gff")
        assert os.path.exists("test.gff") == 1

    @classmethod
    def tearDownClass(cls):
        rm["-rf", 'test.gff']()

if __name__ == '__main__':
    unittest.main()
