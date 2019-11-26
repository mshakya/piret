#! /usr/bin/env python
"""A suite to test file_exist."""
import os
import pytest
import luigi
from luigi.interface import build
from piret.summary import summarize as summ
from plumbum.cmd import rm, cp


def test_summary():
    """Test summary."""
    build([summ.conver2json(gff_file = os.path.join("tests", "data",
                                                    "test_prok.gff"),
                            fasta_file = os.path.join("tests", "data",
                                                     "test_prok.fna"),
                            pathway = True,
                            kingdom = "prokarya",
                            workdir = os.path.join("tests", "data",
                                                     "test_prok"),
                            method = "edgeR",
                            NovelRegions = True)], local_scheduler=True)
    
    assert os.stat(os.path.join("tests", "data", "test_prok",
                                       "prokarya_out.json")).st_size > 1

    rm_cmd = rm["-rf", "tests/data/test_prok/prokarya_out.json"]
    rm_cmd()
    rm_cmd = rm["-rf", "tests/data/test_prok/processes/DESeq2"]
    rm_cmd()
    rm_cmd = rm["-rf", "tests/data/test_prok/processes/edgeR"]
    rm_cmd()