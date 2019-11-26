#! /usr/bin/env python
"""A suite to test file_exist."""
import os
import pytest
import luigi
from luigi.interface import build
from piret.dge import DESeq2
from plumbum.cmd import rm, cp


def test_deseq2():
    """Test deseq3."""
    build([DESeq2.DESeq2(workdir=os.path.join("tests", "data",
                                                     "test_prok"),
                         kingdom="prokarya",
                         exp_design="tests/test_prok.txt",
                         p_value=0.1)], local_scheduler=True)
    
    assert os.stat(os.path.join("tests", "data", "test_prok", "processes",
                                       "DESeq2", "prokarya", "gene",
                                       "liver__over__spleen__gene__sig.csv" )).st_size > 1

    rm_cmd = rm["-rf", "tests/data/test_prok/processes/DESeq2/prokarya/gene/spleen__over__liver__gene__sig.csv"]
    rm_cmd()
    