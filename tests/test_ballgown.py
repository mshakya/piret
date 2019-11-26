#! /usr/bin/env python
"""A suite to test file_exist."""
import os
import pytest
import luigi
from luigi.interface import build
from piret.dge import ballgown
from piret.counts import stringtie
from plumbum.cmd import rm, cp



def test_restringtie():
    """Test restringtie approach."""
    build([stringtie.ReStringTieScoresW(
                        fastq_dic={'samp1':"l",
                                   'samp2':['b'],
                                   'samp3':['b'],
                                   'samp4':['b'],
                                   'samp5':['b'],
                                   'samp6':['b']},
                        num_cpus=2,
                        workdir="tests/data/test_prok/",
                        kingdom="prokarya")], local_scheduler=True)
assert os.stat(os.path.join("tests", "data", "test_prok", "processes",
                                       "ballgown", "prokarya", "samp6",
                                       "e_data.ctab" )).st_size > 1    


def test_ballgown():
    """Test deseq3."""
    build([ballgown.ballgown(workdir=os.path.join("tests", "data",
                                                     "test_prok"),
                         kingdom="prokarya",
                         exp_design="tests/test_prok.txt",
                         p_value=0.1)], local_scheduler=True)
    
    assert os.stat(os.path.join("tests", "data", "test_prok", "processes",
                                       "ballgown", "prokarya",
                                       "liver__spleen__gene__sig.csv" )).st_size > 1

    # rm_cmd = rm["-rf", "tests/data/test_prok/processes/DESeq2/prokarya/gene/spleen__over__liver__gene__sig.csv"]
    # rm_cmd()
    