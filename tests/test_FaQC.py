#! /usr/bin/env python


import sys
import os
import unittest
import shutil
from luigi.interface import build
from piret.qc import FaQC
DIR = os.path.dirname(os.path.realpath(__file__))
LIB = os.path.abspath(os.path.join(DIR, '..'))
sys.path.append(LIB)


def test_faqc_function():
    """Test faqc."""
    fq_path = os.path.join("test_faqc", "processes", "qc", "samp1")
    if os.path.exists(fq_path) is False:
        os.makedirs(fq_path)
    build([FaQC.PairedRunQC(fastqs=["tests/data/fastqs/BTT_test15_R1.fastq.gz",
                                    "tests/data/fastqs/BTT_test15_R2.fastq.gz"],
                            sample="samp1",
                            num_cpus=1,
                            workdir="test_faqc",
                            faqc_min_L=50,
                            avg_q=20,
                            n_cutoff=5)],
          local_scheduler=True, workers=1)

    assert os.stat(os.path.join("test_faqc", "processes", "qc", "samp1",
                                "samp1.1.trimmed.fastq")).st_size > 1

    # rm_cmd = rm["-rf", "tests/data/test_prok/processes/DESeq2/prokarya/gene/spleen__over__liver__gene__sig.csv"]
    # rm_cmd()


def test_fq_out():
    """Test if fq has expected line number"""
    fq_lines = file_len(
        "test_faqc/processes/qc/samp1/samp1.1.trimmed.fastq")
    assert(fq_lines) == 177748


def test_output():
    """Test if 5 files are generated as output"""
    filenames = [x[2] for x in os.walk("test_faqc/processes/qc/samp1")][0]
    assert len(filenames) == 5
    shutil.rmtree('test_faqc')


def file_len(fname):
    """A function to test filename"""
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1
