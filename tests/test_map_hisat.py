#! /usr/bin/env python
"""A suite to test file_exist."""

import sys
import os
import unittest
import shutil
import luigi
from luigi.interface import build
from plumbum.cmd import rm
from piret.maps import hisat2
DIR = os.path.dirname(os.path.realpath(__file__))
LIB = os.path.abspath(os.path.join(DIR, '..'))
sys.path.append(LIB)


def test_hisat2_indexing():
    """Test hisat2."""
    hisat2_path = os.path.join("test_hisat2")
    if os.path.exists(hisat2_path) is False:
        os.makedirs(hisat2_path)
    build([hisat2.HisatIndex(fasta="tests/data/test_prok.fna",
                             hi_index="test_hisat2/prok_index",
                             num_cpus=1)],
          local_scheduler=True, workers=1)

    assert os.stat(os.path.join("test_hisat2", "prok_index.8.ht2")).st_size > 1


def test_hisat2_multindexing():
    """Test hisat2."""
    hisat2_path = os.path.join("test_hisat2")
    if os.path.exists(hisat2_path) is False:
        os.makedirs(hisat2_path)
    build([hisat2.HisatIndex(fasta="tests/data/test_prok.fna,tests/data/test_euk.fna",
                             hi_index="test_hisat2/prok_index",
                             num_cpus=1)],
          local_scheduler=True, workers=1)

    assert os.stat(os.path.join("test_hisat2", "prok_index.8.ht2")).st_size > 1


def test_hisat2_mapping():
    """Test hisat2 mapping."""
    hisat2_path = os.path.join("test_hisat2", "processes", "mapping", "samp1")
    if os.path.exists(hisat2_path) is False:
        os.makedirs(hisat2_path)
    map_dir = os.path.join("test_hisat2", "processes", "mapping")
    sam_file = os.path.join(hisat2_path, "samp1.sam")
    build([hisat2.Hisat(fastqs=["tests/data/BTT_test15_R1.1000.fastq",
                                "tests/data/BTT_test15_R2.1000.fastq"],
                        min_introlen=500,
                        max_introlen=5000,
                        rna_strandness=0,
                        kingdom="prokarya",
                        num_cpus=1,
                        sample="samp1",
                        # map_dir=map_dir,
                        indexfile="test_hisat2/prok_index",
                        # outsam=os.path.join(hisat2_path, "samp1.sam"),
                        workdir="test_hisat2")], local_scheduler=True)
    assert os.path.exists(sam_file)
    shutil.rmtree('test_hisat2')
