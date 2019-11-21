#! /usr/bin/env python
"""A suite to test file_exist."""
import os
import pytest
import luigi
from luigi.interface import build
from piret.counts import featurecounts as fc
from plumbum.cmd import rm, cp


def test_featurecount():
    """Test star index creation and mapping."""

    map_dir = os.path.join("tests/test_count", "processes", "mapping", "samp5")
    if os.path.exists(map_dir) is False:
        os.makedirs(map_dir)
    cp_cmd = ["tests/data/samp5_srt.bam", map_dir]
    cp[cp_cmd]()
    
    build([fc.FeatureCounts(fastq_dic={'samp5':''},
                            kingdom="prokarya",
                            gff_file="tests/data/test_prok.gff",
                            workdir="tests/test_count",
                            indexfile="",
                            num_cpus=2,
                            ref_file="tests/data/test_prok.fna",
                            fid="ID",
                            stranded=0)],
                local_scheduler=True)
    assert os.path.exists("tests/test_count/processes/featureCounts/prokarya/gene_count.tsv") is True

    rm_cmd = rm["-rf", "tests/test_count"]
    rm_cmd()
