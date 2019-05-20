#! /usr/bin/env python
"""A suite to test file_exist."""
import os
import pytest
import luigi
from pypiret.Runs import Summ
from plumbum.cmd import rm



def test_star():
        """Test star index creation and mapping."""
        luigi.interface.build([Summ.FeatureCountsII(gff="tests/data2/chr22_ERCC92.gff3",
                               num_cpus=1, stranded=1, out_dir="tests/test_featurecount",
                      bam_list=["tests/data2/chr22.bam", "tests/data2/chr22.bam"])], local_scheduler=True)
        assert os.path.exists("tests/test_featurecount/gene_count.tsv") is True

        rm_cmd = rm["-rf", "tests/test_star_mapping"]
        rm_cmd()
