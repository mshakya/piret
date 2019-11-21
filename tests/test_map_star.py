#! /usr/bin/env python
"""A suite to test file_exist."""
import os
import pytest
import luigi
from piret.maps import star
from plumbum.cmd import rm



def test_star():
    """Test star index creation, and mapping."""
    star_dir = os.path.join("tests", "test_star_mapping")
    map_dir = os.path.join("tests", "test_star_mapping", "map_dir")
    if os.path.exists(star_dir) is False:
        os.makedirs(star_dir)
    luigi.interface.build([star.STARindex(fasta="tests/data/test_prok.fna",
                                              num_cpus=1,
                                              gff_file="tests/data/test_prok.gff",
                                              stardb_dir=star_dir,
                                              kingdom="prokarya"),
                            star.map_star(fastqs=["tests/data/fastqs/samp5.1.trimmed.fastq", "tests/data/fastqs/samp5.2.trimmed.fastq"],
                                          stardb_dir=star_dir,
                                          map_dir=map_dir,
                                          sample="samp5",
                                          num_cpus=2,
                                          align_intron_max=100,
                                          align_intron_min=2)], local_scheduler=True)
    assert os.path.exists(os.path.join(map_dir, "samp5_srt.bam")) is True
    
    rm_cmd = rm["-rf", "tests/test_star_mapping"]
    rm_cmd()

