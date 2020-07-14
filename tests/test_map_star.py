#! /usr/bin/env python
"""A suite to test file_exist."""
import os
import shutil
from luigi.interface import build
from piret.maps.star import STARindex
from piret.maps.star import map_star


def test_star_multindex_prok():
    """Test star index creation with multiple fasta for prokarytoes."""
    star_dir = os.path.join("tests", "test_multistar_prok")
    build([STARindex(fasta="tests/data/test_prok.fna,tests/data/test_euk.fna",
                     num_cpus=1,
                     gff_file="tests/data/test_prok.gff, tests/data/test_euk.gff",
                     stardb_dir=star_dir,
                     kingdom="prokarya")], local_scheduler=True)
    assert os.path.exists(os.path.join(star_dir, "SAindex")) is True


def test_star_index_prok():
    """Test star index creation."""
    star_dir = os.path.join("tests", "test_star_prok")
    build([STARindex(fasta="tests/data/test_prok.fna",
                     num_cpus=1,
                     gff_file="tests/data/test_prok.gff",
                     stardb_dir=star_dir,
                     kingdom="prokarya")], local_scheduler=True)
    assert os.path.exists(os.path.join(star_dir, "SAindex")) is True
    shutil.rmtree(star_dir)


def test_star_multindex_euk():
    """Test star index creation with multiple fasta for euks."""
    star_dir = os.path.join("tests", "test_multistar_euk")
    build([STARindex(fasta="tests/data/test_prok.fna,tests/data/test_euk.fna",
                     num_cpus=1,
                     gff_file="tests/data/test_prok.gff, tests/data/test_euk.gff",
                     stardb_dir=star_dir,
                     kingdom="eukarya")], local_scheduler=True)
    assert os.path.exists(os.path.join(star_dir, "SAindex")) is True
    shutil.rmtree(star_dir)


def test_star_index_euk():
    """Test star index creation for euks."""
    star_dir = os.path.join("tests", "test_star_euk")
    build([STARindex(fasta="tests/data/test_prok.fna",
                     num_cpus=1,
                     gff_file="tests/data/test_prok.gff",
                     stardb_dir=star_dir,
                     kingdom="eukarya")], local_scheduler=True)
    assert os.path.exists(os.path.join(star_dir, "SAindex")) is True
    shutil.rmtree(star_dir)


def test_star_multindex_both():
    """Test star index creation with multiple fasta for both."""
    star_dir = os.path.join("tests", "test_star_both")
    build([STARindex(fasta="tests/data/test_prok.fna,tests/data/test_euk.fna",
                     num_cpus=1,
                     gff_file="tests/data/test_prok.gff, tests/data/test_euk.gff",
                     stardb_dir=star_dir,
                     kingdom="both")], local_scheduler=True)
    assert os.path.exists(os.path.join(star_dir, "SAindex")) is True
    shutil.rmtree(star_dir)


def test_star_mapping():
    """Test star index mapping."""
    star_dir = os.path.join("tests", "test_multistar_prok")
    map_dir = os.path.join("tests", "test_star_mapping", "map_dir")
    if os.path.exists(star_dir) is False:
        os.makedirs(star_dir)
    build([map_star(fastqs=["tests/data/BTT_test15_R1.1000.fastq",
                            "tests/data/BTT_test15_R2.1000.fastq"],
                    stardb_dir=star_dir,
                    map_dir=map_dir,
                    sample="samp5",
                    num_cpus=2,
                    align_intron_max=100,
                    align_intron_min=2)], local_scheduler=True)
    assert os.path.exists(os.path.join(map_dir, "samp5_srt.bam")) is True
    shutil.rmtree(star_dir)
