#! /usr/bin/env python
"""A suite to test file_exist."""
import os
import shutil
import luigi
from luigi.interface import build
from piret.maps.star import STARindex
from piret.maps.star import map_star
from piret.maps.star import map_starW
from piret.maps.star import SummarizeStarMap


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
    """Test star index creation with just one reference."""
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
    star_dir = os.path.join("tests", "test_star_ind_both")
    build([STARindex(fasta="tests/data/test_prok.fna,tests/data/test_euk.fna",
                     num_cpus=1,
                     gff_file="tests/data/test_prok.gff, tests/data/test_euk.gff",
                     stardb_dir=star_dir,
                     kingdom="both")], local_scheduler=True)
    assert os.path.exists(os.path.join(star_dir, "SAindex")) is True


def test_star_mapping():
    """Test star index mapping."""
    star_dir = os.path.join("tests", "test_star_ind_both")
    map_dir = os.path.join("tests", "test_star_mapping",
                           "processes", "mapping", "samp5")
    build([map_star(fastqs=["tests/data/BTT_test15_R1.1000.fastq",
                            "tests/data/BTT_test15_R2.1000.fastq"],
                    stardb_dir=star_dir,
                    map_dir=map_dir,
                    sample="samp5",
                    num_cpus=2,
                    align_intron_max=100,
                    align_intron_min=2)
           ],
          local_scheduler=True)
    assert os.path.exists(os.path.join(map_dir, "samp5_srt.bam")) is True
    shutil.rmtree("tests/test_star_mapping")


def test_star_mapping_w():
    """Test star index mappingW."""
    star_dir = os.path.join("tests", "test_star_ind_both")
    work_dir = os.path.join("tests", "test_star_mapping")
    fastq_dic = {"samp5": ["tests/data/BTT_test15_R1.1000.fastq",
                           "tests/data/BTT_test15_R2.1000.fastq"],
                 "samp6": ["tests/data/BTT_test15_R1.1000.fastq",
                           "tests/data/BTT_test15_R2.1000.fastq"]}
    build([map_starW(fastq_dic=fastq_dic,
                     stardb_dir=star_dir,
                     workdir=work_dir,
                     num_cpus=1, align_intron_max=100,
                     align_intron_min=2)],
          local_scheduler=True)
    assert os.path.exists(os.path.join(
        work_dir, "processes", "mapping", "samp5", "samp5_srt.bam")) is True
    assert os.path.exists(os.path.join(
        work_dir, "processes", "mapping", "samp6", "samp6_srt.bam")) is True
    shutil.rmtree(star_dir)


def test_star_mapping_summ():
    """Test summary function."""
    fastq_dic = {"samp5": ["tests/data/BTT_test15_R1.1000.fastq",
                           "tests/data/BTT_test15_R2.1000.fastq"]}
    build([SummarizeStarMap(fastq_dic=fastq_dic, workdir="tests/test_star_mapping")],
          local_scheduler=True)
    assert os.path.exists(os.path.join(
        "tests", "test_star_mapping", "processes", "mapping", "MapSummary.csv"))
    shutil.rmtree(os.path.join("tests", "test_star_mapping"))
