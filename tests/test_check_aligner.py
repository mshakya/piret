#! /usr/bin/env python
from piret.checks import aligner
import os
import pytest
import luigi
from piret.counts import featurecounts as fc
from plumbum.cmd import rm, cp




def test_check_aligner_hisat2():

    aln_fdr = os.path.join("tests", "test_check_aligner", "processes")
    if os.path.exists(aln_fdr) is False:
        os.makedirs(aln_fdr)
    with open(os.path.join(aln_fdr, "hisat_index"), 'w') as f:
        f.write("testdata")
    db = aligner.check_aligner(aligner="hisat2",
                               hisat_index="hisat_index",
                               workdir=os.path.join("tests", "test_check_aligner"),
                               star_db="")
    assert os.path.exists(db) is True
    rm_cmd = rm["-rf", "tests/test_check_aligner"]
    rm_cmd()


def test_check_star():

    aln_fdr = os.path.join("tests", "test_check_aligner", "processes", "stardb")
    if os.path.exists(aln_fdr) is False:
        os.makedirs(aln_fdr)
    with open(os.path.join(aln_fdr, "chrName.txt"), 'w') as f:
        f.write("testdata")
    db = aligner.check_aligner(aligner="star",
                               hisat_index="hisat_index",
                               workdir=os.path.join("tests", "test_check_aligner"),
                               star_db="")
    assert os.path.exists(db) is True
    rm_cmd = rm["-rf", "tests/test_check_aligner"]
    rm_cmd()