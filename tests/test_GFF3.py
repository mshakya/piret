#! /usr/bin/env python
"""A suite to test confirm_fasta."""
import os
import sys
import pytest
from plumbum.cmd import rm
DIR = os.path.dirname(os.path.realpath(__file__))
LIB = os.path.abspath(os.path.join(DIR, '..'))
sys.path.append(LIB)
from piret.checks.GFF3 import CheckGFF


def test_size():
    """Test if gff has 9 columns."""
    size_gff = CheckGFF("tests/data/test_prok.gff")
    result = size_gff.size()
    assert result is True

def test_check_id():
    """Test if gff has ID in all rows of 9th column."""
    gff = CheckGFF("tests/data/test_prok.gff")
    check_id_result = gff.check_id()
    assert check_id_result is True

def test_checkuniqueid():
    """Test if gff has unique IDs."""
    gff = CheckGFF("tests/data/eukarya_test.gff3")
    ck_out = gff.check_unique_id()
    assert ck_out is True
