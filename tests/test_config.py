#! /usr/bin/env python
"""A suite to test config reader."""
import argparse
import pytest
from piret.checks.config import check_input, parse_config


def test_prok_incomp_fasta():
    """Test if it can detect configuration for a prok run."""
    param_dicts = parse_config("tests/test_prok.cfg")
    param_dicts["fasta_prok"] = "tests/data/test_prok.fna,"
    parser = argparse.ArgumentParser(prog='piret',
                                     description="""piret""")
    with pytest.raises(SystemExit):
        check_input(param_dicts, parser)


def test_prok_comp_fasta():
    """Test if it can detect configuration for a prok run."""
    param_dicts = parse_config("tests/test_prok.cfg")
    param_dicts["fasta_prok"] = "tests/data/test_prok.fna"
    parser = argparse.ArgumentParser(prog='piret',
                                     description="""piret""")
    assert param_dicts["gff_prok"] == "tests/data/test_prok.gff"


def test_prok_incomp_gff():
    """Test if it can detect configuration for a prok run."""
    param_dicts = parse_config("tests/test_prok.cfg")
    param_dicts["fasta_prok"] = "tests/data/test_prok.fna,tests/data/test_prok.fna"
    param_dicts["gff_prok"] = "tests/data/test_prok.gff"
    parser = argparse.ArgumentParser(prog='piret',
                                     description="""piret""")
    with pytest.raises(SystemExit):
        check_input(param_dicts, parser)


def test_euk_incomp_gff():
    """Test if it can detect configuration for a prok run."""
    param_dicts = parse_config("tests/test_prok.cfg")
    param_dicts["kingdom"] = "eukarya"
    param_dicts["fasta_euk"] = "tests/data/test_euk.fna,tests/data/test_euk.fna"
    param_dicts["gff_euk"] = "tests/data/test_euk.gff"
    parser = argparse.ArgumentParser(prog='piret',
                                     description="""piret""")
    with pytest.raises(SystemExit):
        check_input(param_dicts, parser)
