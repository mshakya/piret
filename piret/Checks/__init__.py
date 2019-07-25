"""Importing."""
import sys, os
dir_path = os.path.dirname(os.path.realpath(__file__))

sys.path.append(dir_path)
from Design import CheckDesign
from Fasta import CheckFasta
from Dependencies import CheckDependencies
