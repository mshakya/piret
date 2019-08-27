"""Importing."""
import sys, os
dir_path = os.path.dirname(os.path.realpath(__file__))

sys.path.append(dir_path)
from design import CheckDesign
from fasta import CheckFasta
from dependencies import CheckDependencies
