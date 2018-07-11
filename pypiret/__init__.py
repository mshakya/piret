"""Import important functions."""

import sys, os
dir_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append(dir_path)

from Checks.Design import CheckDesign
from Checks.Fasta import CheckFasta
from Runs import FaQC
from Runs import srna
from Runs import Map
from Runs import Summ
from Runs import DGE
from Checks.Dependencies import CheckDependencies
from initialize import Initialize
