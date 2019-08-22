"""Import important functions."""

import sys
import os
import shutil
import subprocess
dir_path = os.path.dirname(os.path.realpath(__file__))
script_dir = os.path.abspath(os.path.join(dir_path, "../scripts"))
sys.path.append(dir_path)
os.environ["PATH"] += ":" + script_dir
# from Checks.Design import CheckDesign
# from Checks.Fasta import CheckFasta
# from Checks.Dependencies import CheckDependencies
# from Runs import FaQC
# from Runs import srna
# from Runs import Map
# from Runs import Summ
# from workflows import SingleSeq
# from initialize import Initialize