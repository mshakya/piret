"""Import important functions."""

import sys
import os
import shutil
import subprocess
dir_path = os.path.dirname(os.path.realpath(__file__))
script_dir = os.path.abspath(os.path.join(dir_path, "../scripts"))
sys.path.append(dir_path)
os.environ["PATH"] += ":" + script_dir
