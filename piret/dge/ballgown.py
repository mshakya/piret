#! /usr/bin/env python

"""Check design."""
import os
import sys
import luigi
import shutil
from luigi import LocalTarget
from luigi.util import inherits, requires
import pandas as pd
DIR = os.path.dirname(os.path.realpath(__file__))
script_dir = os.path.abspath(os.path.join(DIR, "../../scripts"))
os.environ["PATH"] += ":" + script_dir
sys.path.insert(0, script_dir)
from plumbum.cmd import Rscript, plot_pathway
import logging


class ballgown(luigi.Task):
    """Find DGE using ballgown."""
    # gff_file=luigi.Parameter()

    exp_design = luigi.Parameter()
    p_value = luigi.FloatParameter()
    workdir = luigi.Parameter()
    kingdom = luigi.Parameter()

    def output(self):
        """Expected output of DGE using edgeR."""
        bg_rdir = os.path.join(self.workdir, "processes", "ballgown", self.kingdom, "migun")
        return LocalTarget(bg_rdir)

    def run(self):
        """Run ballgown."""
        bg_dir = os.path.join(self.workdir, "processes", "ballgown", self.kingdom)
        bg_results = os.path.join(self.workdir, "processes", "ballgown", self.kingdom)
        if os.path.isdir(bg_results) is False:
            os.makedirs(bg_results)

        for name in ["gene", "transcript"]:
            bg_list = [os.path.join(script_dir, "ballgown.R"), "-i", bg_dir, "-e", self.exp_design,
                       "-n", name, "-p", self.p_value,
                       "-o", bg_results]
            bg_cmd = Rscript[bg_list]
            logger = logging.getLogger('luigi-interface')
            logger.info(bg_cmd)
            bg_cmd()

        # self.summ_summ()

    # def summ_summ(self):
    #     """Summarize the summary table to be displayed in edge"""
    #     deseq2_dir = self.workdir + "/DESeq2/" + self.kingdom
    #     all_files = os.listdir(deseq2_dir)
    #     out_file = os.path.join(deseq2_dir, "summary_updown.csv")
    #     summ_files = [pd.read_csv(os.path.join(deseq2_dir, file),
    #                               index_col=0) for file in all_files if "summary.csv" in file ]
    #     summ_df = pd.concat(summ_files)
    #     summ_df.to_csv(out_file)

    # def program_environment(self):
    #     """Environmental variables for this program."""
    #     scriptdir = os.path.join(self.bindir, "/../scripts/")
    #     return {'PATH': scriptdir + ":" + os.environ["PATH"]}
