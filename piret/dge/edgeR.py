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
from plumbum.cmd import RDESeq2, gage_analysis, ballgown_analysis
import logging



# @inherits(Summ.FeatureCounts)
class edgeR(luigi.Task):
    """Find DGE using edgeR."""
    kingdom = luigi.Parameter()
    workdir = luigi.Parameter()
    exp_design = luigi.Parameter()
    p_value = luigi.FloatParameter()
    # GAGE = luigi.BoolParameter()
    # pathway = luigi.BoolParameter()
    # gff_file = luigi.Parameter()

    def output(self):
        """Expected output of DGE using edgeR."""
        edger_dir = os.path.join(self.workdir, "processes", "edgeR", self.kingdom)
        out_filepath = os.path.join(edger_dir, "summary_updown.csv")
        return LocalTarget(out_filepath)

    def run(self):
        """Run edgeR."""
        fcount_dir = os.path.join(self.workdir, "processes", "featureCounts", self.kingdom)
        edger_dir = os.path.join(self.workdir, "processes", "edgeR", self.kingdom)
        if not os.path.exists(edger_dir):
            os.makedirs(edger_dir)
        for file in os.listdir(fcount_dir):
            if file.endswith("tsv"):
                name  = file.split("_")[-2]
                edger_list = [os.path.join(script_dir, "EdgeR"), "-r",
                              os.path.join(fcount_dir, file),
                            "-e", self.exp_design,
                            "-p", self.p_value,
                            "-n", name,
                            "-o", edger_dir]
                edger_cmd = Rscript[edger_list]
                logger = logging.getLogger('luigi-interface')
                logger.info(edger_cmd)
                edger_cmd()
                # if file == "gene_count.tsv":
                # # TODO:convert the first column to locus tag
                #     if self.pathway is True:
                #         path_list = ["-d", edger_dir,
                #                      "-m", "edgeR", "-c", self.org_code] # get pathway information
                #         path_cmd = plot_pathway[path_list]
                #         logger.info(path_cmd)
                #         path_cmd()
                #     if self.GAGE is True:
                #         gage_list = ["-d", edger_dir, "-m",
                #                      "edgeR", "-c", self.org_code]
                #         gage_cmd = gage_analysis[gage_list]
                #         logger.info(gage_cmd)
                #         gage_cmd()
        self.summ_summ()

    def summ_summ(self):
        """Summarize the summary table to be displayed in edge"""
        edger_dir = os.path.join(self.workdir, "processes", "edgeR", self.kingdom)
        all_dirs = os.listdir(edger_dir)
        if all_dirs:
            out_file = os.path.join(edger_dir, "summary_updown.csv")
            summ_files = []
            for root, dirs, files in os.walk(edger_dir):
                for file in files:
                    if "summary.csv" in file:
                        summ_files.append(os.path.join(root, file))
            summ_df = [pd.read_csv(file, index_col=0) for file in summ_files]
            summ_df_ccat = pd.concat(summ_df)
            summ_df_ccat.to_csv(out_file)