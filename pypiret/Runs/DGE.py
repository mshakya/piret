#! /usr/bin/env python

"""Check design."""
from __future__ import print_function
import os
import luigi
from luigi import LocalTarget
from pypiret import Summ
from luigi.util import inherits, requires
import pandas as pd
from plumbum.cmd import Rscript



@requires(Summ.FeatureCounts)
class edgeR(luigi.Task):
    """Find DGE using edgeR."""

    exp_design = luigi.Parameter()
    p_value = luigi.FloatParameter()

    def output(self):
        """Expected output of DGE using edgeR."""
        fcount_dir = self.workdir + "/featureCounts"
        edger_dir = self.workdir + "/edgeR/" + self.kingdom
        for root, dirs, files in os.walk(fcount_dir):
            for file in files:
                if file.endswith("csv"):
                    out_filename = file.split(".tsv")[0] + "_RPKM.csv"
                    out_filepath = os.path.join(edger_dir, out_filename)
                    return LocalTarget(out_filepath)

    def run(self):
        """Run edgeR."""
        fcount_dir = self.workdir + "/featureCounts"
        edger_dir = self.workdir + "/edgeR/" + self.kingdom
        edger_location = os.path.join(self.bindir, "../scripts/edgeR.R")
        if not os.path.exists(edger_dir):
            os.makedirs(edger_dir)
        for root, dirs, files in os.walk(fcount_dir):
            for file in files:
                if file.endswith("tsv"):
                    name = file.split("_")[-2]
                    edger_list = [edger_location,
                                  "-r", os.path.join(root, file),
                                  "-e", self.exp_design,
                                  "-p", self.p_value,
                                  "-n", name,
                                  "-o", edger_dir]
                    edger_cmd = Rscript[edger_list]
                    edger_cmd()
        self.summ_summ()

    def summ_summ(self):
        """Summarize the summary table to be displayed in edge"""
        edger_dir = self.workdir + "/edgeR/" + self.kingdom
        all_files = os.listdir(edger_dir)
        out_file = os.path.join(edger_dir, "summary_updown.csv")
        summ_files = [pd.read_csv(os.path.join(edger_dir, file),
                                  index_col=0) for file in all_files if "summary.csv" in file ]
        summ_df = pd.concat(summ_files)
        summ_df.to_csv(out_file)


@requires(Summ.FeatureCounts)
class DESeq2(luigi.Task):
    """Find DGE using edgeR."""

    exp_design = luigi.Parameter()
    p_value = luigi.FloatParameter()

    def output(self):
        """Expected output of DGE using edgeR."""
        fcount_dir = self.workdir + "/featureCounts"
        edger_dir = self.workdir + "/DESeq2/" + self.kingdom
        for root, dirs, files in os.walk(fcount_dir):
            for file in files:
                if file.endswith("csv"):
                    out_filename = file.split(".tsv")[0] + "_FPKM.csv"
                    out_filepath = os.path.join(edger_dir, out_filename)
                    return LocalTarget(out_filepath)

    def run(self):
        """Run edgeR."""
        fcount_dir = self.workdir + "/featureCounts"
        DESeq2_dir = self.workdir + "/DESeq2/" + self.kingdom
        deseq2_location = os.path.join(self.bindir, "../scripts/DESeq2.R")
        if not os.path.exists(DESeq2_dir):
            os.makedirs(DESeq2_dir)
        for root, dirs, files in os.walk(fcount_dir):
            for file in files:
                if file.endswith("tsv"):
                    name = file.split("_")[-2]
                    deseq2_list = [deseq2_location,
                                  "-r", os.path.join(root, file),
                                  "-e", self.exp_design,
                                  "-p", self.p_value,
                                  "-n", name,
                                  "-o", DESeq2_dir]
                    deseq2_cmd = Rscript[deseq2_list]
                    deseq2_cmd()
        self.summ_summ()

    def summ_summ(self):
        """Summarize the summary table to be displayed in edge"""
        deseq2_dir = self.workdir + "/DESeq2/" + self.kingdom
        all_files = os.listdir(deseq2_dir)
        out_file = os.path.join(deseq2_dir, "summary_updown.csv")
        summ_files = [pd.read_csv(os.path.join(deseq2_dir, file),
                                  index_col=0) for file in all_files if "summary.csv" in file ]
        summ_df = pd.concat(summ_files)
        summ_df.to_csv(out_file)

    def program_environment(self):
        """Environmental variables for this program."""
        scriptdir = os.path.join(self.bindir, "/../scripts/")
        return {'PATH': scriptdir + ":" + os.environ["PATH"]}
