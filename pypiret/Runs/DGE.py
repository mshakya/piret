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
    bindir = luigi.Parameter()

    def output(self):
        """Expected output of DGE using edgeR."""
        fcount_dir = os.path.join(self.workdir, "featureCounts", self.kingdom)
        edger_dir = os.path.join(self.workdir, "edgeR", self.kingdom)
        for root, dirs, files in os.walk(fcount_dir):
            for file in files:
                if file.endswith("csv"):
                    out_filename = file.split(".tsv")[0] + "_RPKM.csv"
                    out_filepath = os.path.join(edger_dir, out_filename)
                    return LocalTarget(out_filepath)

    def run(self):
        """Run edgeR."""
        fcount_dir = os.path.join(self.workdir, "featureCounts", self.kingdom)
        edger_dir = os.path.join(self.workdir, "edgeR", self.kingdom)
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
        if all_files:
            out_file = os.path.join(edger_dir, "summary_updown.csv")
            summ_files = [pd.read_csv(os.path.join(edger_dir, file),
                                  index_col=0) for file in all_files if "summary.csv" in file ]
            summ_df = pd.concat(summ_files)
            summ_df.to_csv(out_file)


@requires(Summ.FeatureCounts)
class DESeq2(luigi.Task):
    """Find DGE using DESeq2."""

    exp_design = luigi.Parameter()
    p_value = luigi.FloatParameter()
    bindir = luigi.Parameter()

    def output(self):
        """Expected output of DGE using edgeR."""
        fcount_dir = os.path.join(self.workdir, "featureCounts", self.kingdom)
        DESeq2_dir = os.path.join(self.workdir, "DESeq2", self.kingdom)
        for root, dirs, files in os.walk(fcount_dir):
            for file in files:
                if file.endswith("csv"):
                    out_filename = file.split(".tsv")[0] + "_FPKM.csv"
                    out_filepath = os.path.join(DESeq2_dir, out_filename)
                    return LocalTarget(out_filepath)

    def run(self):
        """Run edgeR."""
        fcount_dir = os.path.join(self.workdir, "featureCounts", self.kingdom)
        DESeq2_dir = os.path.join(self.workdir, "DESeq2", self.kingdom)
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


@inherits(Summ.ReStringTieScoresW)
class ballgown(luigi.Task):
    """Find DGE using ballgown."""

    exp_design = luigi.Parameter()
    p_value = luigi.FloatParameter()
    bindir = luigi.Parameter()

    def output(self):
        """Expected output of DGE using edgeR."""
        bg_rdir = os.path.join(self.workdir, "bg_results", self.kingdom)
        return LocalTarget(bg_rdir)

    def run(self):
        """Run ballgown."""
        # Rscript scripts/ballgown.R -i tests/test_euk/ballgown/ -e test_euk.txt -o test_ballgown -n exon
        bg_dir = os.path.join(self.workdir, "ballgown", self.kingdom)
        bg_results = os.path.join(self.workdir, "bg_results", self.kingdom)
        if os.path.isdir(bg_results) is False:
            os.makedirs(bg_results)
        bg_loc = os.path.join(self.bindir, "../scripts/ballgown.R")

        for name in ["gene", "transcript"]:
            bg_list = [bg_loc, "-i", bg_dir, "-e", self.exp_design,
                       "-n", name, "-p", self.p_value,
                       "-o", bg_results]
            bg_cmd = Rscript[bg_list]
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

    def program_environment(self):
        """Environmental variables for this program."""
        scriptdir = os.path.join(self.bindir, "/../scripts/")
        return {'PATH': scriptdir + ":" + os.environ["PATH"]}
