#! /usr/bin/env python

"""Check design."""
from __future__ import print_function
import os
import luigi
from luigi import LocalTarget
from pypiret import Summ
from luigi.util import inherits, requires
import pandas as pd
from plumbum.cmd import EdgeR, Rscript, plot_pathway, gage_analysis, ballgown_analysis
import logging


@requires(Summ.FeatureCounts)
class edgeR(luigi.Task):
    """Find DGE using edgeR."""
    exp_design = luigi.Parameter()
    p_value = luigi.FloatParameter()
    org_code = luigi.Parameter()

    def output(self):
        """Expected output of DGE using edgeR."""
        fcount_dir = os.path.join(self.workdir, "featureCounts", self.kingdom)
        edger_dir = os.path.join(self.workdir, "edgeR", self.kingdom)
        for root, dirs, files in os.walk(edger_dir):
            for file in files:
                if file.endswith("__sig.csv"):
                    out_folder = file.split(".csv")[0]
                    out_filepath = os.path.join(edger_dir, out_folder, "greater.csv")
                    return LocalTarget(out_filepath)

    def run(self):
        """Run edgeR."""
        fcount_dir = os.path.join(self.workdir, "featureCounts", self.kingdom)
        edger_dir = os.path.join(self.workdir, "edgeR", self.kingdom)
        if not os.path.exists(edger_dir):
            os.makedirs(edger_dir)
        for root, dirs, files in os.walk(fcount_dir):
            for file in files:
                if file.endswith("tsv"):
                    name = file.split("_")[-2]
                    edger_list = ["-r", os.path.join(root, file),
                                  "-e", self.exp_design,
                                  "-p", self.p_value,
                                  "-n", name,
                                  "-o", edger_dir]
                    edger_cmd = EdgeR[edger_list]
                    logger = logging.getLogger('luigi-interface')
                    logger.info(edger_cmd)
                    edger_cmd()
                    if file == "gene_count.tsv":
                        path_list = ["-d", edger_dir,
                            "-m", "edgeR", "-c", self.org_code] # get pathway information
                        path_cmd = plot_pathway[path_list]
                        logger.info(path_cmd)
                        path_cmd()
                        gage_list = ["-d", edger_dir,
                            "edgeR", "-c", self.org_code]
                        gage_cmd = gage_analysis[gage_list]
                        logger.info(gage_cmd)
                        gage_cmd()
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
    # bindir = luigi.Parameter()
    org_code = luigi.Parameter()

    def output(self):
        """Expected output of DGE using edgeR."""
        fcount_dir = os.path.join(self.workdir, "featureCounts", self.kingdom)
        DESeq2_dir = os.path.join(self.workdir, "DESeq2", self.kingdom)
        for root, dirs, files in os.walk(fcount_dir):
             for file in files:
                if file.endswith("__sig.csv"):
                    out_folder = file.split(".csv")[0]
                    out_filepath = os.path.join(DESeq2_dir, out_folder, "greater.csv")
                    
                    return LocalTarget(out_filepath)

    def run(self):
        """Run edgeR."""
        fcount_dir = os.path.join(self.workdir, "featureCounts", self.kingdom)
        DESeq2_dir = os.path.join(self.workdir, "DESeq2", self.kingdom)

        if not os.path.exists(DESeq2_dir):
            os.makedirs(DESeq2_dir)
        for root, dirs, files in os.walk(fcount_dir):
            for file in files:
                if file.endswith("tsv"):
                    name = file.split("_")[-2]
                    deseq2_list = ["DESeq2.R",
                                   "-r", os.path.join(root, file),
                                   "-e", self.exp_design,
                                   "-p", self.p_value,
                                   "-n", name,
                                   "-o", DESeq2_dir]
                    deseq2_cmd = Rscript[deseq2_list]
                    deseq2_cmd()
                if file == "gene_count.tsv":
                    path_list = ["plot_pathway.R", "-d", DESeq2_dir,
                         "-m", "DESeq2", "-c",
                         self.org_code] # get pathway information
                    path_cmd = Rscript[path_list]
                    path_cmd()
                    gage_list = ["gage_analysis.R", "-d", DESeq2_dir,
                         "DESeq2", "-c", self.org_code]
                    gage_cmd = Rscript[gage_list]

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



@inherits(Summ.ReStringTieScoresW)
class ballgown(luigi.Task):
    """Find DGE using ballgown."""

    exp_design = luigi.Parameter()
    p_value = luigi.FloatParameter()

    def output(self):
        """Expected output of DGE using edgeR."""
        bg_rdir = os.path.join(self.workdir, "bg_results", self.kingdom)
        return LocalTarget(bg_rdir)

    def run(self):
        """Run ballgown."""
        bg_dir = os.path.join(self.workdir, "ballgown", self.kingdom)
        bg_results = os.path.join(self.workdir, "bg_results", self.kingdom)
        if os.path.isdir(bg_results) is False:
            os.makedirs(bg_results)

        for name in ["gene", "transcript"]:
            bg_list = ["-i", bg_dir, "-e", self.exp_design,
                       "-n", name, "-p", self.p_value,
                       "-o", bg_results]
            bg_cmd = ballgown_analysis[bg_list]
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
