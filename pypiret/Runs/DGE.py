#! /usr/bin/env python

"""Check design."""
import os
import sys
import luigi
import shutil
from luigi import LocalTarget
from pypiret.Runs import Summ
from luigi.util import inherits, requires
import pandas as pd
DIR = os.path.dirname(os.path.realpath(__file__))
script_dir = os.path.abspath(os.path.join(DIR, "../../scripts"))
os.environ["PATH"] += ":" + script_dir
sys.path.insert(0, script_dir)
from plumbum.cmd import EdgeR, Rscript, plot_pathway
from plumbum.cmd import RDESeq2, gage_analysis, ballgown_analysis
import logging


@requires(Summ.FeatureCounts)
class edgeR(luigi.Task):
    """Find DGE using edgeR."""
    exp_design = luigi.Parameter()
    p_value = luigi.FloatParameter()
    prok_org_code = luigi.Parameter()
    euk_org_code = luigi.Parameter()
    GAGE = luigi.BoolParameter()
    pathway = luigi.BoolParameter()
    gff_file = luigi.Parameter()

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
        for file in os.listdir(fcount_dir):
            if file.endswith("tsv"):
                name  = file.split("_")[-2]
                edger_list = ["-r", os.path.join(fcount_dir, file),
                          "-e", self.exp_design,
                          "-p", self.p_value,
                          "-n", name,
                          "-o", edger_dir]
                #TODO: get the output that has locus tag
                edger_cmd = EdgeR[edger_list]
                logger = logging.getLogger('luigi-interface')
                logger.info(edger_cmd)
                edger_cmd()
                if file == "gene_count.tsv":
                #TODO:convert the first column to locus tag
                    if self.pathway is True:
                        path_list = ["-d", edger_dir,
                                "-m", "edgeR", "-c", self.org_code] # get pathway information
                        path_cmd = plot_pathway[path_list]
                        logger.info(path_cmd)
                        path_cmd()
                    if self.GAGE is True:
                        gage_list = ["-d", edger_dir, "-m",
                                "edgeR", "-c", self.org_code]
                        gage_cmd = gage_analysis[gage_list]
                        logger.info(gage_cmd)
                        gage_cmd()
        self.summ_summ()

    def summ_summ(self):
        """Summarize the summary table to be displayed in edge"""
        edger_dir = os.path.join(self.workdir, "edgeR", self.kingdom)
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


@requires(Summ.FeatureCounts)
class DESeq2(luigi.Task):
    """Find DGE using DESeq2."""
    exp_design = luigi.Parameter()
    p_value = luigi.FloatParameter()
    prok_org_code = luigi.Parameter()
    euk_org_code = luigi.Parameter()
    GAGE = luigi.BoolParameter()
    pathway = luigi.BoolParameter()

    def output(self):
        """Expected output of DGE using DESeq2."""
        fcount_dir = os.path.join(self.workdir, "featureCounts", self.kingdom)
        DESeq2_dir = os.path.join(self.workdir, "DESeq2", self.kingdom)
        for file in os.listdir(fcount_dir):
            if file.endswith("__sig.csv"):
                out_folder = file.split(".csv")[0]
                out_filepath = os.path.join(DESeq2_dir, out_folder,
                                            "greater.csv")
                return LocalTarget(out_filepath)

    def run(self):
        """Run DESeq2."""
        fcount_dir = os.path.join(self.workdir, "featureCounts", self.kingdom)
        DESeq2_dir = os.path.join(self.workdir, "DESeq2", self.kingdom)
        if not os.path.exists(DESeq2_dir):
            os.makedirs(DESeq2_dir)
        for file in os.listdir(fcount_dir):
            if file.endswith("tsv"):
                feat_name = file.split("_")[-2]
                # only run deseq on gene and CDS
                if any(feat == feat_name for feat in ["CDS", "gene"]) is True:
                    deseq2_list = ["-r", os.path.join(fcount_dir, file),
                               "-e", self.exp_design, "-p", self.p_value,
                               "-n", feat_name,
                               "-o", DESeq2_dir]
                    deseq2_cmd = RDESeq2[deseq2_list]
                    logger = logging.getLogger('luigi-interface')
                    logger.info(deseq2_cmd)
                    deseq2_cmd()
            if file == "gene_count.tsv":
                if self.prok_org_code is None:
                    org_code=self.euk_org_code
                else:
                    org_code=self.prok_org_code
                if self.pathway is True:
                    path_list = ["-d", DESeq2_dir,
                         "-m", "DESeq2", "-c",
                         org_code] # get pathway information
                    path_cmd = plot_pathway[path_list]
                    logger.info(path_cmd)
                    path_cmd()
                if self.GAGE is True:
                    gage_list = ["-d", DESeq2_dir, "-m",
                         "DESeq2", "-c", self.org_code]
                    gage_cmd = gage_analysis[gage_list]
                    logger.info(gage_cmd)
                    gage_cmd()
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
    # gff_file=luigi.Parameter()

    exp_design = luigi.Parameter()
    p_value = luigi.FloatParameter()

    def output(self):
        """Expected output of DGE using edgeR."""
        bg_rdir = os.path.join(self.workdir, "bg_results", self.kingdom)
        return LocalTarget(bg_rdir)

    def run(self):
        """Run ballgown."""
        bg_dir = os.path.join(self.workdir, "bg_results", self.kingdom)
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
