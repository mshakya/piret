#! /usr/bin/env python

"""Luigi Tasks to call opaver functions
"""

import os
import luigi
import sys
import numpy as np
dir_path = os.path.dirname(os.path.realpath(__file__))
lib_path = os.path.abspath(os.path.join(dir_path, '..'))
sys.path.append(lib_path)
from luigi.contrib.external_program import ExternalProgramTask
from luigi import LocalTarget
from plumbum.cmd import perl, time
import pandas as pd
from sys import stderr, exit
import logging
import glob
from piret.qc.FaQC import RefFile


class RunOpaver(luigi.Task):
    """Extract chromosomes of euk and prok fasta."""

    kingdom = luigi.Parameter()
    workdir = luigi.Parameter()
    method = luigi.Parameter()

    def requires(self):
        """de."""
        emapper_file = os.path.join(self.workdir, "processes", "emapper",
                                    self.kingdom,
                                    "emapper.emapper.annotations")
        return RefFile(emapper_file)

    def output(self):
        """Local target."""
        dge_dir = os.path.join(self.workdir, "processes", self.method,
                               self.kingdom, "CDS")
        dge_files = [f for f in glob.glob(dge_dir + "**/*sig.csv",
                     recursive=True)]
        for file in dge_files:
            cond_fdr = os.path.basename(file).split("__sig.csv")[0]
            opaver_outfile = os.path.join(self.workdir, "processes", "opaver",
                                                self.kingdom, self.method,
                                                cond_fdr, "opaver.tsv")
            return LocalTarget(opaver_outfile)

    def run(self):
        """."""
        if self.method == "edgeR":
            opaver_files = self.prep_opaver(method="edgeR")
            # self.run_opaver(method="edgeR", opavers=opaver_files)
        elif self.method == "DESeq2":
            opaver_files = self.prep_opaver(method="DESeq2")
            # self.run_opaver(method="DESeq2", opavers=opaver_files)

    def prep_opaver(self, method):
        out_files = []
        dge_dir = os.path.join(self.workdir, "processes", method,
                               self.kingdom, "CDS")
        dge_files = [f for f in glob.glob(dge_dir + "**/*sig.csv",
                     recursive=True)]
        emapper_file = os.path.join(self.workdir, "processes", "emapper",
                                    self.kingdom,
                                    "emapper.emapper.annotations")
        # read in the emapper file.
        emap = pd.read_csv(emapper_file, sep='\t', skiprows=[0, 1, 2],
                           skipinitialspace=True, skipfooter=3, header=None,
                           engine='python')
        emap1 = emap.reset_index()
        emap1.columns = emap1.iloc[0]
        emap2 = emap1.drop(0).set_index('#query_name')
        if os.path.exists(emapper_file) is False:
            sys.exit("emapper output is not found! No pathway analysis for you!")
        else:
            for file in dge_files:
                dge_df = pd.read_csv(file, sep=",", index_col=0)
                if method == "edgeR":
                    cond_fdr = os.path.basename(file).split("__et.csv")[0]
                    opaver_outdir = os.path.join(self.workdir, "processes",
                                              "opaver", self.kingdom,
                                              method, cond_fdr)
                    if os.path.exists(opaver_outdir) is False:
                        os.makedirs(opaver_outdir)
                    dge_df = dge_df.drop(["Chr", "Start", "End",
                                      "Strand", "Length", "logCPM"],
                                       axis=1)
                    # merge the two file.
                    opaver_df = pd.merge(dge_df, emap2, left_on="Geneid",
                                     right_index=True,
                                     how="inner").loc[:, ['Geneid', 'KEGG_ko',
                                                          'logFC',
                                                          'PValue']]
                    opaver_df.KEGG_ko = opaver_df.KEGG_ko.apply(self.remove_ko)
                    opaver_df.to_csv(os.path.join(opaver_outdir, "opaver.tsv"),
                                     sep="\t", index=False)
                    out_files.append(os.path.join(opaver_outdir, "opaver.tsv"))
                elif method == "DESeq2":
                    cond_fdr = os.path.basename(file).split("__sig.csv")[0]
                    opaver_outdir = os.path.join(self.workdir, "processes", "opaver",
                                            self.kingdom, method,
                                            cond_fdr)
                    if os.path.exists(opaver_outdir) is False:
                        os.makedirs(opaver_outdir)
                    dge_df = dge_df.drop(["baseMean", "lfcSE",
                                      "stat"], axis=1)
                    opaver_df = pd.merge(dge_df, emap2, left_index=True,
                                     right_index=True,
                                     how="inner").loc[:, ['KEGG_ko',
                                                          'log2FoldChange',
                                                          'pvalue', 'padj']]
                    opaver_df.KEGG_ko = opaver_df.KEGG_ko.apply(self.remove_ko)
                    opaver_df.to_csv(os.path.join(opaver_outdir, "opaver.tsv"),
                                 sep="\t")
                    out_files.append(os.path.join(opaver_outdir, "opaver.tsv"))
        return out_files

    def remove_ko(self, x):
        "removes a ko string"
        if x != "nan" and x not in ["01100", "01120", "01130"]:
            return str(x).replace("ko:", "")
        else:
            return None

    def run_opaver(self, method, opavers):
        """ Run opaver.pl script."""
        for opav_file in opavers:
            out_folder = os.path.dirname(opav_file)
            if os.path.exists(out_folder) is False:
                os.makedirs(out_folder)
            opaver_pl = ["opaver.pl", "--tran", opav_file,
                     "-o", out_folder]
            time[opaver_pl]()
            opaver_r = ["opaver.r", "--wd="+out_folder]
            time[opaver_r]()

