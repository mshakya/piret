#! /usr/bin/env python

"""Luigi Tasks to perform various RNA seq functions, from Mapping to Counting.

Mapping is done using hisat2 and counting is done using featurecounts
and stringtie
"""

from __future__ import print_function
import os
import luigi
import sys
import numpy as np
dir_path = os.path.dirname(os.path.realpath(__file__))
lib_path = os.path.abspath(os.path.join(dir_path, '..'))
sys.path.append(lib_path)
from luigi.contrib.external_program import ExternalProgramTask
from luigi import LocalTarget
opaver_dir = os.path.abspath(os.path.join(dir_path, "../../thirdparty/omics-pathway-viewer/scripts"))
from plumbum.cmd import perl
import pandas as pd
from sys import stderr, exit
import logging
import glob
# from Bio import SeqIO


class RunOpaver(luigi.Task):
    """Extract chromosomes of euk and prok fasta."""

    kingdom = luigi.Parameter()
    workdir = luigi.Parameter()
    method = luigi.Parameter()

    def requires(self):
        """de."""


    def output(self):
        """Local target."""
        if self.method == "edgeR":
            opaver_outfile = os.path.join(self.workdir, "processes", "opaver",
                                          self.kingdom, "edgeR", "opaver.tsv")
        elif self.method == "DESeq2":
            opaver_outfile = os.path.join(self.workdir, "processes", "opaver",
                                          self.kingdom, "DESeq2", "opaver.tsv")
        return LocalTarget(opaver_outfile)


    def run(self):
        """."""
        if self.method == "edgeR":
            self.prep_opaver(method="edgeR")
            self.run_opaver(method="edgeR")
        elif self.method == "DESeq2":
            self.prep_opaver(method="DESeq2")
            self.run_opaver(method="DESeq2")

    def prep_opaver(self, method):
        opaver_outdir = os.path.join(self.workdir, "processes", "opaver",
                                     self.kingdom, method)
        if os.path.exists(opaver_outdir) is False:
            os.makedirs(opaver_outdir)
        dge_dir = os.path.join(self.workdir, "processes", method,
                               self.kingdom, "CDS")
        dge_files = [f for f in glob.glob(dge_dir + "**/*et.csv",
                     recursive=True)]
        emapper_file = os.path.join(self.workdir, "processes", "emapper",
                                    "emapper.emapper.annotations")
        # read in the emapper file.
        emap = pd.read_csv(emapper_file, sep='\t', skiprows=[0, 1, 2],
                                   skipinitialspace=True, skipfooter=3,
                                   header=None, engine='python')
        emap1 = emap.reset_index()
        emap1.columns = emap1.iloc[0]
        emap2 = emap1.drop(0).set_index('#query_name')
        if os.path.exists(emapper_file) is False:
            sys.exit("emapper file is not present")
        else:
            for file in dge_files:
                dge_df = pd.read_csv(file, sep=",", index_col=0)
            if method == "edgeR":
                dge_df = dge_df.drop(["Chr", "Start", "End",
                                      "Strand", "Length", "logCPM"],
                                       axis=1)
                # merge the two file.
                opaver_df = pd.merge(dge_df, emap2, left_on="Geneid",
                                     right_on='#query_name',
                                     how="inner").loc[:, ['Geneid', 'KEGG_ko',
                                                          'logFC',
                                                          'PValue']]
                opaver_df.KEGG_ko = opaver_df.KEGG_ko.apply(self.remove_ko)
                opaver_df.to_csv(os.path.join(opaver_outdir, "opaver.tsv"),
                                 sep="\t", index=False)
            elif method == "DESeq2":
                dge_df = dge_df.drop(["baseMean", "lfcSE",
                                      "stat"], axis=1)
                opaver_df = pd.merge(dge_df, emap2, left_index=True,
                                     right_on='#query_name',
                                     how="inner").loc[:, ['KEGG_ko',
                                                          'log2FoldChange',
                                                          'pvalue', 'padj']]
                opaver_df.KEGG_ko = opaver_df.KEGG_ko.apply(self.remove_ko)
                opaver_df.to_csv(os.path.join(opaver_outdir, "opaver.tsv"),
                                 sep="\t")

    def remove_ko(self, x):
        "removes a ko string"
        if x != "nan":
            return str(x).replace("ko:", "")
        else:
            return None

    def run_opaver(self, method):
        """ Run opaver.pl script."""
        in_file = os.path.join(self.workdir, "processes", "opaver", self.kingdom,
                               method, "opaver.tsv")

        out_folder = os.path.join(self.workdir, "processes", "opaver", self.kingdom,
                                  method, "path_ids")
        if os.path.exists(out_folder) is False:
            os.makedirs(out_folder)
        opaver_pl = [os.path.join(opaver_dir, "opaver.pl"), "--tran", in_file,
                     "-o", out_folder]
        perl[opaver_pl]()

