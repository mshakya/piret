
#! /usr/bin/env python

"""Check design."""
from __future__ import print_function
import os
from os.path import basename, splitext
from plumbum.cmd import stringtie, featureCounts
# from piret.runs import Map
import pandas as pd
from luigi.util import inherits, requires
from luigi.contrib.external_program import ExternalProgramTask
from luigi import LocalTarget, Parameter, IntParameter
import luigi
import logging


class FeatureCounts(luigi.Task):
    """Summarize mapped reads classificaion using FeatureCount."""

    fastq_dic = luigi.DictParameter()
    kingdom = luigi.Parameter()
    gff_file = luigi.Parameter()
    workdir = luigi.Parameter()
    indexfile = luigi.Parameter()
    num_cpus = luigi.IntParameter()
    ref_file = luigi.Parameter()
    fid = luigi.Parameter()
    stranded = luigi.IntParameter()
    feat2count = luigi.ListParameter()

    def output(self):
        """Expected output of featureCounts."""
        counts_dir = os.path.join(self.workdir, "processes", "featureCounts",
                                  self.kingdom)
        gff_fp = os.path.abspath(self.gff_file)
        # features = list(set(pd.read_csv(gff_fp, sep="\t", header=None,
        #                                 comment='#')[2].tolist()))
        features = [feat for feat in self.feat2count]
        loc_target = LocalTarget(os.path.join(counts_dir, features[-1] +
                                              "_count.tsv"))
        return loc_target

    def run(self):
        """Running featureCounts on all."""
        map_dir = os.path.join(self.workdir, "processes", "mapping")
        samp_list = list(self.fastq_dic.keys())
        in_srtbam_list = [os.path.join(map_dir, samp, samp + "_srt.bam")
                          for samp in samp_list]
        counts_dir = os.path.join(self.workdir, "processes", "featureCounts",
                                  self.kingdom)
        if not os.path.exists(counts_dir):
            os.makedirs(counts_dir)
        if ',' in self.gff_file:
            gff_list = self.gff_file.split(",")
            gff_full_path = [os.path.abspath(gff) for gff in gff_list]
            for gffs in gff_full_path:
                feature = list(set(pd.read_csv(gffs,
                                               sep="\t", header=None,
                                               comment='#')[2].tolist()))
                for feat in feature:
                    if feat in self.feat2count:
                        fcount_cmd_opt = ["-a", self.gff_file,
                                          "-s", self.stranded,
                                          "-B",
                                          "-p", "-P", "-C",
                                          "-g", self.fid,
                                          "-t", feat,
                                          "-T", self.num_cpus,
                                          "-o", counts_dir + "/" + feat +
                                          "_count.tsv"] + in_srtbam_list
                        fcount_cmd = featureCounts[fcount_cmd_opt]
                        fcount_cmd()
                    elif feat in ['gene']:
                        fcount_cmd_opt = ["-a", self.gff_file,
                                          "-s", self.stranded,
                                          "-B",
                                          "-p", "-P", "-C",
                                          "-g", "ID",
                                          "-t", feat,
                                          "-T", self.num_cpus,
                                          "-o", counts_dir + "/" + feat +
                                                "_count.tsv"] + in_srtbam_list
                        fcount_cmd = featureCounts[fcount_cmd_opt]
                        fcount_cmd()
                    else:
                        pass
        else:
            feature = list(set(pd.read_csv(self.gff_file, sep="\t", header=None,
                                           comment='#')[2].tolist()))
        for feat in feature:
            if feat in self.feat2count:
                fcount_cmd_opt = ["-a", self.gff_file,
                                  "-s", self.stranded,
                                  "-B",
                                  "-p", "-P", "-C",
                                  "-g", self.fid,
                                  "-t", feat,
                                  "-T", self.num_cpus,
                                  "-o", counts_dir + "/" + feat +
                                        "_count.tsv"] + in_srtbam_list
                fcount_cmd = featureCounts[fcount_cmd_opt]
                fcount_cmd()
            if feat in ['gene']:
                fcount_cmd_opt = ["-a", self.gff_file,
                                  "-s", self.stranded,
                                  "-B",
                                  "-p", "-P", "-C",
                                  "-g", "ID",
                                  "-t", feat,
                                  "-T", self.num_cpus,
                                  "-o", counts_dir + "/" + feat +
                                  "_count.tsv"] + in_srtbam_list
                fcount_cmd = featureCounts[fcount_cmd_opt]
                fcount_cmd()
