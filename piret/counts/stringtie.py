#! /usr/bin/env python

"""Luigi Tasks to perform various RNA seq functions, from Mapping to Counting.

Mapping is done using hisat2 and counting is done using featurecounts
and stringtie
"""

from __future__ import print_function
import os
import luigi
import sys
dir_path = os.path.dirname(os.path.realpath(__file__))
lib_path = os.path.abspath(os.path.join(dir_path, '..'))
sys.path.append(lib_path)

from luigi import LocalTarget
from luigi import Parameter, IntParameter, DictParameter, ListParameter
from luigi.util import inherits, requires
import subprocess
from plumbum.cmd import gffread, hisat2
from plumbum.cmd import samtools, stringtie, mv, awk
from plumbum.cmd import STAR
import pandas as pd
from sys import stderr, exit
import logging
from collections import defaultdict as dd, Counter
from os.path import basename, splitext
from plumbum.cmd import stringtie, featureCounts
from piret.Runs import Map
from piret.miscs import RefFile
import luigi
import logging


class StringTieScores(luigi.Task):
    """Calculate scores using string tie."""
    gff_file = luigi.Parameter()
    out_gtf = luigi.Parameter()
    out_cover = luigi.Parameter()
    out_abun = luigi.Parameter()
    in_bam_file = luigi.Parameter()
    num_cpus = Parameter()

    def requires(self):
        """."""
        if ',' in self.gff_file:
            gffs = self.gff_file.split(",")
            return [RefFile(os.path.abspath(gff)) for gff in gffs]
        else:
            return [RefFile(os.path.abspath(self.gff_file))]

    def output(self):
        """Index output."""
        return LocalTarget(self.out_gtf)

    def run(self):
        """Run stringtie."""
        stringtie_opt = ["-o", self.out_gtf,
                         "-p", self.num_cpus,
                         "-G", self.gff_file,
                         "-C", self.out_cover,
                         "-A", self.out_abun,
                         self.in_bam_file]
        stringtie_cmd = stringtie[stringtie_opt]
        logger = logging.getLogger('luigi-interface')
        logger.info(stringtie_cmd)
        stringtie_cmd()

# @requires(StringTieScores)
class StringTieScoresW(luigi.WrapperTask):
    """Wrapper function for stringtie in all samples"""
    fastq_dic = luigi.DictParameter()
    gff_file = Parameter()
    kingdom = Parameter()
    workdir = Parameter()
    num_cpus = IntParameter()

    def requires(self):
        """A wrapper for running Stringtie scores on all samples."""
        for samp, fastq in self.fastq_dic.items():
            map_dir = os.path.join(self.workdir, "processes", "mapping", samp)
            trim_dir = os.path.join(self.workdir, "processes", "qc", samp)
            stng_dir = os.path.join(self.workdir, "processes", "stringtie", samp)
            if os.path.isdir(stng_dir) is False:
                os.makedirs(stng_dir)
            if self.kingdom in ['prokarya', 'eukarya']:
                if self.kingdom == 'prokarya':
                    apd = '_prok'
                elif self.kingdom == 'eukarya':
                    apd = '_euk'
                yield StringTieScores(num_cpus=self.num_cpus,
                                      gff_file=self.gff_file,
                                      out_gtf=os.path.join(stng_dir, samp +  apd + "_sTie.gtf"),
                                      out_cover=os.path.join(stng_dir, samp + apd + "_covered_sTie.gtf"),
                                      out_abun=os.path.join(stng_dir, samp +  apd + "_sTie.tab"),
                                      in_bam_file=os.path.join(map_dir, samp + "_srt.bam"))
            elif self.kingdom == 'both':
                prok_gff = self.gff_file.split(",")[0]
                euk_gff = self.gff_file.split(",")[1]
                yield StringTieScores(gff_file=self.gff_file.split(",")[0],
                                      num_cpus=self.num_cpus,
                                      out_gtf=stng_dir + "/" + samp + "_prok" + "_sTie.gtf",
                                      out_cover=stng_dir + "/" + samp + "_prok" + "_covered_sTie.gtf",
                                      out_abun=stng_dir + "/" + samp + "_prok" + "_sTie.tab",
                                      in_bam_file=map_dir + "/" + samp + "_srt_prok.bam")
                yield StringTieScores(gff_file=self.gff_file.split(",")[1],
                                      num_cpus=self.num_cpus,
                                      out_gtf=os.path.join(stng_dir, samp + "_euk_sTie.gtf"),
                                      out_cover=os.path.join(stng_dir, samp + "_euk_covered_sTie.gtf"),
                                      out_abun=os.path.join(stng_dir, samp + "_euk_sTie.tab"),
                                      in_bam_file=map_dir + "/" + samp + "_srt_euk.bam")

@requires(StringTieScoresW)
class MergeStringTies(luigi.Task):
    """Summarize mapped reads classification using StringTie."""

    def ouput(self):
        """Ouptut of string tie merge."""
        if self.kingdom in ['prokarya', 'eukarya']:
            return [RefFile(os.path.join(self.workdir, "processes",
                                             "stringtie", "stringtie_merged_transcript.gtf"))]
        elif self.kingdom == 'both':
            return [RefFile(self.workdir + "/" +
                                "prok_sTie_merged_transcript.gtf"),
                    RefFile(self.workdir + "/" +
                                "euk_sTie_merged_transcript.gtf")]

    def run(self):
        """Running stringtie merge."""
        samp_list=list(self.fastq_dic.keys())
        if self.kingdom in ['prokarya', 'eukarya']:
            if self.kingdom == 'prokarya':
                append_name = '_prok'
            elif self.kingdom == 'eukarya':
                append_name = '_euk'
            gtf = self.workdir + "/" + \
                self.gff_file.split("/")[-1].split(".gff")[0] + ".gtf"
            out_gtf_list = [os.path.join(self.workdir, "processes",
                            "stringtie", samp, samp + append_name +
                            "_sTie.gtf") for samp in samp_list]
            stie_cmd_opt = ["--merge", "-G", self.gff_file,
                            "-o", os.path.join(self.workdir, "processes",
                                               "stringtie",
                            "stringtie_merged_transcript.gtf")] + out_gtf_list
            stie_cmd = stringtie[stie_cmd_opt]
            stie_cmd()
        elif self.kingdom == 'both':
            prok_gff = self.gff_file.split(",")[0]
            euk_gff = self.gff_file.split(",")[1]

            prokout_gtf_list = [os.path.join(self.workdir, "processes",
                                             "stringtie", samp,
                                             samp + "_prok_sTie.gtf")
                                for samp in samp_list]
            eukout_gtf_list = [os.path.join(self.workdir, "processes",
                                            "stringtie", samp, 
                                            samp + "_euk_sTie.gtf")
                               for samp in samp_list]
            stie_cmd_euk_opt = ["--merge", "-G", euk_gff,
                                "-o", self.workdir +
                                "/" + "/euk_sTie_merged_transcript.gtf"] + eukout_gtf_list
            stie_cmd_euk = stringtie[stie_cmd_euk_opt]
            stie_cmd_euk()

            stie_cmd_prok_opt = ["--merge", "-G", prok_gff,
                                 "-o", self.workdir +
                                 "/" +
                                 "/prok_sTie_merged_transcript.gtf"] + prokout_gtf_list
            stie_cmd_prok = stringtie[stie_cmd_prok_opt]
            stie_cmd_prok()


class ReStringTieScores(luigi.Task):
    """Calculate scores using string tie."""

    gtf = luigi.Parameter()
    out_gtf = luigi.Parameter()
    out_cover = luigi.Parameter()
    out_abun = luigi.Parameter()
    in_bam_file = luigi.Parameter()
    num_cpus = luigi.Parameter()

    def requires(self):
        """Requires gff file and bam file."""
        return [RefFile(self.gtf), RefFile(self.in_bam_file)]

    def output(self):
        """Output gtf."""
        return LocalTarget(self.out_gtf)

    def run(self):
        """Run stringtie."""
        stie_opt = ["-o", self.out_gtf,
                    "-p", self.num_cpus,
                    "-G", self.gtf,
                    "-C", self.out_cover,
                    "-A", self.out_abun,
                    "-B", self.in_bam_file]
        logger = logging.getLogger('luigi-interface')
        stie_cmd = stringtie[stie_opt]
        logger.info(stie_cmd)
        stie_cmd()


class ReStringTieScoresW(luigi.WrapperTask):
    """Recount to only have common transcripts."""
    fastq_dic = luigi.DictParameter()
    num_cpus = IntParameter()
    workdir = luigi.Parameter()
    kingdom = luigi.Parameter()
    
    def requires(self):
        """A wrapper for restringtie processes."""
        for samp, fastq in self.fastq_dic.items():
            map_dir = os.path.join(self.workdir, "processes", "mapping", samp)
            if self.kingdom in ['prokarya', 'eukarya']:
                bg_dir = os.path.join(self.workdir, "processes", "ballgown", self.kingdom, samp)
                if os.path.isdir(bg_dir) is False:
                    os.makedirs(bg_dir)
                yield ReStringTieScores(num_cpus=self.num_cpus,
                                        gtf=os.path.join(self.workdir, "processes",
                                             "stringtie", "stringtie_merged_transcript.gtf"),
                                        out_gtf=os.path.join(bg_dir, samp + "_merged_sTie.gtf"),
                                        out_cover=os.path.join(bg_dir, samp + "_merged_covered_sTie.gtf"),
                                        out_abun=os.path.join(bg_dir, samp + "_merged_sTie.tab"),
                                        in_bam_file=os.path.join(map_dir, samp + "_srt.bam"))
            elif self.kingdom == 'both':
                euk_bg_dir = os.path.join(self.workdir, "processes", "ballgown", "eukarya", samp)
                prok_bg_dir = os.path.join(self.workdir, "processes", "ballgown", "prokarya", samp)
                if os.path.isdir(euk_bg_dir) is False:
                    os.makedirs(euk_bg_dir)
                    os.makedirs(prok_bg_dir)
                prok_gtf = os.path.join(self.workdir,
                                        "prok_sTie_merged_transcript.gtf")
                euk_gtf = os.path.join(self.workdir,
                                       "euk_sTie_merged_transcript.gtf")
                yield ReStringTieScores(num_cpus=self.num_cpus,
                                        gtf=prok_gtf,
                                        out_gtf=prok_bg_dir + "/" + samp + "_prok_sTie.gtf",
                                        out_cover=prok_bg_dir + "/" + samp + "_prok_covered_sTie.gtf",
                                        out_abun=prok_bg_dir + "/" + samp + "_prok_sTie.tab",
                                        in_bam_file=os.path.join(map_dir, samp + "_srt_prok.bam"))
                yield ReStringTieScores(num_cpus=self.num_cpus,
                                        gtf=euk_gtf, out_gtf=euk_bg_dir + "/" + samp + "_euk_sTie.gtf",
                                        out_cover=euk_bg_dir + "/" + samp + "_euk_covered_sTie.gtf",
                                        out_abun=euk_bg_dir + "/" + samp + "_euk_sTie.tab",
                                        in_bam_file=os.path.join(map_dir, samp + "_srt_euk.bam"))
