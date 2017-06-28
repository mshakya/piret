#! /usr/bin/env python

"""Check design."""
from __future__ import print_function
import os
import luigi
from plumbum.cmd import stringtie, featureCounts
from luigi.contrib.external_program import ExternalProgramTask
from luigi import LocalTarget
from pypiret import Map
import pandas as pd
from luigi.util import inherits, requires


class RefFile(luigi.ExternalTask):
    """An ExternalTask like this."""

    path = luigi.Parameter()

    def output(self):
        """Check."""
        return LocalTarget(os.path.abspath(self.path))


@requires(Map.SortBAMfileW)
class FeatureCounts(luigi.Task):
    """Summarize mapped reads classificaion using FeatureCount."""

    kingdom = luigi.Parameter()
    euk_gff = luigi.Parameter()
    prok_gff = luigi.Parameter()
    workdir = luigi.Parameter()
    indexfile = luigi.Parameter()
    bindir = luigi.Parameter()
    numCPUs = luigi.Parameter()
    scriptdir = luigi.Parameter()
    ref_file = luigi.Parameter()

    def output(self):
        """Expected output of featureCounts."""
        if self.kingdom == 'prokarya':
            prok_gtf = self.workdir + "/" + self.prok_gff.split("/")[-1].split(".gff")[0] + ".gtf"
            features = list(set(pd.read_csv(prok_gtf, sep="\t", header=None)[2].tolist()))
            loc_target = [LocalTarget(self.workdir + feat + ".count") for feat in features]
            return loc_target
        elif self.kingdom == 'eukarya':
            euk_gtf = self.workdir + "/" + self.euk_gff.split("/")[-1].split(".gff")[0] + ".gtf"
            features = list(set(pd.read_csv(euk_gtf, sep="\t", header=None)[2].tolist()))
            loc_target = [LocalTarget(self.workdir + feat + ".count") for feat in features]
            return loc_target

    def run(self):
        """Running featureCounts on all."""
        samp_list = list(self.fastq_dic.keys())
        in_srtbam_list = [self.workdir + "/" + samp + "/" +
                          "mapping_results" + "/" + samp + "_srt.bam"
                          for samp in samp_list]
        if self.kingdom == 'prokarya':
            prok_gtf = self.workdir + "/" + self.prok_gff.split("/")[-1].split(".gff")[0] + ".gtf"
            features = list(set(pd.read_csv(prok_gtf, sep="\t", header=None)[2].tolist()))
            for feat in features:
                fcount_cmd_opt = ["-a", prok_gtf,
                                  "-s", 1,
                                  "-B",
                                  "-p", "-P", "-C",
                                  "-g", "transcript_id",
                                  "-t", feat,
                                  "-T", self.numCPUs,
                                  "-o", self.workdir + "/" + feat + ".count"] + in_srtbam_list
                fcount_cmd = featureCounts[fcount_cmd_opt]
                fcount_cmd()
        if self.kingdom == 'eukarya':
            euk_gtf = self.workdir + "/" + self.euk_gff.split("/")[-1].split(".gff")[0] + ".gtf"
            features = list(set(pd.read_csv(euk_gtf, sep="\t", header=None)[2].tolist()))
            for feat in features:
                fcount_cmd_opt = ["-a", euk_gtf,
                                  "-s", 1,
                                  "-B",
                                  "-p", "-P", "-C",
                                  "-g", "transcript_id",
                                  "-t", feat,
                                  "-T", self.numCPUs,
                                  "-o", self.workdir + "/" + feat + ".count"] + in_srtbam_list
                fcount_cmd = featureCounts[fcount_cmd_opt]
                fcount_cmd()

    def program_environment(self):
        """Environmental variables for this program."""
        return {'PATH': os.environ["PATH"] + ":" + self.bindir}


@requires(Map.SortBAMfileW)
class FeatureCountsBoth(luigi.Task):
    """Summarize mapped reads classificaion using FeatureCount."""

    kingdom = luigi.Parameter()
    euk_gff = luigi.Parameter()
    prok_gff = luigi.Parameter()
    workdir = luigi.Parameter()
    indexfile = luigi.Parameter()
    bindir = luigi.Parameter()
    numCPUs = luigi.Parameter()
    scriptdir = luigi.Parameter()
    ref_file = luigi.Parameter()

    def output(self):
        """Index output."""
        prok_gtf = self.workdir + "/" + \
            self.prok_gff.split("/")[-1].split(".gff")[0] + ".gtf"
        prok_features = list(set(pd.read_csv(prok_gtf, sep="\t", header=None)[2].tolist()))
        prok_target = [LocalTarget(self.workdir + "/prok_" + feat + ".count") for feat in prok_features]
        euk_gtf = self.workdir + "/" + \
            self.euk_gff.split("/")[-1].split(".gff")[0] + ".gtf"
        euk_features = list(set(pd.read_csv(euk_gtf, sep="\t", header=None)[2].tolist()))
        euk_target = [LocalTarget(self.workdir + "/euk_" + feat + ".count") for feat in euk_features]
        loc_target = prok_target + euk_target
        return loc_target

    def run(self):
        """Running featureCounts on all."""
        samp_list = list(self.fastq_dic.keys())
        in_srtbam_list = [self.workdir + "/" + samp + "/" +
                          "mapping_results" + "/" + samp + "_srt.bam"
                          for samp in samp_list]
        prok_gtf = self.workdir + "/" + \
            self.prok_gff.split(";")[0].split(
                "/")[-1].split(".gff")[0] + ".gtf"
        prok_features = list(set(pd.read_csv(prok_gtf, sep="\t", header=None)[2].tolist()))
        euk_gtf = self.workdir + "/" + \
            self.euk_gff.split(";")[0].split(
                "/")[-1].split(".gff")[0] + ".gtf"
        euk_features = list(set(pd.read_csv(euk_gtf, sep="\t", header=None)[2].tolist()))
        for feat in euk_features:
            fcount_euk_cmd_opt = ["-a", euk_gtf,
                                  "-s", 1,
                                  "-B",
                                  "-p", "-P", "-C",
                                  "-g", "transcript_id",
                                  "-t", feat,
                                  "-T", self.numCPUs,
                                  "-o", self.workdir + "/euk_" + feat + ".count"] + in_srtbam_list
            fcount_euk_cmd = featureCounts[fcount_euk_cmd_opt]
            fcount_euk_cmd()
        for feat in prok_features:
            fcount_prok_cmd_opt = ["-a", prok_gtf,
                                   "-s", 1,
                                   "-B",
                                   "-p", "-P", "-C",
                                   "-g", "transcript_id",
                                   "-t", feat,
                                   "-T", self.numCPUs,
                                   "-o", self.workdir + "/prok_" + feat + ".count"] + in_srtbam_list
            fcount_prok_cmd = featureCounts[fcount_prok_cmd_opt]
            fcount_prok_cmd()

    def program_environment(self):
        """Environmental variables for this program."""
        return {'PATH': os.environ["PATH"] + ":" + self.bindir}


@inherits(Map.StringTieScoresW)
class MergeStringTies(luigi.Task):
    """Summarize mapped reads classification using StringTie."""

    def requires(self):
        """Required."""
        samp_list = list(self.fastq_dic.keys())
        if self.kingdom in ['prokarya', 'eukarya']:
            out_gtf_list = [self.workdir + "/" + samp + "/" +
                            "mapping_results" + "/" + samp + "_sTie.gtf"
                            for samp in samp_list]
            return[RefFile(out_gtf_file) for out_gtf_file in out_gtf_list]
        elif self.kingdom == 'both':
            prokout_gtf_list = [self.workdir + "/" + samp + "/" +
                                "mapping_results" + "/" + samp + "_prok_sTie.gtf"
                                for samp in samp_list]
            eukout_gtf_list = [self.workdir + "/" + samp + "/" +
                               "mapping_results" + "/" + samp + "_euk_sTie.gtf"
                               for samp in samp_list]
            out_gtf_list = prokout_gtf_list + eukout_gtf_list
            return[Map.RefFile(out_gtf) for out_gtf in out_gtf_list]

    def ouput(self):
        """Ouptut of string tie merge."""
        if self.kingdom in ['prokarya', 'eukarya']:
            return [Map.RefFile(self.workdir + "/" + "merged_transcript.gtf")]
        elif self.kingdom == 'both':
            return [Map.RefFile(self.workdir + "/" + "prok_merged_transcript.gtf"),
                    Map.RefFile(self.workdir + "/" + "euk_merged_transcript.gtf")]

    def run(self):
        """Running stringtie merge."""
        samp_list = list(self.fastq_dic.keys())
        if self.kingdom in ['prokarya', 'eukarya']:
            gtf = self.workdir + "/" + \
                self.gff_file.split("/")[-1].split(".gff")[0] + ".gtf"
            out_gtf_list = [self.workdir + "/" + samp + "/" +
                            "mapping_results" + "/" + samp + "_sTie.gtf"
                            for samp in samp_list]

            stie_cmd_opt = ["--merge", "-G", gtf,
                            "-o", self.workdir +
                            "/" + "/merged_transcript.gtf"] + out_gtf_list
            stie_cmd = stringtie[stie_cmd_opt]
            stie_cmd()
        elif self.kingdom == 'both':
            prok_gtf = self.workdir + "/" + \
                self.gff_file.split(";")[0].split(
                    "/")[-1].split(".gff")[0] + ".gtf"
            euk_gtf = self.workdir + "/" + \
                self.gff_file.split(";")[1].split(
                    "/")[-1].split(".gff")[0] + ".gtf"

            prokout_gtf_list = [self.workdir + "/" + samp + "/" +
                                "mapping_results" + "/" + samp + "_prok_sTie.gtf"
                                for samp in samp_list]
            eukout_gtf_list = [self.workdir + "/" + samp + "/" +
                               "mapping_results" + "/" + samp + "_euk_sTie.gtf"
                               for samp in samp_list]
            stie_cmd_euk_opt = ["--merge", "-G", euk_gtf,
                                "-o", self.workdir +
                                "/" + "/euk_merged_transcript.gtf"] + eukout_gtf_list
            stie_cmd_euk = stringtie[stie_cmd_euk_opt]
            stie_cmd_euk()

            stie_cmd_prok_opt = ["--merge", "-G", prok_gtf,
                                 "-o", self.workdir +
                                 "/" + "/prok_merged_transcript.gtf"] + prokout_gtf_list
            stie_cmd_prok = stringtie[stie_cmd_prok_opt]
            stie_cmd_prok()


@inherits(Map.StringTieScores)
class ReStringTieScores(ExternalProgramTask):
    """Calculate scores using string tie."""

    gtf = luigi.Parameter()
    out_gtf = luigi.Parameter()
    out_cover = luigi.Parameter()
    out_abun = luigi.Parameter()
    in_bam_file = luigi.Parameter()

    def requires(self):
        """Require reference fasta format file."""
        return [RefFile(self.gtf), RefFile(self.in_bam_file)]

    def output(self):
        """Index output."""
        return LocalTarget(self.out_gtf)

    def program_args(self):
        """Run stringtie."""
        return ["stringtie",
                "-o", self.out_gtf,
                "-G", self.gtf,
                "-C", self.out_cover,
                "-A", self.out_abun,
                "-B",
                self.in_bam_file]


@inherits(Map.SortBAMfileW)
class ReStringTieScoresW(luigi.WrapperTask):
    """Recount to only have common transcripts."""

    gff_file = luigi.Parameter()
    kingdom = luigi.Parameter()
    workdir = luigi.Parameter()
    fastq_dic = luigi.DictParameter()

    def requires(self):
        """A wrapper for running the QC."""
        splice_list = [self.workdir + "/" +
                       f for f in os.listdir(self.workdir) if f.endswith('.splice')]
        if len(splice_list) > 1:
            splice_file = ','.join([splice_list])
        elif len(splice_list) == 1:
            splice_file = splice_list[0]
        else:
            splice_file = ''

        for samp, fastq in self.fastq_dic.iteritems():
            map_dir = self.workdir + "/" + samp + "/mapping_results"
            trim_dir = self.workdir + "/" + samp + "/trimming_results"
            if self.kingdom in ['prokarya', 'eukarya']:
                bg_dir = self.workdir + "/ballgown/" + samp
                if os.path.isdir(bg_dir) is False:
                    os.makedirs(bg_dir)
                gtf = self.workdir + "/" + \
                    self.gff_file.split("/")[-1].split(".gff")[0] + ".gtf"
                yield ReStringTieScores(fastq1=trim_dir + "/" + samp + ".1.trimmed.fastq",
                                        fastq2=trim_dir + "/" + samp + ".2.trimmed.fastq",
                                        numCPUs=self.numCPUs,
                                        indexfile=self.indexfile,
                                        spliceFile=splice_file,
                                        mappingLogFile=map_dir + "/mapping.log",
                                        unalned=map_dir + "/unligned.fastq",
                                        outsam=map_dir + "/" + samp + ".sam",
                                        bam_file=map_dir + "/" + samp + ".bam",
                                        sorted_bam_file=map_dir + "/" + samp + "_srt.bam",
                                        scriptdir=self.scriptdir,
                                        ref_file=self.ref_file,
                                        gff_file=self.gff_file,
                                        gtf=gtf,
                                        out_gtf=bg_dir + "/" + samp + "_sTie.gtf",
                                        out_cover=bg_dir + "/" + samp + "_covered_sTie.gtf",
                                        out_abun=bg_dir + "/" + samp + "_sTie.tab",
                                        in_bam_file=map_dir + "/" + samp + "_srt.bam",
                                        bindir=self.bindir,
                                        workdir=self.workdir)
            elif self.kingdom == 'both':
                bg_dir_prok = self.workdir + "/ballgown/prok/" + samp
                if os.path.isdir(bg_dir_prok) is False:
                    os.makedirs(bg_dir_prok)
                bg_dir_euk = self.workdir + "/ballgown/euk/" + samp
                if os.path.isdir(bg_dir_euk) is False:
                    os.makedirs(bg_dir_euk)
                prok_gtf = self.workdir + "/" + \
                    self.gff_file.split(";")[0].split(
                        "/")[-1].split(".gff")[0] + ".gtf"
                euk_gtf = self.workdir + "/" + \
                    self.gff_file.split(";")[1].split(
                        "/")[-1].split(".gff")[0] + ".gtf"

                yield ReStringTieScores(fastq1=trim_dir + "/" + samp + ".1.trimmed.fastq",
                                        fastq2=trim_dir + "/" + samp + ".2.trimmed.fastq",
                                        numCPUs=self.numCPUs,
                                        indexfile=self.indexfile,
                                        spliceFile=splice_file,
                                        mappingLogFile=map_dir + "/mapping.log",
                                        unalned=map_dir + "/unligned.fastq",
                                        outsam=map_dir + "/" + samp + ".sam",
                                        bam_file=map_dir + "/" + samp + ".bam",
                                        sorted_bam_file=map_dir + "/" + samp + "_srt.bam",
                                        scriptdir=self.scriptdir,
                                        ref_file=self.ref_file,
                                        gff_file=self.gff_file,
                                        gtf=prok_gtf,
                                        out_gtf=bg_dir_prok + "/" + samp + "_prok_sTie.gtf",
                                        out_cover=bg_dir_prok + "/" + samp + "_prok_covered_sTie.gtf",
                                        out_abun=bg_dir_prok + "/" + samp + "_prok_sTie.tab",
                                        in_bam_file=map_dir + "/prokarya.bam",
                                        bindir=self.bindir,
                                        workdir=self.workdir)
                yield ReStringTieScores(fastq1=trim_dir + "/" + samp + ".1.trimmed.fastq",
                                        fastq2=trim_dir + "/" + samp + ".2.trimmed.fastq",
                                        numCPUs=self.numCPUs,
                                        indexfile=self.indexfile,
                                        spliceFile=splice_file,
                                        mappingLogFile=map_dir + "/mapping.log",
                                        unalned=map_dir + "/unligned.fastq",
                                        outsam=map_dir + "/" + samp + ".sam",
                                        bam_file=map_dir + "/" + samp + ".bam",
                                        sorted_bam_file=map_dir + "/" + samp + "_srt.bam",
                                        scriptdir=self.scriptdir,
                                        ref_file=self.ref_file,
                                        gff_file=self.gff_file,
                                        gtf=euk_gtf,
                                        out_gtf=bg_dir_euk + "/" + samp + "_euk_sTie.gtf",
                                        out_cover=bg_dir_euk + "/" + samp + "_euk_covered_sTie.gtf",
                                        out_abun=bg_dir_euk + "/" + samp + "_euk_sTie.tab",
                                        in_bam_file=map_dir + "/eukarya.bam",
                                        bindir=self.bindir,
                                        workdir=self.workdir)
