#! /usr/bin/env python

"""Check design."""
from __future__ import print_function
import os
import luigi
from plumbum.cmd import stringtie
from luigi.contrib.external_program import ExternalProgramTask
from luigi import LocalTarget
from luigi.util import inherits
from pypiret import Map


class RefFile(luigi.ExternalTask):
    """An ExternalTask like this."""

    path = luigi.Parameter()

    def output(self):
        """Check."""
        return LocalTarget(os.path.abspath(self.path))


class FeatureCounts(ExternalProgramTask):
    """Summarize mapped reads classificaion using FeatureCount."""

    kingdom = luigi.Parameter()
    numCPUs = luigi.IntParameter()
    prok_gff = luigi.Parameter()
    euk_gff = luigi.Parameter()
    fastq_dic = luigi.DictParameter()
    workdir = luigi.Parameter()
    bindir = luigi.Parameter()
    indexfile = luigi.Parameter()

    def requires(self):
        """Require reference fasta format file."""
        samp_list = list(self.fastq_dic.keys())
        if self.kingdom == 'prokarya':
            bam_filelist = [self.workdir + "/" + samp + "/" +
                            "mapping_results" + "/" +
                            "prokarya.bam" for samp in samp_list]
            return[RefFile(bam) for bam in bam_filelist]
        elif self.kingdom == 'eukarya':
            bam_filelist = [self.workdir + "/" + samp + "/" +
                            "mapping_results" + "/" +
                            "eukarya.bam" for samp in samp_list]
            return[RefFile(bam) for bam in bam_filelist]

    def output(self):
        """Index output."""
        if self.kingdom == 'prokarya':
            return LocalTarget(self.workdir + "/prok.count")
        elif self.kingdom == 'eukarya':
            return LocalTarget(self.workdir + "/euk.count")

    def program_args(self):
        """Run featureCounts."""
        samp_list = list(self.fastq_dic.keys())
        if self.kingdom == 'prokarya':
            prok_gtf = self.workdir + "/" + \
                self.prok_gff.split("/")[-1].split(".gff")[0] + ".gtf"
            bam_filelist = [self.workdir + "/" + samp + "/" +
                            "mapping_results" + "/" +
                            "prokarya.bam" for samp in samp_list]
            return ["featureCounts",
                    "-a", prok_gtf,
                    "-s", 1,
                    "-B",
                    "-p", "-P",
                    "-o", self.workdir + "/prok.count"] + bam_filelist
        elif self.kingdom == 'eukarya':
            euk_gtf = self.workdir + "/" + \
                self.euk_gff.split("/")[-1].split(".gff")[0] + ".gtf"
            bam_filelist = [self.workdir + "/" + samp + "/" +
                            "mapping_results" + "/" +
                            "eukarya.bam" for samp in samp_list]
            return ["featureCounts",
                    "-a", euk_gtf,
                    "-s", 1,
                    "-p", "-P", "-B",
                    "-o", self.workdir + "/euk.count"] + bam_filelist

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
        if self.kingdom in ['prokarya', 'eukarya']:
            gtf = self.workdir + "/" + \
                self.gff.split("/")[-1].split(".gff")[0] + ".gtf"
            samp_list = list(self.fastq_dic.keys())
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
                self.prok_gff.split("/")[-1].split(".gff")[0] + ".gtf"

            euk_gtf = self.workdir + "/" + \
                self.euk_gff.split("/")[-1].split(".gff")[0] + ".gtf"

            stie_cmd_euk = stringtie["--merge", "-G", euk_gtf,
                                     "-o", self.workdir +
                                     "/" + "/euk_merged_transcript.gtf"] + self.eukout_gtf_list
            stie_cmd_euk()

            stie_cmd_prok = stringtie["--merge", "-G", prok_gtf,
                                      "-o", self.workdir +
                                      "/" + "/prok_merged_transcript.gtf"] + self.prokout_gtf_list
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

    gff = luigi.Parameter()
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
            bg_dir = self.workdir + "/ballgown/" + samp
            if os.path.isdir(bg_dir) is False:
                os.makedirs(bg_dir)
            if self.kingdom in ['prokarya', 'eukarya']:
                gtf = self.workdir + "/" + \
                    self.gff.split("/")[-1].split(".gff")[0] + ".gtf"
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
                                        gtf=gtf,
                                        out_gtf=bg_dir + "/" + samp + "_sTie.gtf",
                                        out_cover=bg_dir + "/" + samp + "_covered_sTie.gtf",
                                        out_abun=bg_dir + "/" + samp + "_sTie.tab",
                                        in_bam_file=map_dir + "/" + samp + "_srt.bam",
                                        bindir=self.bindir)
            elif self.kingdom == 'both':
                euk_gtf = self.workdir + "/" + \
                    self.euk_gff.split("/")[-1].split(".gff")[0] + ".gtf"
                prok_gtf = self.workdir + "/" + \
                    self.prok_gff.split("/")[-1].split(".gff")[0] + ".gtf"
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
                                        gtf=prok_gtf,
                                        out_gtf=map_dir + "/" + samp + "_prok_sTie.gtf",
                                        out_cover=map_dir + "/" + samp + "_prok_covered_sTie.gtf",
                                        out_abun=map_dir + "/" + samp + "_prok_sTie.tab",
                                        in_bam_file=map_dir + "/prokarya.bam",
                                        bindir=self.bindir)
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
                                        gtf=euk_gtf,
                                        out_gtf=map_dir + "/" + samp + "_euk_sTie.gtf",
                                        out_cover=map_dir + "/" + samp + "_euk_covered_sTie.gtf",
                                        out_abun=map_dir + "/" + samp + "_euk_sTie.tab",
                                        in_bam_file=map_dir + "/eukarya.bam",
                                        bindir=self.bindir)
