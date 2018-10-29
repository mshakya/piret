#! /usr/bin/env python

"""Check design."""
from __future__ import print_function
import os
from os.path import basename, splitext
from plumbum.cmd import stringtie, featureCounts
from pypiret import Map
import pandas as pd
from luigi.util import inherits, requires
from luigi.contrib.external_program import ExternalProgramTask
from luigi import LocalTarget
import luigi

class RefFile(luigi.ExternalTask):
    """An ExternalTask like this."""

    path = luigi.Parameter()

    def output(self):
        """Check."""
        return LocalTarget(os.path.abspath(self.path))


@requires(Map.HisatMapW)
class FeatureCounts(luigi.Task):
    """Summarize mapped reads classificaion using FeatureCount."""

    kingdom = luigi.Parameter()
    gff = luigi.Parameter()
    workdir = luigi.Parameter()
    indexfile = luigi.Parameter()
    num_cpus = luigi.IntParameter()
    ref_file = luigi.Parameter()
    fid = luigi.Parameter()
    stranded = luigi.IntParameter()

    def output(self):
        """Expected output of featureCounts."""
        counts_dir = os.path.join(self.workdir, "featureCounts",
                                  self.kingdom)
        if ',' in self.gff:
            gff_list = self.gff.split(",")
            gff_full_path = [os.path.abspath(gff) for gff in gff_list]
            all_target = []
            for gffs in gff_full_path:
                feature = list(set(pd.read_csv(gffs, sep="\t", header=None, comment='#')[2].tolist()))
                loc_target = [LocalTarget(counts_dir + feat + "_count.tsv") for feat in feature]
                all_target = loc_target + all_target
                return all_target
        else:
            gff_fp = os.path.abspath(self.gff)
            features = list(set(pd.read_csv(gff_fp, sep="\t", header=None, comment='#')[2].tolist()))
            features = [feat for feat in features if feat in ['CDS', 'rRNA', 'tRNA', 'exon', 'gene', 'transcript']]
            loc_target = LocalTarget(counts_dir + "/" + features[-1] + "_count.tsv")
            return loc_target

    def run(self):
        """Running featureCounts on all."""
        samp_list = list(self.fastq_dic.keys())
        in_srtbam_list = [self.workdir + "/" + samp + "/" +
                          "mapping_results" + "/" + samp + "_srt.bam"
                          for samp in samp_list]
        counts_dir = os.path.join(self.workdir, "featureCounts",
                                  self.kingdom)
        if not os.path.exists(counts_dir):
            os.makedirs(counts_dir)
        if ',' in self.gff:
            gff_list = self.gff.split(",")
            gff_full_path = [os.path.abspath(gff) for gff in gff_list]
            for gffs in gff_full_path:
                feature = list(set(pd.read_csv(gffs,
                                               sep="\t", header=None,
                                               comment='#')[2].tolist()))
                for feat in feature:
                    if feat in ['CDS', 'rRNA', 'tRNA', 'exon',
                                'novel_region', 'transcript', 'mRNA']:
                        fcount_cmd_opt = ["-a", self.gff,
                                          "-s", self.stranded,
                                          "-B",
                                          "-p", "-P", "-C",
                                          "-g", self.fid,
                                          "-t", feat,
                                          "-T", self.num_cpus,
                                          "-o", counts_dir + "/" + feat + "_count.tsv"] + in_srtbam_list
                    elif feat in ['gene']:
                        fcount_cmd_opt = ["-a", self.gff,
                                      "-s", self.stranded,
                                      "-B",
                                      "-p", "-P", "-C",
                                      "-g", "locus_tag",
                                      "-t", feat,
                                      "-T", self.num_cpus,
                                      "-o", counts_dir + "/" + feat + "_count.tsv"] + in_srtbam_list
                    fcount_cmd = featureCounts[fcount_cmd_opt]
                    fcount_cmd()
        else:
            feature = list(set(pd.read_csv(self.gff, sep="\t", header=None,
                                        comment='#')[2].tolist()))
        for feat in feature:
            if feat in ['CDS', 'rRNA', 'tRNA', 'exon', 'transcript',
                        'novel_region']:
                fcount_cmd_opt = ["-a", self.gff,
                                    "-s", self.stranded,
                                    "-B",
                                    "-p", "-P", "-C",
                                    "-g", self.fid,
                                    "-t", feat,
                                    "-T", self.num_cpus,
                                    "-o", counts_dir + "/" + feat + "_count.tsv"] + in_srtbam_list
                fcount_cmd = featureCounts[fcount_cmd_opt]
                fcount_cmd()
            if feat in ['gene']:
                fcount_cmd_opt = ["-a", self.gff,
                                    "-s", self.stranded,
                                    "-B",
                                    "-p", "-P", "-C",
                                    "-g", "locus_tag",
                                    "-t", feat,
                                    "-T", self.num_cpus,
                                    "-o", counts_dir + "/" + feat + "_count.tsv"] + in_srtbam_list
                fcount_cmd = featureCounts[fcount_cmd_opt]
                fcount_cmd()


@requires(Map.StringTieScoresW)
class MergeStringTies(luigi.Task):
    """Summarize mapped reads classification using StringTie."""

    def ouput(self):
        """Ouptut of string tie merge."""
        if self.kingdom in ['prokarya', 'eukarya']:
            return [Map.RefFile(self.workdir + "/" + "stringtie_merged_transcript.gtf")]
        elif self.kingdom == 'both':
            return [Map.RefFile(self.workdir + "/" + "prok_sTie_merged_transcript.gtf"),
                    Map.RefFile(self.workdir + "/" + "euk_sTie_merged_transcript.gtf")]

    def run(self):
        """Running stringtie merge."""
        samp_list = list(self.fastq_dic.keys())
        if self.kingdom in ['prokarya', 'eukarya']:
            if self.kingdom == 'prokarya':
                append_name = '_prok'
            elif self.kingdom == 'eukarya':
                append_name = '_euk'
            gtf = self.workdir + "/" + \
                self.gff_file.split("/")[-1].split(".gff")[0] + ".gtf"
            out_gtf_list = [self.workdir + "/" + samp + "/" + "stie_results" +
                            "/" + samp + append_name +
                            "_sTie.gtf" for samp in samp_list]
            stie_cmd_opt = ["--merge", "-G", self.gff_file,
                            "-o", self.workdir +
                            "/" + "stringtie_merged_transcript.gtf"] + out_gtf_list
            stie_cmd = stringtie[stie_cmd_opt]
            stie_cmd()
        elif self.kingdom == 'both':
            prok_gff = self.gff_file.split(",")[0]
            euk_gff = self.gff_file.split(",")[1]

            prokout_gtf_list = [os.path.join(self.workdir, samp, 
                                             "stie_results",
                                             samp + "_prok_sTie.gtf")
                                for samp in samp_list]
            eukout_gtf_list = [os.path.join(self.workdir, samp, 
                                             "stie_results",
                                             samp + "_euk_sTie.gtf")
                               for samp in samp_list]
            stie_cmd_euk_opt = ["--merge", "-G", euk_gff,
                                "-o", self.workdir +
                                "/" + "/euk_sTie_merged_transcript.gtf"] + eukout_gtf_list
            stie_cmd_euk = stringtie[stie_cmd_euk_opt]
            stie_cmd_euk()

            stie_cmd_prok_opt = ["--merge", "-G", prok_gff,
                                 "-o", self.workdir +
                                 "/" + "/prok_sTie_merged_transcript.gtf"] + prokout_gtf_list
            stie_cmd_prok = stringtie[stie_cmd_prok_opt]
            stie_cmd_prok()


class ReStringTieScores(ExternalProgramTask):
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

    def program_args(self):
        """Run stringtie."""
        return ["stringtie",
                "-o", self.out_gtf,
                "-p", self.num_cpus,
                "-G", self.gtf,
                "-C", self.out_cover,
                "-A", self.out_abun,
                "-B",
                self.in_bam_file]


@inherits(MergeStringTies)
class ReStringTieScoresW(luigi.WrapperTask):
    """Recount to only have common transcripts."""

    gff_file = luigi.Parameter()
    kingdom = luigi.Parameter()
    workdir = luigi.Parameter()
    fastq_dic = luigi.DictParameter()
    num_cpus = luigi.Parameter()
    indexfile = luigi.Parameter()
    ref_file = luigi.Parameter()

    def requires(self):
        """A wrapper for running the QC."""
        for samp, fastq in self.fastq_dic.items():
            map_dir = self.workdir + "/" + samp + "/mapping_results"
            if self.kingdom in ['prokarya', 'eukarya']:
                bg_dir = os.path.join(self.workdir, "ballgown", self.kingdom, samp)
                if os.path.isdir(bg_dir) is False:
                    os.makedirs(bg_dir)
                yield ReStringTieScores(num_cpus=self.num_cpus,
                                        gtf=self.workdir + "/" + "stringtie_merged_transcript.gtf",
                                        out_gtf=bg_dir + "/" + samp + "_merged_sTie.gtf",
                                        out_cover=bg_dir + "/" + samp + "_merged_covered_sTie.gtf",
                                        out_abun=bg_dir + "/" + samp + "_merged_sTie.tab",
                                        in_bam_file=map_dir + "/" + samp + "_srt.bam")
            elif self.kingdom == 'both':
                euk_bg_dir = os.path.join(self.workdir, "ballgown", "eukarya", samp)
                prok_bg_dir = os.path.join(self.workdir, "ballgown", "prokarya", samp)
                if os.path.isdir(euk_bg_dir) is False:
                    os.makedirs(euk_bg_dir)
                    os.makedirs(prok_bg_dir)
                prok_gtf = os.path.join(self.workdir, "prok_sTie_merged_transcript.gtf")
                euk_gtf = os.path.join(self.workdir, "euk_sTie_merged_transcript.gtf")
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
