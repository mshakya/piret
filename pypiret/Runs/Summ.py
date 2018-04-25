#! /usr/bin/env python

"""Check design."""
from __future__ import print_function
import os
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
    bindir = luigi.Parameter()
    numCPUs = luigi.IntParameter()
    ref_file = luigi.Parameter()
    fid = luigi.Parameter()
    stranded = luigi.IntParameter()

    def output(self):
        """Expected output of featureCounts."""
        counts_dir = self.workdir + "/featureCounts"
        if ',' in self.gff:
            gff_list = self.gff.split(",")
            gff_full_path = [os.path.abspath(gff) for gff in gff_list]
            gtfs = [self.workdir + "/" + gff.split("/")[-1].split(".gff")[0] + ".gtf" for gff in gff_full_path]
            all_target = []
            for gffs in gff_full_path:
                feature = list(set(pd.read_csv(gffs, sep="\t", header=None, comment='#')[2].tolist()))
                loc_target = [LocalTarget(counts_dir +  os.path.basename(gffs) + "_" + feat + "_count.tsv") for feat in feature]
                all_target = loc_target + all_target
                return all_target
        else:
            gff_fp = os.path.abspath(self.gff)
            features = list(set(pd.read_csv(gff_fp, sep="\t", header=None, comment='#')[2].tolist()))
            features = [feat for feat in features if feat in ['CDS', 'rRNA', 'tRNA', 'exon', 'gene', 'transcript']]
            loc_target = LocalTarget(counts_dir + "/" + os.path.basename(self.gff).split(".gff")[0] + "_" + features[-1] + "_count.tsv")
            return loc_target

    def run(self):
        """Running featureCounts on all."""
        samp_list = list(self.fastq_dic.keys())
        in_srtbam_list = [self.workdir + "/" + samp + "/" +
                          "mapping_results" + "/" + samp + "_srt.bam"
                          for samp in samp_list]
        counts_dir = self.workdir + "/featureCounts"
        if not os.path.exists(counts_dir):
            os.makedirs(counts_dir)
        if ',' in self.gff:
            gff_list = self.gff.split(",")
            gff_full_path = [os.path.abspath(gff) for gff in gff_list]
            for gffs in gff_full_path:
                feature = list(set(pd.read_csv(gffs, sep="\t", header=None, comment='#')[2].tolist()))
                for feat in feature:
                    if feat in ['CDS', 'rRNA', 'tRNA', 'exon', 'gene', 'transcript']:
                        fcount_cmd_opt = ["-a", self.gff,
                                          "-s", self.stranded,
                                          "-B",
                                          "-p", "-P", "-C",
                                          "-g", self.fid,
                                          "-t", feat,
                                          "-T", self.numCPUs,
                                          "-o", counts_dir + "/" + gffs.split("/")[-1].split("gff")[0] + "_" + feat + "_count.tsv"] + in_srtbam_list
                    fcount_cmd = featureCounts[fcount_cmd_opt]
                    fcount_cmd()
        else:
            feature = list(set(pd.read_csv(self.gff, sep="\t", header=None, comment='#')[2].tolist()))
            for feat in feature:
                if feat in ['CDS', 'rRNA', 'tRNA', 'exon', 'gene', 'transcript']:
                    fcount_cmd_opt = ["-a", self.gff,
                                      "-s", self.stranded,
                                      "-B",
                                      "-p", "-P", "-C",
                                      "-g", self.fid,
                                      "-t", feat,
                                      "-T", self.numCPUs,
                                      "-o", counts_dir + "/" + self.gff.split("/")[-1].split(".gff")[0] + "_" + feat + "_count.tsv"] + in_srtbam_list
                    fcount_cmd = featureCounts[fcount_cmd_opt]
                    print(fcount_cmd)
                    fcount_cmd()


@requires(Map.HisatMapW)
class FeatureCountsBoth(luigi.Task):
    """Summarize mapped reads classificaion using FeatureCount."""

    kingdom = luigi.Parameter()
    gff_file = luigi.Parameter()
    workdir = luigi.Parameter()
    indexfile = luigi.Parameter()
    bindir = luigi.Parameter()
    numCPUs = luigi.IntParameter()
    ref_file = luigi.Parameter()
    fid = luigi.Parameter()
    stranded = luigi.IntParameter()

    def output(self):
        """Index output."""
        counts_dir = self.workdir + "/featureCounts"
        prok_features = list(set(pd.read_csv(self.gff_file.split(",")[0],
                                             sep="\t",
                                             header=None,
                                             comment="#")[2].tolist()))
        prok_target = [LocalTarget(counts_dir + "/prok_" +
                                   feat + ".count") for feat in prok_features]
        euk_features = list(set(pd.read_csv(self.gff_file.split(",")[1],
                                            sep="\t",
                                            header=None,
                                            comment="#")[2].tolist()))
        euk_target = [LocalTarget(counts_dir + "/euk_" +
                                  feat + ".count") for feat in euk_features]
        loc_target = prok_target + euk_target
        return loc_target

    def run(self):
        """Running featureCounts on all."""
        counts_dir = self.workdir + "/featureCounts"
        samp_list = list(self.fastq_dic.keys())
        in_srtbam_list = [self.workdir + "/" + samp + "/" +
                          "mapping_results" + "/" + samp + "_srt.bam"
                          for samp in samp_list]
        if not os.path.exists(self.workdir + "/featureCounts"):
            os.makedirs(self.workdir + "/featureCounts")
        prok_features = list(set(pd.read_csv(self.gff_file.split(",")[0],
                                             sep="\t",
                                             header=None,
                                             comment="#")[2].tolist()))
        euk_features = list(set(pd.read_csv(self.gff_file.split(",")[1],
                                            sep="\t",
                                            header=None,
                                            comment="#")[2].tolist()))
        for feat in euk_features:
            if feat in ['CDS', 'rRNA', 'tRNA', 'exon', 'gene', 'transcript']:
                fcount_euk_cmd_opt = ["-a", self.gff_file.split(",")[1],
                                      "-s", self.stranded,
                                      "-B",
                                      "-p", "-P", "-C",
                                      "-g", self.fid,
                                      "-t", feat,
                                      "-T", self.numCPUs,
                                      "-o",
                                      counts_dir + "/euk_" + feat +
                                      "_count.tsv"] + in_srtbam_list
            fcount_euk_cmd = featureCounts[fcount_euk_cmd_opt]
            fcount_euk_cmd()
        for feat in prok_features:
            if feat in ['CDS', 'rRNA', 'tRNA', 'exon', 'gene', 'transcript']:
                fcount_prok_cmd_opt = ["-a", self.gff_file.split(",")[0],
                                       "-s", self.stranded,
                                       "-B",
                                       "-p", "-P", "-C",
                                       "-g", self.fid,
                                       "-t", feat,
                                       "-T", self.numCPUs,
                                       "-o",
                                       counts_dir + "/prok_" + feat +
                                       "_count.tsv"] + in_srtbam_list
            fcount_prok_cmd = featureCounts[fcount_prok_cmd_opt]
            fcount_prok_cmd()

    def program_environment(self):
        """Environmental variables for this program."""
        return {'PATH': os.environ["PATH"] + ":" + self.bindir}


@requires(Map.StringTieScoresW)
class MergeStringTies(luigi.Task):
    """Summarize mapped reads classification using StringTie."""

    def ouput(self):
        """Ouptut of string tie merge."""
        if self.kingdom in ['prokarya', 'eukarya']:
            return [Map.RefFile(self.workdir + "/" + "stringtie_merged_transcript.gtf")]
        elif self.kingdom == 'both':
            return [Map.RefFile(self.workdir + "/" + "prok_stringtie_merged_transcript.gtf"),
                    Map.RefFile(self.workdir + "/" + "euk_stringtie_merged_transcript.gtf")]

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
            gff_name = gtf.split(".gtf")[0].split("/")[-1]
            out_gtf_list = [self.workdir + "/" + samp + "/" + "stie_results" +
                            "/" + samp + "_" + gff_name + append_name +
                            "_sTie.gtf" for samp in samp_list]
            stie_cmd_opt = ["--merge", "-G", self.gff_file,
                            "-o", self.workdir +
                            "/" + "stringtie_merged_transcript.gtf"] + out_gtf_list
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
                                "/" + "/euk_stringtie_merged_transcript.gtf"] + eukout_gtf_list
            stie_cmd_euk = stringtie[stie_cmd_euk_opt]
            stie_cmd_euk()

            stie_cmd_prok_opt = ["--merge", "-G", prok_gtf,
                                 "-o", self.workdir +
                                 "/" + "/prok_stringtie_merged_transcript.gtf"] + prokout_gtf_list
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


@inherits(Map.HisatMapW)
class ReStringTieScoresW(luigi.WrapperTask):
    """Recount to only have common transcripts."""

    gff_file = luigi.Parameter()
    kingdom = luigi.Parameter()
    workdir = luigi.Parameter()
    fastq_dic = luigi.DictParameter()

    def requires(self):
        """A wrapper for running the QC."""
        splist = [self.workdir + "/" +
                  f for f in os.listdir(self.workdir) if f.endswith('.splice')]
        if len(splist) > 1:
            splice_file = ','.join(splist)
        elif len(splist) == 1:
            splice_file = splist[0]
        else:
            splice_file = ''

        for samp, fastq in self.fastq_dic.items():
            map_dir = self.workdir + "/" + samp + "/mapping_results"
            trim_dir = self.workdir + "/" + samp + "/trimming_results"
            if self.kingdom in ['prokarya', 'eukarya']:
                bg_dir = self.workdir + "/ballgown/" + samp
                if os.path.isdir(bg_dir) is False:
                    os.makedirs(bg_dir)
                yield ReStringTieScores(fastqs=[trim_dir + "/" + samp +
                                                ".1.trimmed.fastq",
                                                trim_dir + "/" + samp +
                                                ".2.trimmed.fastq"],
                                        numCPUs=self.numCPUs,
                                        indexfile=self.indexfile,
                                        spliceFile=splice_file,
                                        outsam=map_dir + "/" + samp + ".sam",
                                        ref_file=self.ref_file,
                                        gff_file=self.gff_file,
                                        gtf=self.workdir + "/" + "stringtie_merged_transcript.gtf",
                                        out_gtf=bg_dir + "/" + samp + "_merged_sTie.gtf",
                                        out_cover=bg_dir + "/" + samp + "_merged_covered_sTie.gtf",
                                        out_abun=bg_dir + "/" + samp + "_merged_sTie.tab",
                                        in_bam_file=map_dir + "/" + samp + "_srt.bam",
                                        bindir=self.bindir,
                                        workdir=self.workdir,
                                        sample=samp,
                                        qc_outdir=trim_dir,
                                        map_dir=map_dir)
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
                                        outsam=self.map_dir + "/" + samp + ".sam",
                                        ref_file=self.ref_file,
                                        gff_file=self.gff_file,
                                        gtf=prok_gtf,
                                        out_gtf=bg_dir_prok + "/" + samp + "_prok_sTie.gtf",
                                        out_cover=bg_dir_prok + "/" + samp + "_prok_covered_sTie.gtf",
                                        out_abun=bg_dir_prok + "/" + samp + "_prok_sTie.tab",
                                        in_bam_file=self.map_dir + "/prokarya.bam",
                                        bindir=self.bindir,
                                        workdir=self.workdir)
                yield ReStringTieScores(fastq1=trim_dir + "/" + samp + ".1.trimmed.fastq",
                                        fastq2=trim_dir + "/" + samp + ".2.trimmed.fastq",
                                        numCPUs=self.numCPUs,
                                        indexfile=self.indexfile,
                                        spliceFile=splice_file,
                                        outsam=self.map_dir + "/" + samp + ".sam",
                                        ref_file=self.ref_file,
                                        gff_file=self.gff_file,
                                        gtf=euk_gtf,
                                        out_gtf=bg_dir_euk + "/" + samp + "_euk_sTie.gtf",
                                        out_cover=bg_dir_euk + "/" + samp + "_euk_covered_sTie.gtf",
                                        out_abun=bg_dir_euk + "/" + samp + "_euk_sTie.tab",
                                        in_bam_file=self.map_dir + "/eukarya.bam",
                                        bindir=self.bindir,
                                        workdir=self.workdir)

# class SampleInfo(luigi.task):
#     luigi.Parameter() = workdir

#     def run(self):
#         for root, dirs, files in os.walk(workdir):

            # if file.endswith
