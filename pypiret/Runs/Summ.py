#! /usr/bin/env python

"""Check design."""
from __future__ import print_function
import os
import luigi
from luigi.contrib.external_program import ExternalProgramTask
from luigi import LocalTarget
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
            prok_gtf = self.workdir + "/" + self.prok_gff.split("/")[-1].split(".gff")[0] + ".gtf"
            bam_filelist = [self.workdir + "/" + samp + "/" +
                            "mapping_results" + "/" +
                            "prokarya.bam" for samp in samp_list]
            return ["featureCounts",
                    "-a", prok_gtf,
                    "-s", 1,
                    "-B"
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


# class MergeStringTies(ExternalProgramTask):
#     """Summarize mapped reads classification using StringTie."""

#     kingdom = luigi.Parameter()
#     prok_gtf = luigi.Parameter()
#     euk_gtf = luigi.Parameter()
#     samp_list = luigi.listParameter()


#     def requires(self):
#         """Requ."""

#     def ouput(self):
#         """Mig."""

#     def program_args(self):
#         return ["stringtie",
#                 "-o", self.out_gtf,
#                 "-G", euk_gff,
#                 "-C", self.out_cover,
#                 "-A", self.out_abun,
#                 self.bam_file]


#     def program_environment(self):
