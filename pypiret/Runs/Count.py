#! /usr/bin/env python

"""Check design."""
from __future__ import print_function
import os
import luigi
from luigi.contrib.external_program import ExternalProgramTask
from luigi import LocalTarget


class RefFile(luigi.ExternalTask):
    """An ExternalTask like this."""

    path = luigi.Parameter()

    def output(self):
        """Check."""
        return LocalTarget(os.path.abspath(self.path))


class FeatureCount(ExternalProgramTask):
    """Create HisatIndex."""

    gff = luigi.Parameter()
    count_file = luigi.Parameter()
    bam_file = luigi.Parameter()
    bindir = luigi.Parameter()

    def requires(self):
        """Require reference fasta format file."""
        return [RefFile(self.gff), RefFile(self.bam_file)]

    def output(self):
        """Index output."""
        hisat_index_ht8l = self.hi_index + ".8.ht2l"
        return LocalTarget(hisat_index_ht8l)

    def program_args(self):
        """Run hisat2-build command."""
        return ["bin/featureCounts",
                "-a", self.gff,
                "-s", 1,
                "-p", "-P",
                "-o", self.count_file,
                self.bam_file]

    def program_environment(self):
        """Environmental variables for this program."""
        return {'PATH': os.environ["PATH"] + ":" + self.bindir}


class CountAll(luigi.WrapperTask):
    """Run all QC."""

    # fastq_dic = luigi.DictParameter()
    gff = luigi.Parameter()
    bam_file = luigi.Parameter()
    bindir = luigi.Parameter()
    kingdom = luigi.Parameter()

    def requires(self):
        """A wrapper for counting."""
        for samp in self.fastq_dic.iteritems():
            map_dir = self.workdir + "/" + samp + "/mapping_results"
            if self.kingdom == 'prokarya':
                yield FeatureCount(gff=self.gff,
                                   count_file=map_dir + "/" + samp + "_prok.count",
                                   bam_file=self.bam_file,
                                   bindir=self.bindir)
            elif self.kingdom == 'eukarya':
                yield FeatureCount(gff=self.gff,
                                   count_file=map_dir + "/" + samp + "_euk.count",
                                   bam_file=self.bam_file,
                                   bindir=self.bindir)
