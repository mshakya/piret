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
    mapp_file = luigi.Parameter()
    bindir = luigi.Parameter()

    def requires(self):
        """Require reference fasta format file."""
        return [RefFile(self.gff), RefFile(self.mapp_file)]

    def output(self):
        """Index output."""
        hisat_index_ht8l = self.hi_index + ".8.ht2l"
        return LocalTarget(hisat_index_ht8l)

    def program_args(self):
        """Run hisat2-build command."""
        return ["featureCounts",
                "-a", self.gff,
                "-s", 1,
                "-p", "-P",
                self.mapp_file,
                "-o", self.count_file,
                self.mapp_file]

    def program_environment(self):
        """Environmental variables for this program."""
        return {'PATH': os.environ["PATH"] + ":" + self.bindir}
