#! /usr/bin/env python

"""Check design."""
from __future__ import print_function
import os
import luigi
from luigi.contrib.external_program import ExternalProgramTask
from luigi.contrib.external_program import ExternalPythonProgramTask
from luigi import ExternalTask
from luigi import LocalTarget
from luigi import Parameter


class RefFile(ExternalTask):
    """An ExternalTask like this."""

    path = Parameter()

    def output(self):
        """Check."""
        return LocalTarget(os.path.abspath(self.path))


class ConcatF(ExternalProgramTask):
    """Create HisatIndex."""

    prok_ref = luigi.Parameter()
    euk_ref = luigi.Parameter()
    ref = luigi.Parameter()

    def requires(self):
        """Require reference fasta format file."""
        return [RefFile(self.prok_ref), RefFile(self.euk_ref)]

    def output(self):
        """Index output."""
        return luigi.LocalTarget(self.ref)

    def program_args(self):
        """Run hisat2-build command."""
        return ["cat", self.prok_ref, self.euk_ref, '>', self.ref]


class HisatIndex(ExternalProgramTask):
    """Create HisatIndex."""

    fasta = Parameter()
    hi_index = Parameter()
    bindir = Parameter()

    def requires(self):
        """Require reference fasta format file."""
        if ',' in self.fasta():
            fas = self.fasta.split(",")
            return [RefFile(fa) for fa in fas]
        else:
            return [RefFile(self.fasta)]

    def output(self):
        """Index output."""
        hisat_index_ht8l = self.hi_index + ".8.ht2l"
        return LocalTarget(hisat_index_ht8l)

    def program_args(self):
        """Run hisat2-build command."""
        return ["hisat2-build",
                "--large-index",
                "-q", self.fasta, self.hi_index]

    def program_environment(self):
        """Environmental variables for this program."""
        return {'PATH': os.environ["PATH"] + ":" + self.bindir}


class GFF2GTF(ExternalProgramTask):
    """Convert to gtf format."""

    gff_file = Parameter()
    out_gtf = Parameter()
    bindir = Parameter()

    def requires(self):
        """Require reference gff file."""
        return [RefFile(self.gff_file)]

    def output(self):
        """Output converted gtf file."""
        return LocalTarget(self.out_gtf)

    def program_args(self):
        """The gffread conversion."""
        return ["gffread", self.gff_file,
                "-T", "-o", self.out_gtf]

    def program_environment(self):
        """Environmental variables for this program."""
        return {'PATH': os.environ["PATH"] + ":" + self.bindir}


class CreateSplice(ExternalPythonProgramTask):
    """Find splice sites off gtf file."""

    gff_file = Parameter()
    euk_gtf = Parameter()
    out_splice = Parameter()
    bindir = Parameter()

    def requires(self):
        """Require pair of fastq."""
        return [GFF2GTF(self.gff_file, self.euk_gtf, self.bindir)]

    def output(self):
        """Splice site output."""
        return LocalTarget(self.out_splice)

    def program_args(self):
        """Calling hisat2_extract_splice_sites."""
        return ["scripts/hisat2_extract_splice_sites.py",
                self.euk_gtf, ">", self.out_splice]


class MapHisat(ExternalProgramTask):
    """Running the perl script for QC."""

    fastq1 = Parameter()
    fastq2 = Parameter()
    numCPUs = luigi.IntParameter()
    indexfile = Parameter()
    spliceFile = Parameter()
    mappingLogFile = Parameter()
    outsam = Parameter()
    unalned = Parameter()

    def requires(self):
        """Require pair of fastq and index files."""
        return [RefFile(self.fastq1), RefFile(self.fastq2),
                RefFile(self.indexfile + ".8.ht2l")]

    def output(self):
        """QC output."""
        return LocalTarget(self.outsam)

    def program_args(self):
        """Run hisat2."""
        if self.spliceFile is None:
            return ["hisat2",
                    "-p", self.numCPUs,
                    "-x", self.indexfile,
                    "-1", self.fastq1,
                    "-2", self.fastq2,
                    "--un-conc", self.unalned,
                    "2>", self.mappingLogFile,
                    ">", self.outsam]
        else:
            return ["hisat2",
                    "--known-splicesite-infile", self.spliceFile,
                    "-p", self.numCPUs,
                    "-x", self.indexfile,
                    "-1", self.fastq1,
                    "-2", self.fastq2,
                    "--un-conc", self.unalned,
                    "2>", self.mappingLogFile,
                    ">", self.outsam]


class RunAllProkMap(luigi.WrapperTask):
    """Run all Mapping steps from Prokarya genome."""

    fastq_dic = luigi.DictParameter()
    gff_file = Parameter()
    gtf_file = Parameter()
    numCPUs = luigi.IntParameter()
    indexfile = Parameter()
    bindir = Parameter()
    workdir = Parameter()

    def requires(self):
        """A wrapper for running mapping of proks."""
        yield GFF2GTF(gff_file=self.gff_file, out_gtf=self.gtf_file,
                      bindir=self.bindir)
        for samp, fastq in self.fastq_dic.iteritems():
            trim_dir = self.workdir + "/" + samp + "/trimming_results"
            map_dir = self.workdir + "/" + samp + "/mapping_results"
            if os.path.isdir(map_dir) is False:
                os.makedirs(map_dir)
            yield MapHisat(fastq1=trim_dir + "/" + samp +
                           ".1.trimmed.fastq",
                           fastq2=trim_dir + "/" + samp +
                           ".2.trimmed.fastq",
                           numCPUs=self.numCPUs,
                           spliceFile=None,
                           indexfile=self.indexfile,
                           unalned=map_dir + "/unligned.fastq",
                           mappingLogFile=map_dir + "/mapping.log",
                           outsam=map_dir + "/" + samp + ".mapped.sam")


class RunAllEukMap(luigi.WrapperTask):
    """Run all Mapping steps."""

    gff_file = Parameter()
    gtf_file = Parameter()
    fastq_dic = luigi.DictParameter()
    numCPUs = luigi.IntParameter()
    indexfile = Parameter()
    bindir = Parameter()
    workdir = Parameter()

    def requires(self):
        """A wrapper for running the QC."""
        yield GFF2GTF(gff_file=self.gff_file,
                      out_gtf=self.gtf_file,
                      bindir=self.bindir)
        yield CreateSplice(gff_file=self.gff_file,
                           euk_gtf=self.gtf_file,
                           bindir=self.bindir,
                           out_splice=self.workdir + "/splice_sites.txt")
        for samp, fastq in self.fastq_dic.iteritems():
            trim_dir = self.workdir + "/" + samp + "/trimming_results"
            map_dir = self.workdir + "/" + samp + "/mapping_results"
            yield MapHisat(fastq1=trim_dir + "/" + samp + ".1.trimmed.fastq",
                           fastq2=trim_dir + "/" + samp + ".2.trimmed.fastq",
                           numCPUs=self.numCPUs,
                           indexfile=self.indexfile,
                           spliceFile=self.workdir + "/splice_sites.txt",
                           mappingLogFile=map_dir + "/mapping.log",
                           unalned=map_dir + "/unligned.fastq",
                           outsam=map_dir + "/" + samp + ".mapped.sam")
