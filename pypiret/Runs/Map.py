#! /usr/bin/env python

"""Luigi Tasks to perform various RNA seq functions, from Mapping to Counting.

Mapping is done using hisat2 and counting is done using featurecounts and stringtie
"""
from __future__ import print_function
import os
import luigi
from luigi.contrib.external_program import ExternalProgramTask
from luigi import ExternalTask
from luigi import LocalTarget
from luigi import Parameter
from luigi.util import inherits, requires
import subprocess
from pypiret import FastQC
# from Bio import SeqIO


class RefFile(ExternalTask):
    """An ExternalTask like this."""

    path = Parameter()

    def output(self):
        """Check."""
        return LocalTarget(os.path.abspath(self.path))


class GetChromName(luigi.Task):
    """Extract chromosomes of euk and prok fasta."""

    prok_ref = Parameter()
    euk_ref = Parameter()
    workdir = Parameter()

    def requires(self):
        """Require slist of rerence sequence to be made."""
        return [RefFile(self.euk_ref), RefFile(self.prok_ref)]

    def output(self):
        """Split Files."""
        prok_chrom_file = self.workdir + "/" + "prok.chroms"
        euk_chrom_file = self.workdir + "/" + "euk.chroms"
        return [luigi.LocalTarget(prok_chrom_file),
                luigi.LocalTarget(euk_chrom_file)]

    def run(self):
        """Run the subprocess for grep and sed."""
        prok_chrom_file = self.workdir + "/" + "prok.chroms"
        prok_cmd = "grep '>' %s | sed 's/>//g' | sed 's/ .*$//g' > %s" % (self.prok_ref,
                                                                          prok_chrom_file)
        subprocess.Popen(prok_cmd, shell=True)

        euk_chrom_file = self.workdir + "/" + "euk.chroms"
        euk_cmd = "grep '>' %s | sed 's/>//g' | sed 's/ .*$//g' > %s" % (self.euk_ref,
                                                                         euk_chrom_file)
        subprocess.Popen(euk_cmd, shell=True)


class HisatIndex(ExternalProgramTask):
    """Create Hisat Indeices from given fasta file."""

    fasta = Parameter()
    hi_index = Parameter()
    bindir = Parameter()
    numCPUs = Parameter()

    def requires(self):
        """Require reference fasta format file."""
        if ',' in self.fasta:
            fas = self.fasta.split(",")
            return [RefFile(fa) for fa in fas]
        else:
            return [RefFile(self.fasta)]

    def output(self):
        """Expected index output."""
        hisat_index_ht8l = self.hi_index + ".8.ht2l"
        return LocalTarget(hisat_index_ht8l)

    def program_args(self):
        """Run hisat2-build command."""
        return ["hisat2-build",
                "--large-index",
                "-q", "-p", self.numCPUs,
                self.fasta, self.hi_index]

    def program_environment(self):
        """Environmental variables for this program."""
        return {'PATH': os.environ["PATH"] + ":" + self.bindir}


class GFF2GTF(ExternalProgramTask):
    """Convert GFF to GTF format."""

    gff_file = Parameter()
    workdir = Parameter()
    bindir = Parameter()

    def requires(self):
        """Require reference gff(s) file."""
        if ',' in self.gff_file:
            gffs = self.gff_file.split(",")
            return [RefFile(gff) for gff in gffs]
        else:
            return [RefFile(self.gff_file)]

    def output(self):
        """Output converted gtf file."""
        if ',' in self.gff_file:
            gffs = self.gff_file.split(",")
            return [LocalTarget(self.workdir + "/" +
                                gff.split("/")[-1].split(".")[0] +
                                ".gtf") for gff in gffs]
        else:
            return LocalTarget(self.workdir + "/" +
                               self.gff_file.split("/")[-1].split(".")[0] +
                               ".gtf")

    def program_args(self):
        """The gffread conversion."""
        if ',' in self.gff_file:
            gffs = self.gff_file.split(",")
            for gff in gffs:
                out_file = self.workdir + "/" +\
                    self.gff_file.split("/")[-1].split(".")[0] + ".gtf"
                return ["gffread", gff, "-T", "-o", out_file]
        else:
            out_file = self.workdir + "/" +\
                self.gff_file.split("/")[-1].split(".")[0] + ".gtf"
            return ["gffread", self.gff_file, "-T", "-o", out_file]

    def program_environment(self):
        """Environmental variables for this program."""
        return {'PATH': os.environ["PATH"] + ":" + self.bindir}


@inherits(GFF2GTF)
class CreateSplice(ExternalProgramTask):
    """Find splice sites off gtf file."""

    scriptdir = Parameter()

    def requires(self):
        """Require GTF file."""
        if ',' in self.gff_file:
            gffs = self.gff_file.split(",")
            return [GFF2GTF(gff_file=gff,
                            workdir=self.workdir,
                            bindir=self.bindir) for gff in gffs]
        else:
            return [GFF2GTF(gff_file=self.gff_file,
                            workdir=self.workdir,
                            bindir=self.bindir)]

    def output(self):
        """Splice site output."""
        if ',' in self.gff_file:
            gffs = self.gff_file.split(",")
            return [LocalTarget(self.workdir + "/" +
                                gff.split("/")[-1].split(".")[0] +
                                ".splice") for gff in gffs]
        else:
            return LocalTarget(self.workdir + "/" +
                               self.gff_file.split("/")[-1].split(".")[0] +
                               ".splice")

    def program_args(self):
        """Extracting splice sites."""
        if ',' in self.gff_file:
            gffs = self.gff_file.split(",")
            for gff in gffs:
                gtf_file = self.workdir + "/" +\
                    gff.split("/")[-1].split(".")[0] + ".gtf"
                out_file = self.workdir + "/" +\
                    gff.split("/")[-1].split(".")[0] + ".splice"
            return [self.scriptdir + "/hisat2_extract_splice_sites.py",
                    "-i", gtf_file, "-o", out_file]
        else:
            gtf_file = self.workdir + "/" +\
                self.gff_file.split("/")[-1].split(".")[0] + ".gtf"

            return [self.scriptdir + "/hisat2_extract_splice_sites.py",
                    "-i", gtf_file, "-o", self.workdir + "/" +
                    self.gff_file.split("/")[-1].split(".")[0] + ".splice"]

    def program_environment(self):
        """Environmental variables for this program."""
        return {'PATH': os.environ["PATH"] + ":" + self.bindir + "../scripts"}


# @inherits(FastQC.PairedRunQC)
class Hisat(ExternalProgramTask):
    """Mapping the QCed sequences to reference."""

    fastq1 = Parameter()
    fastq2 = Parameter()
    numCPUs = luigi.IntParameter()
    indexfile = Parameter()
    spliceFile = Parameter()
    mappingLogFile = Parameter()
    outsam = Parameter()
    unalned = Parameter()
    scriptdir = Parameter()
    ref_file = Parameter()
    bindir = Parameter()

    def requires(self):
        """Require pair of fastq and index files."""
        return[FastQC.PairedRunQC(fastqs=[self.fastq1, self.fastq2],
                                  sample=self.outsam.split(
                                      ".sam")[0].split("/")[-1],
                                  numCPUs=self.numCPUs,
                                  outdir=self.mappingLogFile.split("mapping_results")[0] +
                                  "trimming_results",
                                  scriptdir=self.scriptdir),
               HisatIndex(fasta=self.ref_file,
                          hi_index=self.indexfile,
                          bindir=self.bindir,
                          numCPUs=self.numCPUs)]

    def output(self):
        """SAM file output of the mapping."""
        return LocalTarget(self.outsam)

    def program_args(self):
        """Run hisat2."""
        if self.spliceFile is "":
            return ["hisat2",
                    "-p", self.numCPUs,
                    "-x", self.indexfile,
                    "-1", self.fastq1,
                    "-2", self.fastq2,
                    "-S", self.outsam,
                    "--un-conc", self.unalned,
                    "2>", self.mappingLogFile]
        else:
            return ["hisat2",
                    "--known-splicesite-infile", self.spliceFile,
                    "-p", self.numCPUs,
                    "-x", self.indexfile,
                    "-1", self.fastq1,
                    "-2", self.fastq2,
                    "-S", self.outsam,
                    "--un-conc", self.unalned,
                    "2>", self.mappingLogFile]


@inherits(FastQC.RunAllQC)
class HisatMapW(luigi.WrapperTask):
    """A wrapper task for running mapping."""

    ref_file = luigi.Parameter()
    indexfile = luigi.Parameter()
    bindir = luigi.Parameter()

    def requires(self):
        """A wrapper task for running mapping."""
        splice_list = [self.workdir + "/" +
                       f for f in os.listdir(self.workdir) if f.endswith('.splice')]
        if len(splice_list) > 1:
            splice_file = ','.join([splice_list])
        elif len(splice_list) == 1:
            splice_file = splice_list[0]
        else:
            splice_file = ''
        for samp, fastq in self.fastq_dic.iteritems():
            trim_dir = self.workdir + "/" + samp + "/trimming_results"
            map_dir = self.workdir + "/" + samp + "/mapping_results"
            if os.path.isdir(map_dir) is False:
                os.makedirs(map_dir)
            yield Hisat(fastq1=trim_dir + "/" + samp + ".1.trimmed.fastq",
                        fastq2=trim_dir + "/" + samp + ".2.trimmed.fastq",
                        numCPUs=self.numCPUs,
                        indexfile=self.indexfile,
                        spliceFile=splice_file,
                        mappingLogFile=map_dir + "/mapping.log",
                        unalned=map_dir + "/unligned.fastq",
                        outsam=map_dir + "/" + samp + ".sam",
                        scriptdir=self.scriptdir,
                        ref_file=self.ref_file,
                        bindir=self.bindir)


@inherits(HisatMapW)
class HiSatBoth(luigi.WrapperTask):
    """Mapping."""

    def requires(self):
        """Mapping reads."""
        splice_list = [self.workdir + "/" +
                       f for f in os.listdir(self.workdir) if f.endswith('.splice')]
        print(splice_list)
        if len(splice_list) > 1:
            splice_file = ','.join(splice_list)
        elif len(splice_list) == 1:
            splice_file = splice_list[0]
        else:
            splice_file = ''

        for samp, fastq in self.fastq_dic.iteritems():
            trim_dir = self.workdir + "/" + samp + "/trimming_results"
            map_dir = self.workdir + "/" + samp + "/mapping_results"
            if os.path.isdir(map_dir) is False:
                os.makedirs(map_dir)
            yield Hisat(fastq1=trim_dir + "/" + samp + ".1.trimmed.fastq",
                        fastq2=trim_dir + "/" + samp + ".2.trimmed.fastq",
                        numCPUs=self.numCPUs,
                        indexfile=self.indexfile,
                        spliceFile=splice_file,
                        mappingLogFile=map_dir + "/mapping.log",
                        unalned=map_dir + "/unligned.fastq",
                        outsam=map_dir + "/" + samp + ".sam",
                        scriptdir=self.scriptdir,
                        ref_file=self.ref_file,
                        bindir=self.bindir)


@requires(Hisat)
class SAM2BAMfile(ExternalProgramTask):
    """Convert SAM file to BAM file with only mapped reads."""

    bam_file = Parameter()

    def output(self):
        """Output BAM file."""
        return LocalTarget(self.bam_file)

    def program_args(self):
        """Output BAM file with only mapped reads."""
        return ["samtools", "view", "-bSh", "-F",
                "4", self.outsam, "-o", self.bam_file]


@inherits(HisatMapW)
class SAM2BAMfileW(luigi.WrapperTask):
    """From Mapping to Counting step for Eukaryotic reference."""

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
            trim_dir = self.workdir + "/" + samp + "/trimming_results"
            map_dir = self.workdir + "/" + samp + "/mapping_results"
            if os.path.isdir(map_dir) is False:
                os.makedirs(map_dir)
            yield SAM2BAMfile(fastq1=trim_dir + "/" + samp + ".1.trimmed.fastq",
                              fastq2=trim_dir + "/" + samp + ".2.trimmed.fastq",
                              numCPUs=self.numCPUs,
                              indexfile=self.indexfile,
                              spliceFile=splice_file,
                              mappingLogFile=map_dir + "/mapping.log",
                              unalned=map_dir + "/unligned.fastq",
                              outsam=map_dir + "/" + samp + ".sam",
                              bam_file=map_dir + "/" + samp + ".bam",
                              scriptdir=self.scriptdir,
                              ref_file=self.ref_file,
                              bindir=self.bindir)


@requires(SAM2BAMfile)
class SortBAMfile(ExternalProgramTask):
    """Sort BAM file."""

    sorted_bam_file = Parameter()

    def output(self):
        """Sorted BAM file."""
        return LocalTarget(self.sorted_bam_file)

    def program_args(self):
        """Sort BAM file."""
        return ["samtools", "sort",
                self.bam_file, "-o", self.sorted_bam_file]


@inherits(SAM2BAMfileW)
class SortBAMfileW(luigi.WrapperTask):
    """Sort all bam files."""

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
            trim_dir = self.workdir + "/" + samp + "/trimming_results"
            map_dir = self.workdir + "/" + samp + "/mapping_results"
            if os.path.isdir(map_dir) is False:
                os.makedirs(map_dir)

            yield SortBAMfile(fastq1=trim_dir + "/" + samp + ".1.trimmed.fastq",
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
                              bindir=self.bindir)


@requires(SortBAMfile)
class RefNames(luigi.Task):
    """Extract the name of chromosomes where reads were mapped in BAM file."""

    map_ref = Parameter()

    def output(self):
        """Split Files."""
        return luigi.LocalTarget(self.map_ref)

    def run(self):
        """Extract names of chromosome from SAM file."""
        grep_cmd = "samtools view -h %s |\
                    grep -v -e '^.SQ' |\
                    grep -v -e 'PN:' |\
                    grep -v -e 'SO:' |\
                    awk -F '\t' '{print $3}' |\
                    sort | uniq |\
                    grep -v '*' > %s" % (self.sorted_bam_file, self.map_ref)
        subprocess.Popen(grep_cmd, shell=True)


@inherits(SortBAMfileW)
class GetRefNames(luigi.WrapperTask):
    """From Mapping to Counting step for Eukaryotic and Prokaryotic."""

    def requires(self):
        """A pipeline that runs from mapping to count for euk and prok."""
        splice_list = [self.workdir + "/" +
                       f for f in os.listdir(self.workdir) if f.endswith('.splice')]
        if len(splice_list) > 1:
            splice_file = ','.join([splice_list])
        else:
            splice_file = splice_list[0]
        for samp, fastq in self.fastq_dic.iteritems():
            trim_dir = self.workdir + "/" + samp + "/trimming_results"
            map_dir = self.workdir + "/" + samp + "/mapping_results"
            yield RefNames(fastq1=trim_dir + "/" + samp + ".1.trimmed.fastq",
                           fastq2=trim_dir + "/" + samp + ".2.trimmed.fastq",
                           numCPUs=self.numCPUs,
                           indexfile=self.indexfile,
                           spliceFile=splice_file,
                           mappingLogFile=map_dir + "/mapping.log",
                           unalned=map_dir + "/unligned.fastq",
                           outsam=map_dir + "/" + samp + ".sam",
                           bam_file=map_dir + "/" + samp + ".bam",
                           sorted_bam_file=map_dir + "/" + samp + "_srt.bam",
                           map_ref=map_dir + "/" + samp + ".reflist",
                           scriptdir=self.scriptdir,
                           ref_file=self.ref_file,
                           bindir=self.bindir)


@inherits(SortBAMfile)
class StringTieScores(ExternalProgramTask):
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
                self.in_bam_file]


@inherits(SortBAMfileW)
class StringTieScoresW(luigi.WrapperTask):
    """From Mapping to Counting step for Eukaryotic reference."""

    gff = Parameter()
    kingdom = Parameter()

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

        gtf = self.workdir + "/" + \
            self.gff.split("/")[-1].split(".gff")[0] + ".gtf"
        for samp, fastq in self.fastq_dic.iteritems():
            map_dir = self.workdir + "/" + samp + "/mapping_results"
            trim_dir = self.workdir + "/" + samp + "/trimming_results"
            if os.path.isdir(map_dir) is False:
                os.makedirs(map_dir)
            if self.kingdom in ['prokarya', 'eukarya']:

                yield StringTieScores(fastq1=trim_dir + "/" + samp + ".1.trimmed.fastq",
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
                                      out_gtf=map_dir + "/" + samp + "_sTie.gtf",
                                      out_cover=map_dir + "/" + samp + "_covered_sTie.gtf",
                                      out_abun=map_dir + "/" + samp + "_sTie.tab",
                                      in_bam_file=map_dir + "/" + samp + "_srt.bam",
                                      bindir=self.bindir)
            elif self.kingdom == 'both':
                yield StringTieScores(fastq1=trim_dir + "/" + samp + ".1.trimmed.fastq",
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
                                      out_gtf=map_dir + "/" + samp + "_prok_sTie.gtf",
                                      out_cover=map_dir + "/" + samp + "_prok_covered_sTie.gtf",
                                      out_abun=map_dir + "/" + samp + "_prok_sTie.tab",
                                      in_bam_file=map_dir + "/prokarya.bam",
                                      bindir=self.bindir)
                yield StringTieScores(fastq1=trim_dir + "/" + samp + ".1.trimmed.fastq",
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
                                      out_gtf=map_dir + "/" + samp + "_euk_sTie.gtf",
                                      out_cover=map_dir + "/" + samp + "_euk_covered_sTie.gtf",
                                      out_abun=map_dir + "/" + samp + "_euk_sTie.tab",
                                      in_bam_file=map_dir + "/eukarya.bam",
                                      bindir=self.bindir)


@requires(RefNames)
class SplitBAMfile(luigi.Task):
    """Split BAM file to individual chromosomes."""

    split_logfile = Parameter()

    def get_file_path(self):
        """Get bam filename."""
        with open(self.map_ref, 'r') as mrf:
            bam = self.sorted_bam_file.split(
                ".")[0] + ".REF_" + str(mrf.readline()) + ".bam"
        return bam

    def output(self):
        """Split Files."""
        return luigi.LocalTarget(self.split_logfile)

    def run(self):
        """Split BAM file based on chromosome."""
        bam_cmd = "bamtools split -reference -in %s 2> %s" % (
            self.sorted_bam_file, self.split_logfile)
        subprocess.Popen(bam_cmd, shell=True)


@inherits(GetRefNames)
class SplitBAMBoth(luigi.WrapperTask):
    """From Mapping to Counting step for Eukaryotic and Prokaryotic."""

    def requires(self):
        """A pipeline that runs from mapping to count for euk and prok."""
        splice_list = [self.workdir + "/" +
                       f for f in os.listdir(self.workdir) if f.endswith('.splice')]
        if len(splice_list) > 1:
            splice_file = ','.join([splice_list])
        else:
            splice_file = splice_list[0]
        for samp, fastq in self.fastq_dic.iteritems():
            trim_dir = self.workdir + "/" + samp + "/trimming_results"
            map_dir = self.workdir + "/" + samp + "/mapping_results"

            yield SplitBAMfile(fastq1=trim_dir + "/" + samp + ".1.trimmed.fastq",
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
                               map_ref=map_dir + "/" + samp + ".reflist",
                               ref_file=self.ref_file,
                               bindir=self.bindir,
                               split_logfile=map_dir + "/" + samp + ".splitlog")


@requires(RefNames)
class SplitRefProkEuk(luigi.Task):
    """Create a input BAM file list of prok and euk."""

    map_dir = Parameter()
    workdir = Parameter()

    def output(self):
        """Two files, with full paths to prok and euk."""
        euk_chroms_fp = self.map_dir + "/" + "eukarya_chromos.fullpath"
        prok_chroms_fp = self.map_dir + "/" + "prokarya_chromos.fullpath"
        return [luigi.LocalTarget(euk_chroms_fp), luigi.LocalTarget(prok_chroms_fp)]

    def run(self):
        """Split BAM file based on chromosome."""
        ref_bams = [f for f in os.listdir(self.map_dir) if 'REF_' in f]
        prok_ref_file = self.workdir + '/' + 'prok.chroms'
        prok_refs = [line.rstrip() for line in open(prok_ref_file, 'r')]
        euk_ref_file = self.workdir + '/' + 'euk.chroms'
        euk_refs = [line.rstrip() for line in open(euk_ref_file, 'r')]
        euk_list = []
        prok_list = []
        for ref in ref_bams:
            if ref.split('REF_')[1].split('.bam')[0] in prok_refs:
                ref1 = self.map_dir + "/" + ref
                prok_list.append(ref1)
            elif ref.split('REF_')[1].split('.bam')[0] in euk_refs:
                ref1 = self.map_dir + "/" + ref
                euk_list.append(ref1)
        euk_chroms_fp = self.map_dir + "/" + "eukarya_chromos.fullpath"
        prok_chroms_fp = self.map_dir + "/" + "prokarya_chromos.fullpath"
        with open(euk_chroms_fp, 'w') as ecf:
            for f in euk_list:
                ecf.write("%s\n" % f)
        with open(prok_chroms_fp, 'w') as pcf:
            for f in prok_list:
                pcf.write("%s\n" % f)


@inherits(GetRefNames)
class SplitProkEukBoth(luigi.WrapperTask):
    """Group chromosomes to prokaryotic and eukaryotic."""

    def requires(self):
        """A pipeline that runs from mapping to count for euk and prok."""
        splice_list = [self.workdir + "/" +
                       f for f in os.listdir(self.workdir) if f.endswith('.splice')]
        if len(splice_list) > 1:
            splice_file = ','.join([splice_list])
        else:
            splice_file = splice_list[0]
        for samp, fastq in self.fastq_dic.iteritems():
            trim_dir = self.workdir + "/" + samp + "/trimming_results"
            map_dir = self.workdir + "/" + samp + "/mapping_results"
            yield SplitRefProkEuk(fastq1=trim_dir + "/" + samp + ".1.trimmed.fastq",
                                  fastq2=trim_dir + "/" + samp + ".2.trimmed.fastq",
                                  numCPUs=self.numCPUs,
                                  indexfile=self.indexfile,
                                  spliceFile=splice_file,
                                  mappingLogFile=map_dir + "/mapping.log",
                                  unalned=map_dir + "/unligned.fastq",
                                  outsam=map_dir + "/" + samp + ".sam",
                                  bam_file=map_dir + "/" + samp + ".bam",
                                  sorted_bam_file=map_dir + "/" + samp + "_srt.bam",
                                  map_ref=map_dir + "/" + samp + ".reflist",
                                  map_dir=map_dir,
                                  scriptdir=self.scriptdir,
                                  ref_file=self.ref_file,
                                  bindir=self.bindir,
                                  workdir=self.workdir)


@requires(SplitRefProkEuk)
class MergeBAMfile(ExternalProgramTask):
    """Merge BAM file to prokaryotic and Eukaryotic."""

    kingdom = Parameter()

    def output(self):
        """Split Files."""
        if self.kingdom == 'prokarya':
            prok_bam = self.map_dir + "/" + "prokarya.bam"
            return luigi.LocalTarget(prok_bam)
        elif self.kingdom == 'eukarya':
            euk_bam = self.map_dir + "/" + "eukarya.bam"
            return luigi.LocalTarget(euk_bam)

    def program_args(self):
        """Merge BAM files to euk and prok."""
        if self.kingdom == 'prokarya':
            prok_list = self.map_dir + "/" + "prokarya_chromos.fullpath"
            prok_bam = self.map_dir + "/" + "prokarya.bam"
            return["bamtools", "merge", "-list", prok_list,
                   "-out", prok_bam]
        elif self.kingdom == 'eukarya':
            euk_list = self.map_dir + "/" + "eukarya_chromos.fullpath"
            euk_bam = self.map_dir + "/" + "eukarya.bam"
            return["bamtools", "merge", "-list",
                   euk_list, "-out", euk_bam]


@inherits(SplitProkEukBoth)
class MergeBAMfileBoth(luigi.WrapperTask):
    """From Mapping to Counting step for Eukaryotic and Prokaryotic."""

    def requires(self):
        """A pipeline that runs from mapping to count for euk and prok."""
        splice_list = [self.workdir + "/" +
                       f for f in os.listdir(self.workdir) if f.endswith('.splice')]
        if len(splice_list) > 1:
            splice_file = ','.join([splice_list])
        else:
            splice_file = splice_list[0]
        for samp, fastq in self.fastq_dic.iteritems():
            trim_dir = self.workdir + "/" + samp + "/trimming_results"
            map_dir = self.workdir + "/" + samp + "/mapping_results"
            yield MergeBAMfile(fastq1=trim_dir + "/" + samp + ".1.trimmed.fastq",
                               fastq2=trim_dir + "/" + samp + ".2.trimmed.fastq",
                               numCPUs=self.numCPUs,
                               indexfile=self.indexfile,
                               spliceFile=splice_file,
                               mappingLogFile=map_dir + "/mapping.log",
                               unalned=map_dir + "/unligned.fastq",
                               outsam=map_dir + "/" + samp + ".sam",
                               bam_file=map_dir + "/" + samp + ".bam",
                               sorted_bam_file=map_dir + "/" + samp + "_srt.bam",
                               map_ref=map_dir + "/" + samp + ".reflist",
                               map_dir=map_dir,
                               scriptdir=self.scriptdir,
                               ref_file=self.ref_file,
                               bindir=self.bindir,
                               workdir=self.workdir,
                               kingdom='prokarya')
            yield MergeBAMfile(fastq1=trim_dir + "/" + samp + ".1.trimmed.fastq",
                               fastq2=trim_dir + "/" + samp + ".2.trimmed.fastq",
                               numCPUs=self.numCPUs,
                               indexfile=self.indexfile,
                               spliceFile=splice_file,
                               mappingLogFile=map_dir + "/mapping.log",
                               unalned=map_dir + "/unligned.fastq",
                               outsam=map_dir + "/" + samp + ".sam",
                               bam_file=map_dir + "/" + samp + ".bam",
                               sorted_bam_file=map_dir + "/" + samp + "_srt.bam",
                               map_ref=map_dir + "/" + samp + ".reflist",
                               map_dir=map_dir,
                               scriptdir=self.scriptdir,
                               ref_file=self.ref_file,
                               bindir=self.bindir,
                               workdir=self.workdir,
                               kingdom='eukarya'
                               )


@inherits(SortBAMfileW)
class StringTieScoresBoth(luigi.WrapperTask):
    """From Mapping to Counting step for Eukaryotic and Prokaryotic."""

    euk_gff = Parameter()
    prok_gff = Parameter()

    def requires(self):
        """A pipeline that runs from mapping to count for euk and prok."""
        splice_list = [self.workdir + "/" +
                       f for f in os.listdir(self.workdir) if f.endswith('.splice')]
        if len(splice_list) > 1:
            splice_file = ','.join([splice_list])
        else:
            splice_file = splice_list[0]

        euk_gtf = self.workdir + "/" + \
            self.euk_gff.split("/")[-1].split(".gff")[0] + ".gtf"
        prok_gtf = self.workdir + "/" + \
            self.prok_gff.split("/")[-1].split(".gff")[0] + ".gtf"
        for samp, fastq in self.fastq_dic.iteritems():
            map_dir = self.workdir + "/" + samp + "/mapping_results"
            trim_dir = self.workdir + "/" + samp + "/trimming_results"
            yield StringTieScores(fastq1=trim_dir + "/" + samp + ".1.trimmed.fastq",
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
                                  in_bam_file=map_dir + "/" + "eukarya.bam",
                                  bindir=self.bindir)
            yield StringTieScores(fastq1=trim_dir + "/" + samp + ".1.trimmed.fastq",
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
                                  in_bam_file=map_dir + "/" + "prokarya.bam",
                                  bindir=self.bindir)
