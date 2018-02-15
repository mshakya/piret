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
from plumbum.cmd import touch, bamtools, hisat2, gffread, python, cp, rm, cut, samtools, stringtie
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


class GFF2GTF(luigi.Task):
    """Converts GFF to GTF format."""

    gff_file = Parameter()
    workdir = Parameter()
    bindir = Parameter()

    def requires(self):
        """Require reference gff(s) file."""
        if ',' in self.gff_file:
            gffs = self.gff_file.split(",")
            return [RefFile(os.path.abspath(gff)) for gff in gffs]
        else:
            return [RefFile(os.path.abspath(self.gff_file))]

    def output(self):
        """Output converted gtf file."""
        if ',' in self.gff_file:
            gffs = self.gff_file.split(",")
            return [LocalTarget(self.workdir + "/" +
                                os.path.abspath(gff).split("/")[-1].rsplit(".", 1)[0] +
                                ".gtf") for gff in gffs]
        else:
            gff = os.path.abspath(self.gff_file)
            return LocalTarget(self.workdir + "/" +
                               gff.split("/")[-1].rsplit(".", 1)[0] +
                               ".gtf")

    def run(self):
        """Main."""
        if ',' in self.gff_file:
            gffs = self.gff_file.split(",")
            for gff in gffs:
                gff_abs = os.path.abspath(gff)
                out_file = self.workdir + "/" +\
                    gff_abs.split("/")[-1].rsplit(".", 1)[0] + ".gtf"
                gffread_option = [gff_abs, "-O", "-o", out_file]
                gffread_cmd = gffread[gffread_option]
                gffread_cmd()
        else:
            gff = os.path.abspath(self.gff_file)
            out_file = self.workdir + "/" +\
                gff.split("/")[-1].rsplit(".", 1)[0] + ".gtf"
            gffread_option = [self.gff_file, "-O","-o", out_file]
            gffread_cmd = gffread[gffread_option]
            gffread_cmd()


@requires(GFF2GTF)
class CreateSplice(ExternalProgramTask):
    """Find splice sites off gtf file."""

    def output(self):
        """Splice site output."""
        if ',' in self.gff_file:
            gffs = self.gff_file.split(",")
            return [LocalTarget(self.workdir + "/" +
                                os.path.abspath(gff).split("/")[-1].split(".")[0] +
                                ".splice") for gff in gffs]
        else:
            gff = os.path.abspath(self.gff_file)
            return LocalTarget(self.workdir + "/" +
                               gff.split("/")[-1].split(".")[0] +
                               ".splice")

    def run(self):
        """Main."""
        if ',' in self.gff_file:
            gffs = self.gff_file.split(",")
            for gff in gffs:
                gtf_file = self.workdir + "/" +\
                    gff.split("/")[-1].rsplit(".", 1)[0] + ".gtf"
                out_file = self.workdir + "/" +\
                    gff.split("/")[-1].rsplit(".", 1)[0] + ".splice"
                hess_fpath = self.bindir + "/../scripts/hisat2_extract_splice_sites.py"
                hess_opt = [hess_fpath, "-i", gtf_file, "-o", out_file]
                hess_cmd = python[hess_opt]
                hess_cmd()
        else:
            gff = os.path.abspath(self.gff_file)
            gtf_file = self.workdir + "/" +\
                gff.split("/")[-1].rsplit(".", 1)[0] + ".gtf"
            out_file = self.workdir + "/" +\
                self.gff_file.split("/")[-1].rsplit(".", 1)[0] + ".splice"
            hess_fpath = self.bindir + "/../scripts/hisat2_extract_splice_sites.py"
            hess_opt = [hess_fpath, "-i", gtf_file, "-o", out_file]
            hess_cmd = python[hess_opt]
            hess_cmd()


# @requires(GFF2GTF)
# class GetIntGenRegions(ExternalProgramTask):
#     """Find splice sites off gtf file."""

#     def output(self):
#         """Splice site output."""
#         if ',' in self.gff_file:
#             gffs = self.gff_file.split(",")
#             return [LocalTarget(self.workdir + "/" +
#                                 gff.split("/")[-1].split(".")[0] +
#                                 ".ints") for gff in gffs]
#         else:
#             return LocalTarget(self.workdir + "/" +
#                                self.gff_file.split("/")[-1].split(".")[0] +
#                                ".ints")

#     def run(self):
#         """Main."""
#         if ',' in self.ref_file:
#             refs = self.refs.split(",")
#             for ref in refs:
#                 target_fa = self.workdir + "/" + ref
#                 cp_opt = [ref, target_fa]
#                 cp_cmd = cp[cp_opt]
#                 cp_cmd()
#                 sam_opt = ["faidx", target_fa]
#                 sam_cmd = samtools[sam_opt]
#                 sam_cmd()
#                 rm_opt = [target_fa]
#                 rm_cmd = rm[rm_opt]
#                 rm_cmd()
#                 index_file = target_fa + ".fai"
#                 cut_opt = ["-f1,2", index_file]
#                 cut_cmd = cut[cut_opt]
#                 cut_cmd()



#         else:
#             gtf_file = self.workdir + "/" +\
#                 self.gff_file.split("/")[-1].split(".")[0] + ".gtf"
#             hess_fpath = self.bindir + "/../scripts/hisat2_extract_splice_sites.py"
#             hess_opt = [hess_fpath, "-i", gtf_file, "-o", out_file]
#             hess_cmd = python[hess_opt]


# @inherits(FastQC.PairedRunQC)
class Hisat(luigi.Task):
    """Mapping the QCed sequences to reference."""

    fastq1 = Parameter()
    fastq2 = Parameter()
    numCPUs = luigi.IntParameter()
    indexfile = Parameter()
    spliceFile = Parameter()
    mappingLogFile = Parameter()
    outsam = Parameter()
    unalned = Parameter()
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
                                  bindir=self.bindir),
               HisatIndex(fasta=self.ref_file,
                          hi_index=self.indexfile,
                          bindir=self.bindir,
                          numCPUs=self.numCPUs)]

    def output(self):
        """SAM file output of the mapping."""
        return LocalTarget(self.outsam)

    def run(self):
        """Run hisat2."""
        if self.spliceFile is "":
            hisat2_nosplice_option = ["-p", self.numCPUs,
                                      "-x", self.indexfile,
                                      "-1", self.fastq1,
                                      "-2", self.fastq2,
                                      "-S", self.outsam,
                                      "--un-conc", self.unalned,
                                      "2>", self.mappingLogFile]
            hisat2_cmd = hisat2[hisat2_nosplice_option]
            hisat2_cmd()
        else:
            hisat2_splice_option = ["--known-splicesite-infile", self.spliceFile,
                                    "-p", self.numCPUs,
                                    "-x", self.indexfile,
                                    "-1", self.fastq1,
                                    "-2", self.fastq2,
                                    "-S", self.outsam,
                                    "--un-conc", self.unalned,
                                    "2>", self.mappingLogFile]
            hisat2_cmd = hisat2[hisat2_splice_option]
            hisat2_cmd()


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
            yield SAM2BAMfile(fastq1=trim_dir + "/" + samp + ".1.trimmed.fastq",
                              fastq2=trim_dir + "/" + samp + ".2.trimmed.fastq",
                              numCPUs=self.numCPUs,
                              indexfile=self.indexfile,
                              spliceFile=splice_file,
                              mappingLogFile=map_dir + "/mapping.log",
                              unalned=map_dir + "/unligned.fastq",
                              outsam=map_dir + "/" + samp + ".sam",
                              bam_file=map_dir + "/" + samp + ".bam",
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
        """A wrapper task for converting sam to bam."""
        splice_list = [self.workdir + "/" +
                       f for f in os.listdir(self.workdir) if f.endswith('.splice')]
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
            splice_file = ','.join(splice_list)
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
                           ref_file=self.ref_file,
                           bindir=self.bindir)


@inherits(GFF2GTF)
@inherits(SortBAMfile)
class StringTieScores(luigi.Task):
    """Calculate scores using string tie."""

    in_gtf = luigi.Parameter()
    out_gtf = luigi.Parameter()
    out_cover = luigi.Parameter()
    out_abun = luigi.Parameter()
    in_bam_file = luigi.Parameter()

    def requires(self):
        """Require reference fasta format file."""
        if ',' in self.gff_file:
            gffs = self.gff_file.split(",")
            gff_depen = [GFF2GTF(gff_file=self.gff_file,
                                 workdir=self.workdir,
                                 bindir=self.bindir) for gff in gffs]
        else:
            gff_depen = [GFF2GTF(gff_file=self.gff_file,
                                 workdir=self.workdir,
                                 bindir=self.bindir)]
        return [
            SortBAMfile(fastq1=self.fastq1,
                        fastq2=self.fastq2,
                        numCPUs=self.numCPUs,
                        indexfile=self.indexfile,
                        spliceFile=self.spliceFile,
                        mappingLogFile=self.mappingLogFile,
                        unalned=self.unalned,
                        outsam=self.outsam,
                        bam_file=self.bam_file,
                        sorted_bam_file=self.sorted_bam_file,
                        ref_file=self.ref_file,
                        bindir=self.bindir)] + gff_depen

    def output(self):
        """Index output."""
        return LocalTarget(self.out_gtf)

    def run(self):
        """Run stringtie."""
        stringtie_opt = ["-o", self.out_gtf,
                         "-G", self.in_gtf,
                         "-C", self.out_cover,
                         "-A", self.out_abun,
                         self.in_bam_file]
        stringtie_cmd = stringtie[stringtie_opt]
        stringtie_cmd()


@inherits(GFF2GTF)
@inherits(SortBAMfileW)
class StringTieScoresW(luigi.WrapperTask):
    """From Mapping to Counting step for Eukaryotic reference."""

    gff_file = Parameter()
    kingdom = Parameter()

    def requires(self):
        """A wrapper for running Stringtie scores on all samples."""
        splice_list = [self.workdir + "/" +
                       f for f in os.listdir(self.workdir) if f.endswith('.splice')]
        if len(splice_list) > 1:
            splice_file = ','.join(splice_list)
        elif len(splice_list) == 1:
            splice_file = splice_list[0]
        for samp, fastq in self.fastq_dic.iteritems():
            map_dir = self.workdir + "/" + samp + "/mapping_results"
            trim_dir = self.workdir + "/" + samp + "/trimming_results"
            if os.path.isdir(map_dir) is False:
                os.makedirs(map_dir)
            if self.kingdom in ['prokarya', 'eukarya']:
                if self.kingdom == 'prokarya':
                    append_name = '_prok'
                elif self.kingdom == 'eukarya':
                    append_name = '_euk'
                if ',' in self.gff_file:
                    gff_list = [os.path.abspath(gff) for gff in self.gff_file.split(",")]
                    for gff in gff_list:
                        gtf = self.workdir + "/" + gff.split("/")[-1].split(".gff")[0] + ".gtf"
                        gff_name = gtf.split(".gtf")[0].split("/")[-1]
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
                                          ref_file=self.ref_file,
                                          in_gtf=gtf,
                                          gff_file=self.gff_file,
                                          out_gtf=map_dir + "/" + samp + "_" + gff_name + append_name + "_sTie.gtf",
                                          out_cover=map_dir + "/" + samp + "_" + gff_name + append_name + "_covered_sTie.gtf",
                                          out_abun=map_dir + "/" + samp + "_" + gff_name + append_name + "_sTie.tab",
                                          in_bam_file=map_dir + "/" + samp + "_srt.bam",
                                          bindir=self.bindir,
                                          workdir=self.workdir)

                else:
                    gtf = self.workdir + "/" + self.gff_file.split("/")[-1].split(".gff")[0] + ".gtf"
                    gff_name = gtf.split(".gtf")[0].split("/")[-1]
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
                                          ref_file=self.ref_file,
                                          in_gtf=gtf,
                                          gff_file=self.gff_file,
                                          out_gtf=map_dir + "/" + samp + "_" + gff_name + append_name + "_sTie.gtf",
                                          out_cover=map_dir + "/" + samp + "_" + gff_name + append_name + "_covered_sTie.gtf",
                                          out_abun=map_dir + "/" + samp + "_" + gff_name + append_name + "_sTie.tab",
                                          in_bam_file=map_dir + "/" + samp + "_srt.bam",
                                          bindir=self.bindir,
                                          workdir=self.workdir)
            elif self.kingdom == 'both':
                prok_gtf = self.workdir + "/" + \
                    self.gff_file.split(";")[0].split("/")[-1].split(".gff")[0] + ".gtf"
                euk_gtf = self.workdir + "/" + \
                    self.gff_file.split(";")[1].split("/")[-1].split(".gff")[0] + ".gtf"
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
                                      ref_file=self.ref_file,
                                      gtf=prok_gtf,
                                      out_gtf=map_dir + "/" + samp + "_prok_sTie.gtf",
                                      out_cover=map_dir + "/" + samp + "_prok_covered_sTie.gtf",
                                      out_abun=map_dir + "/" + samp + "_prok_sTie.tab",
                                      in_bam_file=map_dir + "/prokarya.bam",
                                      bindir=self.bindir,
                                      workdir=self.workdir,
                                      gff_file=self.gff_file)
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
                                      ref_file=self.ref_file,
                                      gtf=euk_gtf,
                                      out_gtf=map_dir + "/" + samp + "_euk_sTie.gtf",
                                      out_cover=map_dir + "/" + samp + "_euk_covered_sTie.gtf",
                                      out_abun=map_dir + "/" + samp + "_euk_sTie.tab",
                                      in_bam_file=map_dir + "/eukarya.bam",
                                      bindir=self.bindir,
                                      workdir=self.workdir,
                                      gff_file=self.gff_file)


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
        bam_cmd_opt = ["split", "-reference", "-in", self.sorted_bam_file]
        bam_cmd = bamtools[bam_cmd_opt]
        bam_cmd()
        touch[self.split_logfile]()


@inherits(GetRefNames)
class SplitBAMBoth(luigi.WrapperTask):
    """From Mapping to Counting step for Eukaryotic and Prokaryotic."""

    def requires(self):
        """A pipeline that runs from mapping to count for euk and prok."""
        splice_list = [self.workdir + "/" +
                       f for f in os.listdir(self.workdir) if f.endswith('.splice')]
        if len(splice_list) > 1:
            splice_file = ','.join(splice_list)
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
            splice_file = ','.join(splice_list)
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
            splice_file = ','.join(splice_list)
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
            splice_file = ','.join(splice_list)
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
                                  ref_file=self.ref_file,
                                  gtf=prok_gtf,
                                  out_gtf=map_dir + "/" + samp + "_prok_sTie.gtf",
                                  out_cover=map_dir + "/" + samp + "_prok_covered_sTie.gtf",
                                  out_abun=map_dir + "/" + samp + "_prok_sTie.tab",
                                  in_bam_file=map_dir + "/" + "prokarya.bam",
                                  bindir=self.bindir)
