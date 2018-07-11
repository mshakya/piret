#! /usr/bin/env python

"""Luigi Tasks to perform various RNA seq functions, from Mapping to Counting.

Mapping is done using hisat2 and counting is done using featurecounts and stringtie
"""

from __future__ import print_function
import os
import luigi
import sys
dir_path = os.path.dirname(os.path.realpath(__file__))
lib_path = os.path.abspath(os.path.join(dir_path, '..'))
sys.path.append(lib_path)
from luigi.contrib.external_program import ExternalProgramTask
from luigi import ExternalTask
from luigi import LocalTarget
from luigi import Parameter, IntParameter, DictParameter, ListParameter
from luigi.util import inherits, requires
import subprocess
from pypiret import FaQC
from plumbum.cmd import gffread, hisat2, python
from plumbum.cmd import samtools, stringtie, mv, awk
import pandas as pd
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
    num_cpus = IntParameter()

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
                "-q", "-p", self.num_cpus,
                self.fasta, self.hi_index]

    def program_environment(self):
        """Environmental variables for this program."""
        return {'PATH': os.environ["PATH"] + ":" + self.bindir}


class SAMindex(luigi.Task):
    """Create Hisat Indeices from given fasta file."""

    fasta = Parameter()
    workdir = Parameter()

    def output(self):
        """Require reference fasta format file."""
        if ',' in self.fasta:
            fas = self.fasta.split(",")
            return [LocalTarget(self.workdir + "/" + os.path.basename(fa) +
                                ".bedfile") for fa in fas]
        else:
            return [LocalTarget(self.workdir + "/" +
                                os.path.basename(self.fasta) +
                                ".bedfile")]

    def make_index(self, ref, workdir):
        """A function to make index from sam files."""
        index_options = ["faidx", ref]
        mv_options = [ref + ".fai", workdir]
        samtools_cmd = samtools[index_options]
        mv_cmd = mv[mv_options]
        samtools_cmd()
        mv_cmd()
        fa_name = os.path.basename(ref)
        return os.path.abspath(self.workdir + "/" + fa_name + ".fai")

    def create_bedfile(self, index_file):
        """Makes bed file from the sam index file."""
        out_bedfile = index_file.split(".fai")[0] + ".bedfile"
        bedfile_opt = ['BEGIN {FS="\t"}; {print $1 FS "0" FS $2}',
                       index_file]
        awk_cmd = ((awk[bedfile_opt]) > out_bedfile)
        awk_cmd()

    def requires(self):
        """Expected index output."""
        if ',' in self.fasta:
            fas = self.fasta.split(",")
            return [RefFile(fa) for fa in fas]
        else:
            return [RefFile(self.fasta)]

    def run(self):
        """Run hisat2-build command."""
        if ',' in self.fasta:
            fas = self.fasta.split(",")
            for fa in fas:
                ind_file = self.make_index(fa, self.workdir)
                self.create_bedfile(ind_file)


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

@requires(FaQC.PairedRunQC)
class Hisat(luigi.Task):
    """Mapping the QCed sequences to reference."""

    fastqs = ListParameter()
    indexfile = Parameter()
    spliceFile = Parameter()
    outsam = Parameter()
    ref_file = Parameter()
    map_dir = Parameter()

    def output(self):
        """SAM file output of the mapping."""
        bam_file = self.outsam.split(".sam")[0] + ".bam"
        return luigi.LocalTarget(bam_file)

    def run(self):
        """Run hisat2."""
        if self.spliceFile is "":
            hisat2_nosplice_option = ["-p", self.num_cpus,
                                      "-x", self.indexfile,
                                      "-1", self.fastqs[0],
                                      "-2", self.fastqs[1],
                                      "-S", self.outsam,
                                      "--no-unal",
                                      "--un-conc",
                                      os.path.join(self.map_dir,
                                                   "unaligned.fastq"),
                                      "2>", os.path.join(self.map_dir,
                                                         "mapping.log")]
            hisat2_cmd = hisat2[hisat2_nosplice_option]
            hisat2_cmd()
            self.sam2bam()
            self.sort_bam()
        else:
            h2_splice_option = ["--known-splicesite-infile", self.spliceFile,
                                "-p", self.num_cpus,
                                "-x", self.indexfile,
                                "-1", self.fastqs[0],
                                "-2", self.fastqs[1],
                                "-S", self.outsam,
                                "--no-unal",
                                "--un-conc",
                                os.path.join(self.map_dir,
                                             "unaligned.fastq"),
                                "2>", os.path.join(self.map_dir,
                                                   "mapping.log")]
            hisat2_cmd = hisat2[h2_splice_option]
            hisat2_cmd()
            self.sam2bam()
            self.sort_bam()

    def sam2bam(self):
        """Convert SAM to BAM file."""
        bam_file = self.outsam.split(".sam")[0] + ".bam"
        options = ["view", "-bS", "-F",
                   "4", self.outsam, "-o", bam_file]
        samtools_cmd = samtools[options]
        samtools_cmd()

    def sort_bam(self):
        """Sort BAM file."""
        bam_file = self.outsam.split(".sam")[0] + ".bam"
        sorted_bam_file = bam_file.split(".bam")[0] + "_srt.bam"
        options = ["sort", bam_file,
                   "-o", sorted_bam_file]
        samtools_cmd = samtools[options]
        samtools_cmd()


class HisatMapW(luigi.WrapperTask):
    """A wrapper task for mapping."""

    fastq_dic = DictParameter()
    ref_file = luigi.Parameter()
    indexfile = luigi.Parameter()
    bindir = luigi.Parameter()
    workdir = luigi.Parameter()
    num_cpus = luigi.IntParameter()

    def requires(self):
        """A wrapper task for running mapping."""
        splist = [self.workdir + "/" +
                  f for f in os.listdir(self.workdir) if f.endswith('.splice')]
        if len(splist) > 1:
            splice_file = ','.join(splist)
        elif len(splist) == 1:
            splice_file = splist[0]
        else:
            splice_file = ''
        for samp, fastq in self.fastq_dic.items():
            trim_dir = self.workdir + "/" + samp + "/trimming_results"
            map_dir = self.workdir + "/" + samp + "/mapping_results"
            if os.path.isdir(map_dir) is False:
                os.makedirs(map_dir)
            yield Hisat(fastqs=[trim_dir + "/" + samp + ".1.trimmed.fastq",
                                trim_dir + "/" + samp + ".2.trimmed.fastq"],
                        qc_outdir=trim_dir,
                        num_cpus=self.num_cpus,
                        indexfile=self.indexfile,
                        spliceFile=splice_file,
                        outsam=map_dir + "/" + samp + ".sam",
                        ref_file=self.ref_file,
                        bindir=self.bindir,
                        map_dir=map_dir,
                        sample=samp)


@requires(Hisat)
class RefNames(luigi.Task):
    """Extract the name of chromosomes where reads were mapped in BAM file."""

    map_ref = Parameter()

    def output(self):
        """Split Files."""
        return luigi.LocalTarget(self.map_ref)

    def run(self):
        """Extract names of chromosome from SAM file."""
        bam_file = self.outsam.split(".sam")[0] + ".bam"
        sorted_bam_file = bam_file.split(".bam")[0] + "_srt.bam"
        grep_cmd = "samtools view -h %s |\
                    grep -v -e '^.SQ' |\
                    grep -v -e 'PN:' |\
                    grep -v -e 'SO:' |\
                    awk -F '\t' '{print $3}' |\
                    sort | uniq |\
                    grep -v '*' > %s" % (sorted_bam_file, self.map_ref)
        subprocess.Popen(grep_cmd, shell=True)


@inherits(HisatMapW)
class GetRefNames(luigi.WrapperTask):
    """From Mapping to Counting step for Eukaryotic and Prokaryotic."""

    def requires(self):
        """A pipeline that runs from mapping to count for euk and prok."""
        splist = [self.workdir + "/" +
                  f for f in os.listdir(self.workdir) if f.endswith('.splice')]
        if len(splist) > 1:
            splice_file = ','.join(splist)
        else:
            splice_file = splist[0]
        for samp, fastq in self.fastq_dic.items():
            trim_dir = self.workdir + "/" + samp + "/trimming_results"
            map_dir = self.workdir + "/" + samp + "/mapping_results"
            yield RefNames(fastqs=[trim_dir + "/" + samp + ".1.trimmed.fastq",
                                   trim_dir + "/" + samp + ".2.trimmed.fastq"],
                           num_cpus=self.num_cpus,
                           indexfile=self.indexfile,
                           spliceFile=splice_file,
                           outsam=map_dir + "/" + samp + ".sam",
                           map_ref=map_dir + "/" + samp + ".reflist",
                           ref_file=self.ref_file,
                           bindir=self.bindir,
                           sample=samp,
                           qc_outdir=trim_dir,
                           map_dir=map_dir)


@inherits(GFF2GTF)
@requires(Hisat)
class StringTieScores(luigi.Task):
    """Calculate scores using string tie."""

    out_gtf = luigi.Parameter()
    out_cover = luigi.Parameter()
    out_abun = luigi.Parameter()
    in_bam_file = luigi.Parameter()

    def output(self):
        """Index output."""
        return LocalTarget(self.out_gtf)

    def run(self):
        """Run stringtie."""
        stringtie_opt = ["-o", self.out_gtf,
                         "-G", self.gff_file,
                         "-C", self.out_cover,
                         "-A", self.out_abun,
                         self.in_bam_file]
        stringtie_cmd = stringtie[stringtie_opt]
        stringtie_cmd()


# @inherits(GFF2GTF)
@inherits(HisatMapW)
class StringTieScoresW(luigi.WrapperTask):
    """From Mapping to Counting step for Eukaryotic reference."""

    gff_file = Parameter()
    kingdom = Parameter()

    def requires(self):
        """A wrapper for running Stringtie scores on all samples."""
        splst = [self.workdir + "/" +
                 f for f in os.listdir(self.workdir) if f.endswith('.splice')]
        if len(splst) > 1:
            splice_file = ','.join(splst)
        elif len(splst) == 1:
            splice_file = splst[0]
        else:
            splice_file = ""
        for samp, fastq in self.fastq_dic.items():
            map_dir = self.workdir + "/" + samp + "/mapping_results"
            trim_dir = self.workdir + "/" + samp + "/trimming_results"
            stng_dir = self.workdir + "/" + samp + "/stie_results"
            if os.path.isdir(stng_dir) is False:
                os.makedirs(stng_dir)
            if self.kingdom in ['prokarya', 'eukarya']:
                if self.kingdom == 'prokarya':
                    apd = '_prok'
                elif self.kingdom == 'eukarya':
                    apd = '_euk'
                gff_name = self.gff_file.split(".gff")[0].split("/")[-1]

                yield StringTieScores(fastqs=[trim_dir + "/" + samp + ".1.trimmed.fastq",
                                              trim_dir + "/" + samp + ".2.trimmed.fastq"],
                                      num_cpus=self.num_cpus,
                                      indexfile=self.indexfile,
                                      spliceFile=splice_file,
                                      outsam=map_dir + "/" + samp + ".sam",
                                      ref_file=self.ref_file,
                                      gff_file=self.gff_file,
                                      out_gtf=stng_dir + "/" + samp + "_" + gff_name + apd + "_sTie.gtf",
                                      out_cover=stng_dir + "/" + samp + "_" + gff_name + apd + "_covered_sTie.gtf",
                                      out_abun=stng_dir + "/" + samp + "_" + gff_name + apd + "_sTie.tab",
                                      in_bam_file=map_dir + "/" + samp + "_srt.bam",
                                      bindir=self.bindir,
                                      workdir=self.workdir,
                                      sample=samp,
                                      qc_outdir=trim_dir,
                                      map_dir=map_dir)
            elif self.kingdom == 'both':
                prok_gff = os.path.basename(self.gff_file.split(";")[0]).split(".gff")[0]
                euk_gff = os.path.basename(self.gff_file.split(";")[1]).split(".gff")[0]
                yield StringTieScores(fastqs=[trim_dir + "/" + samp +
                                              ".1.trimmed.fastq",
                                              trim_dir + "/" + samp +
                                              ".2.trimmed.fastq"],
                                      num_cpus=self.num_cpus,
                                      indexfile=self.indexfile,
                                      spliceFile=splice_file,
                                      outsam=map_dir + "/" + samp + ".sam",
                                      ref_file=self.ref_file,
                                      out_gtf=stng_dir + "/" + samp + "_" + prok_gff + "_prok" + "_sTie.gtf",
                                      out_cover=stng_dir + "/" + samp + "_" + prok_gff + "_prok" + "_covered_sTie.gtf",
                                      out_abun=stng_dir + "/" + samp + "_" + prok_gff + "_prok" + "_sTie.tab",
                                      in_bam_file=map_dir + "/" + samp + "_srt_prok.bam",
                                      bindir=self.bindir,
                                      workdir=self.workdir,
                                      gff_file=self.gff_file.split(";")[0],
                                      sample=samp,
                                      qc_outdir=trim_dir,
                                      map_dir=map_dir)
                yield StringTieScores(fastqs=[trim_dir + "/" + samp + ".1.trimmed.fastq",
                                      trim_dir + "/" + samp + ".2.trimmed.fastq"],
                                      num_cpus=self.num_cpus,
                                      indexfile=self.indexfile,
                                      spliceFile=splice_file,
                                      outsam=map_dir + "/" + samp + ".sam",
                                      ref_file=self.ref_file,
                                      out_gtf=stng_dir + "/" + samp + "_" + euk_gff + "_euk" + "_sTie.gtf",
                                      out_cover=stng_dir + "/" + samp + "_" + euk_gff + "_euk" + "_covered_sTie.gtf",
                                      out_abun=stng_dir + "/" + samp + "_" + euk_gff + "_euk" + "_sTie.tab",
                                      in_bam_file=map_dir + "/" + samp + "_srt_prok.bam",
                                      bindir=self.bindir,
                                      workdir=self.workdir,
                                      gff_file=self.gff_file.split(";")[0],
                                      sample=samp,
                                      qc_outdir=trim_dir,
                                      map_dir=map_dir)

@requires(Hisat)
class Split2ProkEuk(luigi.Task):
    """ Split BAM file to prok and euk"""

    ref_file = Parameter()
    workdir = Parameter()

    def output(self):
        """Split Files."""
        bam_file = self.outsam.split(".sam")[0] + "_srt.bam"
        prok_out = bam_file.split(".bam")[0] + "_prok.bam"
        euk_out = bam_file.split(".bam")[0] + "_euk.bam"
        return [luigi.LocalTarget(prok_out), luigi.LocalTarget(euk_out)]

    def split_aln_file(self, workdir, fasta, bam_file, out_bamfile):
        """Split"""
        bed_file = workdir + "/" + os.path.basename(fasta) + ".bedfile"
        samtools_opt = ["view", "-b", "-L", bed_file, bam_file]
        samtools_cmd = (samtools[samtools_opt]) > out_bamfile
        samtools_cmd()

    def run(self):
        """ Split BAM file to proks and euks"""
        fastas = self.ref_file.split(",")
        bam_file = self.outsam.split(".sam")[0] + "_srt.bam"
        prok_out = bam_file.split(".bam")[0] + "_prok.bam"
        euk_out = bam_file.split(".bam")[0] + "_euk.bam"
        self.split_aln_file(self.workdir, fastas[0], bam_file, prok_out)
        self.split_aln_file(self.workdir, fastas[1], bam_file, euk_out)


@inherits(HisatMapW)
class Split2ProkEukW(luigi.WrapperTask):
    """Group chromosomes to prokaryotic and eukaryotic."""
    ref_file = Parameter()
    workdir = Parameter()

    def requires(self):
        """A pipeline that runs from mapping to count for euk and prok."""
        splst = [self.workdir + "/" +
                 f for f in os.listdir(self.workdir) if f.endswith('.splice')]
        if len(splst) > 1:
            splice_file = ','.join(splst)
        else:
            splice_file = splst[0]
        for samp, fastq in self.fastq_dic.items():
            trim_dir = self.workdir + "/" + samp + "/trimming_results"
            map_dir = self.workdir + "/" + samp + "/mapping_results"
            yield Split2ProkEuk(fastqs=[trim_dir + "/" + samp +
                                        ".1.trimmed.fastq",
                                        trim_dir + "/" + samp +
                                        ".2.trimmed.fastq"],
                                num_cpus=self.num_cpus,
                                indexfile=self.indexfile,
                                spliceFile=splice_file,
                                outsam=map_dir + "/" + samp + ".sam",
                                map_dir=map_dir,
                                ref_file=self.ref_file,
                                bindir=self.bindir,
                                workdir=self.workdir,
                                sample=samp,
                                qc_outdir=trim_dir)


# @inherits(HisatMapW)
# class StringTieScoresBoth(luigi.WrapperTask):
#     """From Mapping to Counting step for Eukaryotic and Prokaryotic."""

#     euk_gff = Parameter()
#     prok_gff = Parameter()

#     def requires(self):
#         """A pipeline that runs from mapping to count for euk and prok."""
#         splice_list = [self.workdir + "/" +
#                        f for f in os.listdir(self.workdir) if f.endswith('.splice')]
#         if len(splice_list) > 1:
#             splice_file = ','.join(splice_list)
#         else:
#             splice_file = splice_list[0]

#         euk_gtf = self.workdir + "/" + \
#             self.euk_gff.split("/")[-1].split(".gff")[0] + ".gtf"
#         prok_gtf = self.workdir + "/" + \
#             self.prok_gff.split("/")[-1].split(".gff")[0] + ".gtf"
#         for samp, fastq in self.fastq_dic.items():
#             map_dir = self.workdir + "/" + samp + "/mapping_results"
#             trim_dir = self.workdir + "/" + samp + "/trimming_results"
#             yield StringTieScores(fastq1=trim_dir + "/" + samp + ".1.trimmed.fastq",
#                                   fastq2=trim_dir + "/" + samp + ".2.trimmed.fastq",
#                                   num_cpus=self.num_cpus,
#                                   indexfile=self.indexfile,
#                                   spliceFile=splice_file,
#                                   # mappingLogFile=map_dir + "/mapping.log",
#                                   outsam=map_dir + "/" + samp + ".sam",
#                                   bam_file=map_dir + "/" + samp + ".bam",
#                                   sorted_bam_file=map_dir + "/" + samp + "_srt.bam",
#                                   ref_file=self.ref_file,
#                                   gtf=euk_gtf,
#                                   out_gtf=map_dir + "/" + samp + "_euk_sTie.gtf",
#                                   out_cover=map_dir + "/" + samp + "_euk_covered_sTie.gtf",
#                                   out_abun=map_dir + "/" + samp + "_euk_sTie.tab",
#                                   in_bam_file=map_dir + "/" + "eukarya.bam",
#                                   bindir=self.bindir)
#             yield StringTieScores(fastq1=trim_dir + "/" + samp + ".1.trimmed.fastq",
#                                   fastq2=trim_dir + "/" + samp + ".2.trimmed.fastq",
#                                   num_cpus=self.num_cpus,
#                                   indexfile=self.indexfile,
#                                   spliceFile=splice_file,
#                                   # mappingLogFile=map_dir + "/mapping.log",
#                                   outsam=map_dir + "/" + samp + ".sam",
#                                   bam_file=map_dir + "/" + samp + ".bam",
#                                   sorted_bam_file=map_dir + "/" + samp + "_srt.bam",
#                                   ref_file=self.ref_file,
#                                   gtf=prok_gtf,
#                                   out_gtf=map_dir + "/" + samp + "_prok_sTie.gtf",
#                                   out_cover=map_dir + "/" + samp + "_prok_covered_sTie.gtf",
#                                   out_abun=map_dir + "/" + samp + "_prok_sTie.tab",
#                                   in_bam_file=map_dir + "/" + "prokarya.bam",
#                                   bindir=self.bindir)
