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
import subprocess
# from Bio import SeqIO


class RefFile(ExternalTask):
    """An ExternalTask like this."""

    path = Parameter()

    def output(self):
        """Check."""
        return LocalTarget(os.path.abspath(self.path))


class HisatIndex(ExternalProgramTask):
    """Create Hisat Indeices from given fasta file."""

    fasta = Parameter()
    hi_index = Parameter()
    bindir = Parameter()

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
                "-q", self.fasta, self.hi_index]

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


class CreateSplice(ExternalProgramTask):
    """Find splice sites off gtf file."""

    gff_file = Parameter()
    bindir = Parameter()
    workdir = Parameter()

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
            return ["scripts/hisat2_extract_splice_sites.py",
                    "-i", gtf_file, "-o", out_file]
        else:
            gtf_file = self.workdir + "/" +\
                self.gff_file.split("/")[-1].split(".")[0] + ".gtf"

            return ["scripts/hisat2_extract_splice_sites.py",
                    "-i", gtf_file, "-o", self.workdir + "/" +
                    self.gff_file.split("/")[-1].split(".")[0] + ".splice"]


class MapHisat(ExternalProgramTask):
    """Mapping the QCed sequences to reference."""

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
        """SAM file output of the mapping."""
        return LocalTarget(self.outsam)

    def program_args(self):
        """Run hisat2."""
        if self.spliceFile is None:
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


class SAM2BAMfile(ExternalProgramTask):
    """Convert SAM file to BAM file with only mapped reads."""

    fastq1 = Parameter()
    fastq2 = Parameter()
    numCPUs = luigi.IntParameter()
    indexfile = Parameter()
    spliceFile = Parameter()
    mappingLogFile = Parameter()
    outsam = Parameter()
    unalned = Parameter()
    bam_file = Parameter()

    def requires(self):
        """Required Mapping step to be complete."""
        return [MapHisat(fastq1=self.fastq1,
                         fastq2=self.fastq2,
                         numCPUs=self.numCPUs,
                         indexfile=self.indexfile,
                         spliceFile=self.spliceFile,
                         mappingLogFile=self.mappingLogFile,
                         outsam=self.outsam,
                         unalned=self.unalned)]

    def output(self):
        """Output BAM file."""
        return LocalTarget(self.bam_file)

    def program_args(self):
        """Output BAM file with only mapped reads."""
        return ["samtools", "view", "-bSh", "-F",
                "4", self.outsam, "-o", self.bam_file]


class SortBAMfile(ExternalProgramTask):
    """Sort BAM file."""

    fastq1 = Parameter()
    fastq2 = Parameter()
    numCPUs = luigi.IntParameter()
    indexfile = Parameter()
    spliceFile = Parameter()
    mappingLogFile = Parameter()
    outsam = Parameter()
    unalned = Parameter()
    bam_file = Parameter()
    sorted_bam_file = Parameter()

    def requires(self):
        """Require SAM2BAMfile to finish."""
        return [SAM2BAMfile(fastq1=self.fastq1,
                            fastq2=self.fastq2,
                            numCPUs=self.numCPUs,
                            indexfile=self.indexfile,
                            spliceFile=self.spliceFile,
                            mappingLogFile=self.mappingLogFile,
                            outsam=self.outsam,
                            bam_file=self.bam_file,
                            unalned=self.unalned)]

    def output(self):
        """Sorted BAM file."""
        return LocalTarget(self.sorted_bam_file)

    def program_args(self):
        """Sort BAM file."""
        return ["samtools", "sort",
                self.bam_file, "-o", self.sorted_bam_file]


class RefNames(luigi.Task):
    """Extract the name of chromosomes where reads mapped from BAM file."""

    fastq1 = Parameter()
    fastq2 = Parameter()
    numCPUs = luigi.IntParameter()
    indexfile = Parameter()
    spliceFile = Parameter()
    mappingLogFile = Parameter()
    outsam = Parameter()
    unalned = Parameter()
    bam_file = Parameter()
    workdir = Parameter()
    sorted_bam_file = Parameter()
    map_ref = Parameter()

    def requires(self):
        """Require sort step to be finished."""
        return [SortBAMfile(fastq1=self.fastq1,
                            fastq2=self.fastq2,
                            numCPUs=self.numCPUs,
                            indexfile=self.indexfile,
                            spliceFile=self.spliceFile,
                            mappingLogFile=self.mappingLogFile,
                            outsam=self.outsam,
                            bam_file=self.bam_file,
                            unalned=self.unalned,
                            sorted_bam_file=self.sorted_bam_file)]

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


class SplitBAMfile(luigi.Task):
    """Split BAM file to individual chromosomes."""

    fastq1 = Parameter()
    fastq2 = Parameter()
    numCPUs = luigi.IntParameter()
    indexfile = Parameter()
    spliceFile = Parameter()
    mappingLogFile = Parameter()
    outsam = Parameter()
    unalned = Parameter()
    bam_file = Parameter()
    workdir = Parameter()
    sorted_bam_file = Parameter()
    map_ref = Parameter()
    split_logfile = Parameter()

    def get_file_path(self):
        """Get bam filename."""
        with open(self.map_ref, 'r') as mrf:
            bam = self.sorted_bam_file.split(
                ".")[0] + ".REF_" + str(mrf.readline()) + ".bam"
        return bam

    def requires(self):
        """Require sort step to be finished."""
        return [RefNames(fastq1=self.fastq1,
                         fastq2=self.fastq2,
                         numCPUs=self.numCPUs,
                         indexfile=self.indexfile,
                         spliceFile=self.spliceFile,
                         mappingLogFile=self.mappingLogFile,
                         outsam=self.outsam,
                         bam_file=self.bam_file,
                         unalned=self.unalned,
                         sorted_bam_file=self.sorted_bam_file,
                         workdir=self.workdir,
                         map_ref=self.map_ref)]

    def output(self):
        """Split Files."""
        return luigi.LocalTarget(self.split_logfile)

    def run(self):
        """Split BAM file based on chromosome."""
        bam_cmd = "bamtools split -reference -in %s 2> %s" % (
            self.sorted_bam_file, self.split_logfile)
        subprocess.Popen(bam_cmd, shell=True)


class SplitRefProkEuk(luigi.Task):
    """Create a input BAM file list of prok and euk."""

    fastq1 = Parameter()
    fastq2 = Parameter()
    numCPUs = luigi.IntParameter()
    indexfile = Parameter()
    spliceFile = Parameter()
    mappingLogFile = Parameter()
    outsam = Parameter()
    unalned = Parameter()
    bam_file = Parameter()
    workdir = Parameter()
    sorted_bam_file = Parameter()
    map_ref = Parameter()
    map_dir = Parameter()

    def requires(self):
        """Require generating list of references to be finished."""
        return [RefNames(fastq1=self.fastq1,
                         fastq2=self.fastq2,
                         numCPUs=self.numCPUs,
                         indexfile=self.indexfile,
                         spliceFile=self.spliceFile,
                         mappingLogFile=self.mappingLogFile,
                         outsam=self.outsam,
                         bam_file=self.bam_file,
                         unalned=self.unalned,
                         sorted_bam_file=self.sorted_bam_file,
                         workdir=self.workdir,
                         map_ref=self.map_ref)]

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
        print(prok_list)
        euk_chroms_fp = self.map_dir + "/" + "eukarya_chromos.fullpath"
        prok_chroms_fp = self.map_dir + "/" + "prokarya_chromos.fullpath"
        with open(euk_chroms_fp, 'w') as ecf:
            for f in euk_list:
                ecf.write("%s\n" % f)
        with open(prok_chroms_fp, 'w') as pcf:
            for f in prok_list:
                pcf.write("%s\n" % f)


class MergeBAMfile(ExternalProgramTask):
    """Merge BAM file to prokaryotic and Eukaryotic."""

    fastq1 = Parameter()
    fastq2 = Parameter()
    numCPUs = luigi.IntParameter()
    indexfile = Parameter()
    spliceFile = Parameter()
    mappingLogFile = Parameter()
    outsam = Parameter()
    unalned = Parameter()
    bam_file = Parameter()
    sorted_bam_file = Parameter()
    map_ref = Parameter()
    workdir = Parameter()
    map_dir = Parameter()
    split_logfile = Parameter()
    kingdom = Parameter()

    def requires(self):
        """Require sort step to be finished."""
        return [SplitBAMfile(fastq1=self.fastq1,
                             fastq2=self.fastq2,
                             numCPUs=self.numCPUs,
                             indexfile=self.indexfile,
                             spliceFile=self.spliceFile,
                             mappingLogFile=self.mappingLogFile,
                             outsam=self.outsam,
                             bam_file=self.bam_file,
                             unalned=self.unalned,
                             sorted_bam_file=self.sorted_bam_file,
                             workdir=self.workdir,
                             map_ref=self.map_ref,
                             split_logfile=self.split_logfile),
                SplitRefProkEuk(fastq1=self.fastq1,
                                fastq2=self.fastq2,
                                numCPUs=self.numCPUs,
                                indexfile=self.indexfile,
                                spliceFile=self.spliceFile,
                                mappingLogFile=self.mappingLogFile,
                                outsam=self.outsam,
                                bam_file=self.bam_file,
                                unalned=self.unalned,
                                sorted_bam_file=self.sorted_bam_file,
                                workdir=self.workdir,
                                map_ref=self.map_ref,
                                map_dir=self.map_dir)]

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


class FeatureCount(ExternalProgramTask):
    """Count Features."""

    gff = luigi.Parameter()
    count_file = luigi.Parameter()
    bam_file = luigi.Parameter()
    bindir = luigi.Parameter()

    def requires(self):
        """Require reference fasta format file."""
        return [RefFile(self.gff), RefFile(self.bam_file)]

    def output(self):
        """Index output."""
        return LocalTarget(self.count_file)

    def program_args(self):
        """Run hisat2-build command."""
        return ["featureCounts",
                "-a", self.gff,
                "-s", 1,
                "-p", "-P",
                "-o", self.count_file,
                self.bam_file]

    def program_environment(self):
        """Environmental variables for this program."""
        return {'PATH': os.environ["PATH"] + ":" + self.bindir}


class AllProk(luigi.WrapperTask):
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
        yield GFF2GTF(gff_file=self.gff_file,
                      bindir=self.bindir, workdir=self.bindir)
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
            yield SAM2BAMfile(fastq1=trim_dir + "/" + samp + ".1.trimmed.fastq",
                              fastq2=trim_dir + "/" + samp + ".2.trimmed.fastq",
                              numCPUs=self.numCPUs,
                              indexfile=self.indexfile,
                              spliceFile=None,
                              mappingLogFile=map_dir + "/mapping.log",
                              unalned=map_dir + "/unligned.fastq",
                              outsam=map_dir + "/" + samp + ".sam",
                              bam_file=map_dir + "/" + samp + "_map.bam")
            yield SortBAMfile(fastq1=trim_dir + "/" + samp + ".1.trimmed.fastq",
                              fastq2=trim_dir + "/" + samp + ".2.trimmed.fastq",
                              numCPUs=self.numCPUs,
                              indexfile=self.indexfile,
                              spliceFile=None,
                              mappingLogFile=map_dir + "/mapping.log",
                              unalned=map_dir + "/unligned.fastq",
                              outsam=map_dir + "/" + samp + ".sam",
                              bam_file=map_dir + "/" + samp + "_map.bam",
                              sorted_bam_file=map_dir + "/" + samp + "_sort_map.bam")
            yield FeatureCount(gff=self.gtf_file,
                               count_file=map_dir + "/" + samp + "_prok.count",
                               bam_file=map_dir + "/" + samp + "_sort_map.bam",
                               bindir=self.bindir)


class AllEuk(luigi.WrapperTask):
    """From Mapping to Counting step for Eukaryotic reference."""

    gff_file = Parameter()
    gtf_file = Parameter()
    fastq_dic = luigi.DictParameter()
    numCPUs = luigi.IntParameter()
    indexfile = Parameter()
    bindir = Parameter()
    workdir = Parameter()

    def requires(self):
        """A wrapper for running the QC."""
        splice_list = [self.workdir + "/" +
                       f for f in os.listdir(self.workdir) if f.endswith('.splice')]
        if len(splice_list) > 1:
            splice_file = ','.join([splice_list])
        else:
            splice_file = splice_list[0]
        for samp, fastq in self.fastq_dic.iteritems():
            trim_dir = self.workdir + "/" + samp + "/trimming_results"
            map_dir = self.workdir + "/" + samp + "/mapping_results"
            if os.path.isdir(map_dir) is False:
                os.makedirs(map_dir)
            yield MapHisat(fastq1=trim_dir + "/" + samp + ".1.trimmed.fastq",
                           fastq2=trim_dir + "/" + samp + ".2.trimmed.fastq",
                           numCPUs=self.numCPUs,
                           indexfile=self.indexfile,
                           spliceFile=splice_file,
                           mappingLogFile=map_dir + "/mapping.log",
                           unalned=map_dir + "/unligned.fastq",
                           outsam=map_dir + "/" + samp + ".sam")
            yield SAM2BAMfile(fastq1=trim_dir + "/" + samp + ".1.trimmed.fastq",
                              fastq2=trim_dir + "/" + samp + ".2.trimmed.fastq",
                              numCPUs=self.numCPUs,
                              indexfile=self.indexfile,
                              spliceFile=splice_file,
                              mappingLogFile=map_dir + "/mapping.log",
                              unalned=map_dir + "/unligned.fastq",
                              outsam=map_dir + "/" + samp + ".sam",
                              bam_file=map_dir + "/" + samp + "_map.bam")
            yield SortBAMfile(fastq1=trim_dir + "/" + samp + ".1.trimmed.fastq",
                              fastq2=trim_dir + "/" + samp + ".2.trimmed.fastq",
                              numCPUs=self.numCPUs,
                              indexfile=self.indexfile,
                              spliceFile=splice_file,
                              mappingLogFile=map_dir + "/mapping.log",
                              unalned=map_dir + "/unligned.fastq",
                              outsam=map_dir + "/" + samp + ".sam",
                              bam_file=map_dir + "/" + samp + "_map.bam",
                              sorted_bam_file=map_dir + "/" + samp + "_sort_map.bam")
            yield FeatureCount(gff=self.gtf_file,
                               count_file=map_dir + "/" + samp + "_euk.count",
                               bam_file=map_dir + "/" + samp + "_sort_map.bam",
                               bindir=self.bindir)


class AllBoth(luigi.WrapperTask):
    """From Mapping to Counting step for Eukaryotic and Prokaryotic."""

    fastq_dic = luigi.DictParameter()
    numCPUs = luigi.IntParameter()
    indexfile = Parameter()
    prok_gff = Parameter()
    euk_gff = Parameter()
    bindir = Parameter()
    workdir = Parameter()

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
            yield MapHisat(fastq1=trim_dir + "/" + samp + ".1.trimmed.fastq",
                           fastq2=trim_dir + "/" + samp + ".2.trimmed.fastq",
                           numCPUs=self.numCPUs,
                           indexfile=self.indexfile,
                           spliceFile=splice_file,
                           mappingLogFile=map_dir + "/mapping.log",
                           unalned=map_dir + "/unligned.fastq",
                           outsam=map_dir + "/" + samp + ".sam")
            yield SAM2BAMfile(fastq1=trim_dir + "/" + samp + ".1.trimmed.fastq",
                              fastq2=trim_dir + "/" + samp + ".2.trimmed.fastq",
                              numCPUs=self.numCPUs,
                              indexfile=self.indexfile,
                              spliceFile=splice_file,
                              mappingLogFile=map_dir + "/mapping.log",
                              unalned=map_dir + "/unligned.fastq",
                              outsam=map_dir + "/" + samp + ".sam",
                              bam_file=map_dir + "/" + samp + "_map.bam")
            yield SortBAMfile(fastq1=trim_dir + "/" + samp + ".1.trimmed.fastq",
                              fastq2=trim_dir + "/" + samp + ".2.trimmed.fastq",
                              numCPUs=self.numCPUs,
                              indexfile=self.indexfile,
                              spliceFile=splice_file,
                              mappingLogFile=map_dir + "/mapping.log",
                              unalned=map_dir + "/unligned.fastq",
                              outsam=map_dir + "/" + samp + ".sam",
                              bam_file=map_dir + "/" + samp + "_map.bam",
                              sorted_bam_file=map_dir + "/" + samp + "_sort_map.bam")
            yield RefNames(fastq1=trim_dir + "/" + samp + ".1.trimmed.fastq",
                           fastq2=trim_dir + "/" + samp + ".2.trimmed.fastq",
                           numCPUs=self.numCPUs,
                           indexfile=self.indexfile,
                           spliceFile=splice_file,
                           mappingLogFile=map_dir + "/mapping.log",
                           unalned=map_dir + "/unligned.fastq",
                           outsam=map_dir + "/" + samp + ".sam",
                           bam_file=map_dir + "/" + samp + "_map.bam",
                           sorted_bam_file=map_dir + "/" + samp + "_sort_map.bam",
                           workdir=self.workdir,
                           map_ref=map_dir + "/" + samp + ".reflist")
            yield SplitBAMfile(fastq1=trim_dir + "/" + samp + ".1.trimmed.fastq",
                               fastq2=trim_dir + "/" + samp + ".2.trimmed.fastq",
                               numCPUs=self.numCPUs,
                               indexfile=self.indexfile,
                               spliceFile=splice_file,
                               mappingLogFile=map_dir + "/mapping.log",
                               unalned=map_dir + "/unligned.fastq",
                               outsam=map_dir + "/" + samp + ".sam",
                               bam_file=map_dir + "/" + samp + "_map.bam",
                               sorted_bam_file=map_dir + "/" + samp + "_sort_map.bam",
                               workdir=self.workdir,
                               map_ref=map_dir + "/" + samp + ".reflist",
                               split_logfile=map_dir + "/" + samp + ".splitlog")
            yield SplitRefProkEuk(fastq1=trim_dir + "/" + samp + ".1.trimmed.fastq",
                                  fastq2=trim_dir + "/" + samp + ".2.trimmed.fastq",
                                  numCPUs=self.numCPUs,
                                  indexfile=self.indexfile,
                                  spliceFile=splice_file,
                                  mappingLogFile=map_dir + "/mapping.log",
                                  unalned=map_dir + "/unligned.fastq",
                                  outsam=map_dir + "/" + samp + ".sam",
                                  bam_file=map_dir + "/" + samp + "_map.bam",
                                  sorted_bam_file=map_dir + "/" + samp + "_sort_map.bam",
                                  workdir=self.workdir,
                                  map_ref=map_dir + "/" + samp + ".reflist",
                                  map_dir=map_dir)
            yield MergeBAMfile(fastq1=trim_dir + "/" + samp + ".1.trimmed.fastq",
                               fastq2=trim_dir + "/" + samp + ".2.trimmed.fastq",
                               numCPUs=self.numCPUs,
                               indexfile=self.indexfile,
                               spliceFile=splice_file,
                               mappingLogFile=map_dir + "/mapping.log",
                               unalned=map_dir + "/unligned.fastq",
                               outsam=map_dir + "/" + samp + ".sam",
                               bam_file=map_dir + "/" + samp + "_map.bam",
                               sorted_bam_file=map_dir + "/" + samp + "_sort_map.bam",
                               workdir=self.workdir,
                               map_ref=map_dir + "/" + samp + ".reflist",
                               map_dir=map_dir,
                               split_logfile=map_dir + "/" + samp + ".splitlog",
                               kingdom='prokarya')
            yield MergeBAMfile(fastq1=trim_dir + "/" + samp + ".1.trimmed.fastq",
                               fastq2=trim_dir + "/" + samp + ".2.trimmed.fastq",
                               numCPUs=self.numCPUs,
                               indexfile=self.indexfile,
                               spliceFile=splice_file,
                               mappingLogFile=map_dir + "/mapping.log",
                               unalned=map_dir + "/unligned.fastq",
                               outsam=map_dir + "/" + samp + ".sam",
                               bam_file=map_dir + "/" + samp + "_map.bam",
                               sorted_bam_file=map_dir + "/" + samp + "_sort_map.bam",
                               workdir=self.workdir,
                               map_ref=map_dir + "/" + samp + ".reflist",
                               map_dir=map_dir,
                               split_logfile=map_dir + "/" + samp + ".splitlog",
                               kingdom='eukarya')
            yield FeatureCount(gff=self.prok_gff,
                               count_file=map_dir + "/" + samp + "_prok.count",
                               bam_file=map_dir + "/" + "prokarya.bam",
                               bindir=self.bindir)
            yield FeatureCount(gff=self.euk_gff,
                               count_file=map_dir + "/" + samp + "_euk.count",
                               bam_file=map_dir + "/" + "eukarya.bam",
                               bindir=self.bindir)
