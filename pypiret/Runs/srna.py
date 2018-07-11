#! /usr/bin/env python

"""Check design."""
from __future__ import print_function
import os
from os.path import basename, isdir
from os import makedirs, listdir
import sys
from plumbum.cmd import samtools, grep, sed, bedtools, cat, awk
from pypiret.Runs import Map
from luigi.util import requires, inherits
import luigi
import pandas as pd


@requires(Map.Hisat)
class ExtractPP(luigi.Task):
    """Extract properly paried reads."""
    

    def output(self):
        """SAM file output of the mapping."""
        bam_file = self.map_dir + "/" + self.sample + ".bam"
        fbam = bam_file.split(".bam")[0] + ".fw_srt.bam"
        return luigi.LocalTarget(fbam)

    def run(self):
        """Run split and merges."""
        bam_file = self.map_dir + "/" + self.sample + ".bam"
        fbam = bam_file.split(".bam")[0] + ".fw.bam"
        bbam = bam_file.split(".bam")[0] + ".bw.bam"
        self.prop_paired(self.sample, bam_file)
        self.merge_prop_paired(self.sample, fbam, bbam)
        self.sort_bam(fbam.split(".bam")[0] + ".bam")
        self.sort_bam(bbam.split(".bam")[0] + ".bam")

    def prop_paired(self, samp, bam_file):
        """Extract properly paired reads."""

        for flag in ["-f99", "-f147", "-f163", "-f83"]:
            options = ["view", "-h", "-b", flag, bam_file, "-o",
                       "/tmp/" + str(samp) + "_" + flag.split("-")[1] + ".bam"]
            samtools_cmd = samtools[options]
            samtools_cmd()

    def merge_prop_paired(self, samp, fbam, bbam):
        """Merge properly paired files."""
        options = ["merge", fbam,
                   "/tmp/" + samp + "_" + "f99.bam",
                   "/tmp/" + samp + "_" + "f147.bam", "-f"]
        merge_cmd = samtools[options]
        merge_cmd()        
        options = ["merge", bbam,
                   "/tmp/" + samp + "_" + "f163.bam",
                   "/tmp/" + samp + "_" + "f83.bam",
                   "-f"]
        merge_cmd = samtools[options]
        merge_cmd()

    def sort_bam(self, bam_file):
        """Sort BAM file."""
        sorted_bam_file = bam_file.split(".bam")[0] + "_srt.bam"
        options = ["sort", bam_file,
                   "-o", sorted_bam_file]
        samtools_cmd = samtools[options]
        samtools_cmd()


class ExtractPPW(luigi.WrapperTask):
    """A wrapper task for mapping."""

    fastq_dic = luigi.DictParameter()
    ref_file = luigi.Parameter()
    indexfile = luigi.Parameter()
    bindir = luigi.Parameter()
    workdir = luigi.Parameter()
    num_cpus = luigi.IntParameter()

    def requires(self):
        """A wrapper task for running mapping."""
        splist = [self.workdir + "/" +
                  f for f in listdir(self.workdir) if f.endswith('.splice')]
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
            yield ExtractPP(fastqs=[trim_dir + "/" + samp + ".1.trimmed.fastq",
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


@requires(ExtractPPW)
class SummarizeMap(luigi.Task):
    """Summarizes FaQC results of all samples into a table"""

    def output(self):
        """Mapping Summary Output."""
        out_file = self.workdir + "/" + "MapSummary.csv"
        return luigi.LocalTarget(out_file)

    def run(self):
        """Parse the FaQC stats."""
        summ_dic = {}
        for samp, fastq in self.fastq_dic.items():
            map_dir = self.workdir + "/" + samp + "/mapping_results"
            filename = map_dir + "/" + "mapping.log"
            with open(filename, 'r') as file:
                lines = file.readlines()
                total_reads = lines[0].split("reads")[0].strip()
                con_unaligned = lines[2].split("(")[0].strip()
                con_aligned = lines[3].split("(")[0].strip()
                multi_aligned = lines[4].split("(")[0].strip()
                summ_dic[samp] = [total_reads,
                                  con_unaligned,
                                  con_aligned,
                                  multi_aligned]
        summ_table = pd.DataFrame.from_dict(summ_dic, orient='index')
        summ_table.columns = ["Paired reads", "Concordantly unaligned",
                              "Concordantly aligned", "Multi aligned"]
        out_file = self.workdir + "/" + "MapSummary.csv"
        summ_table.to_csv(out_file)


@requires(ExtractPP)
class FindNovelRegions(luigi.Task):
    """create a new gff file based on coverage of intergenic regions"""

    gff = luigi.Parameter()
    workdir = luigi.Parameter()

    def output(self):
        novels = self.map_dir + "/" + self.sample + "_bw_novel.bedfile"
        return luigi.LocalTarget(novels)

    def run(self):
        srt_bam_file = self.map_dir + "/" + self.sample + "_srt.bam"
        chrom_sizes = self.map_dir + "/" + self.sample + "_size.txt"
        fw_srt_bam_file = self.map_dir + "/" + self.sample + ".fw_srt.bam"
        bw_srt_bam_file = self.map_dir + "/" + self.sample + ".bw_srt.bam"
        fw_gcov = self.map_dir + "/" + self.sample + "_fw_gcov.bedfile"
        bw_gcov = self.map_dir + "/" + self.sample + "_bw_gcov.bedfile"
        fw_novels = self.map_dir + "/" + self.sample + "_fw_novel.bedfile"
        bw_novels = self.map_dir + "/" + self.sample + "_bw_novel.bedfile"
        out_gff = os.path.join(self.workdir, "no_region.gff")

        self.get_genome_ref(srt_bam_file, chrom_sizes)
        self.remove_region(self.gff, out_gff)
        self.genome_coverage(fw_srt_bam_file, chrom_sizes, fw_gcov)
        self.genome_coverage(bw_srt_bam_file, chrom_sizes, bw_gcov)
        self.novel_regions(out_gff, fw_gcov, fw_novels)
        self.novel_regions(out_gff, bw_gcov, bw_novels)

    def get_genome_ref(self, sorted_bam_file, chrom_sizes):
        """Calculate the size of each chromosome/contig."""
        (samtools["view", "-H", sorted_bam_file] | grep["@SQ"] |
         sed["s/@SQ.*SN://g"] | sed["s/LN://g"] > chrom_sizes)()

    def genome_coverage(self, srt_bam, chrom_sizes, gcov):
        """Get the coverage info."""
        (bedtools["genomecov", "-split", "-bg", "-ibam", srt_bam, "-g",
                  chrom_sizes] > gcov)()

    def remove_region(self, in_gff, out_gff):
        """Remove region feature from gff table."""
        
        with open(out_gff, 'w') as o:
            with open(in_gff, 'r') as gff:
                for line in gff:
                    if not line.startswith("#"):
                        if line.split("\t")[2] != "region":
                            o.write(line)

    def novel_regions(self, gff, gcov, novels):
        """Get regions that are novel, not found in gff"""
        (bedtools["subtract", "-A", "-a", gcov, "-b", gff] | bedtools["merge"] >
         novels)()


class FindNovelRegionsW(luigi.Task):
    fastq_dic = luigi.DictParameter()
    ref_file = luigi.Parameter()
    indexfile = luigi.Parameter()
    bindir = luigi.Parameter()
    workdir = luigi.Parameter()
    num_cpus = luigi.IntParameter()
    gff_file = luigi.Parameter()
    workdir = luigi.Parameter()

    def requires(self):
        """A wrapper task for running mapping."""
        splist = [self.workdir + "/" +
                  f for f in listdir(self.workdir) if f.endswith('.splice')]
        if len(splist) > 1:
            splice_file = ','.join(splist)
        elif len(splist) == 1:
            splice_file = splist[0]
        else:
            splice_file = ''
        for samp, fastq in self.fastq_dic.items():
            trim_dir = self.workdir + "/" + samp + "/trimming_results"
            map_dir = self.workdir + "/" + samp + "/mapping_results"
            yield FindNovelRegions(fastqs=[trim_dir + "/" + samp + ".1.trimmed.fastq",
                                           trim_dir + "/" + samp + ".2.trimmed.fastq"],
                                   qc_outdir=trim_dir,
                                   num_cpus=self.num_cpus,
                                   indexfile=self.indexfile,
                                   spliceFile=splice_file,
                                   outsam=map_dir + "/" + samp + ".sam",
                                   ref_file=self.ref_file,
                                   bindir=self.bindir,
                                   map_dir=map_dir,
                                   sample=samp,
                                   workdir=self.workdir,
                                   gff=self.gff_file)


@inherits(FindNovelRegionsW)
class CompileGFF(luigi.Task):
    """Compile novel regions and create a gff."""

    def output(self):
        bw_novel = self.workdir + "/" + "updated.gff"
        return luigi.LocalTarget(bw_novel)


    def run(self):

        fw_novel = os.path.join(self.workdir, "fw_all_novel.txt")
        bw_novel = os.path.join(self.workdir, "bw_all_novel.txt")
        fw_beds = [os.path.join(self.workdir, samp, "mapping_results", samp +
                                "_fw_novel.bedfile") for samp in self.fastq_dic.keys()]
        bw_beds = [os.path.join(self.workdir, samp, "mapping_results", samp +
                                "_bw_novel.bedfile") for samp in self.fastq_dic.keys()]
        fw_gff = os.path.join(self.workdir, "fw_all_novel.gff")
        bw_gff = os.path.join(self.workdir, "bw_all_novel.gff")
        final_gff = os.path.join(self.workdir, "updated.gff")

        self.compile_novel_regions(fw_beds, fw_novel)
        self.compile_novel_regions(bw_beds, bw_novel)
        self.make_gff(fw_novel, "+", fw_gff)
        self.make_gff(bw_novel, "-", bw_gff)
        self.concat_gff(self.gff_file, fw_gff, bw_gff, final_gff)


    def compile_novel_regions(self, bedfiles, all_novel):
        """get the novel regions that has at least 1 coverage and compile as 1"""
        (cat[bedfiles] | bedtools["sort"] |bedtools["merge"] > all_novel)()


    def make_gff(self, novel, strand, out_gff):
        """convert the novel regions to gff"""
        out = open(out_gff, 'w')
        out.write("##gff-version 3\n")
        out.write("##bed_to_gff_converter.py\n")
        i = 0
        for i, line in enumerate(open(novel)):
            line = line.rstrip('\r\n')
            if line and not line.startswith('#'):
                elems = line.split('\t')
                start = int(elems[1]) + 1
                end = int(elems[2])
                feature = "novel_region"
                score = "."
                if strand == "+":
                    group = 'ID=%s;Name=%s' % ("novel_regions_plus_" + str(i), "novel_regions_" + str(i))
                elif strand == "-":
                    group = 'ID=%s;Name=%s' % ("novel_regions_nega_" + str(i), "novel_regions_" + str(i))
                out.write('%s\tbed2gff\t%s\t%d\t%d\t%s\t%s\t.\t%s;\n'
                    % (elems[0], feature, start, end, score, strand, group))
        out.close()

    def concat_gff(self, orig_gff, novelplus_gff, novelneg_gff, new_gff):
        """concatenate gff."""
        (cat[orig_gff, novelplus_gff, novelneg_gff] > new_gff)()

