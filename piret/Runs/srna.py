#! /usr/bin/env python

"""Check design."""
from __future__ import print_function
import os
from os.path import basename, isdir
from os import makedirs, listdir
import sys
from plumbum.cmd import samtools, grep, sed, bedtools, cat, awk
from piret.Runs import Map
from luigi.util import requires, inherits
import luigi
import logging
import pandas as pd


class ExtractPP(luigi.Task):
    """Extract properly paried reads."""
    kingdom = luigi.Parameter()
    workdir = luigi.Parameter()
    map_dir = luigi.Parameter()
    sample = luigi.Parameter()
    num_cpus = luigi.IntParameter()

    def output(self):
        """SAM file output of the mapping."""
        if self.kingdom in ['prokarya', 'eukarya']:
            bam_file = self.map_dir + "/" + self.sample + "_srt.bam"
            fbam = bam_file.split(".bam")[0] + ".fw_srt.bam"
            return luigi.LocalTarget(fbam)
        elif self.kingdom in ['both']:
            prok_bam_srt = self.map_dir + "/" + self.sample + "_srt_prok.fw_srt.bam"
            return luigi.LocalTarget(prok_bam_srt)

    def run(self):
        """Run split and merges."""
        if self.kingdom in ['prokarya', 'eukarya']:
            bam_file = self.map_dir + "/" + self.sample + "_srt.bam"
            fbam = bam_file.split(".bam")[0] + ".fw.bam"
            bbam = bam_file.split(".bam")[0] + ".bw.bam"
            self.prop_paired(bam_file)
            self.merge_prop_paired(bam_file, fbam, bbam)
            self.sort_bam(fbam.split(".bam")[0] + ".bam")
            self.sort_bam(bbam.split(".bam")[0] + ".bam")
        elif self.kingdom == "both":
            euk_bam_file = self.map_dir + "/" + self.sample + "_srt_euk.bam"
            prok_bam_file = self.map_dir + "/" + self.sample + "_srt_prok.bam"
            euk_fbam = euk_bam_file.split(".bam")[0] + ".fw.bam"
            prok_fbam = prok_bam_file.split(".bam")[0] + ".fw.bam"
            euk_bbam = euk_bam_file.split(".bam")[0] + ".bw.bam"
            prok_bbam = prok_bam_file.split(".bam")[0] + ".bw.bam"
            self.prop_paired(euk_bam_file)
            self.prop_paired(prok_bam_file)
            self.merge_prop_paired(euk_bam_file, euk_fbam, euk_bbam)
            self.merge_prop_paired(prok_bam_file, prok_fbam, prok_bbam)
            self.sort_bam(euk_fbam.split(".bam")[0] + ".bam")
            self.sort_bam(prok_fbam.split(".bam")[0] + ".bam")
            self.sort_bam(euk_bbam.split(".bam")[0] + ".bam")
            self.sort_bam(prok_bbam.split(".bam")[0] + ".bam")

    def prop_paired(self, bam_file):
        """Extract properly paired reads."""

        tmp_dir = os.path.join(self.workdir, "tmp")
        if os.path.isdir(tmp_dir) is False:
            os.makedirs(tmp_dir)
        for flag in ["-f99", "-f147", "-f163", "-f83"]:
            fn = basename(bam_file).split(".bam")[0] + "_" + flag.split("-")[1] + ".bam"  # name of the file
            options = ["view", "-h", "-b", flag, bam_file, "-o",
                       os.path.join(tmp_dir, fn)]
            samtools_cmd = samtools[options]
            logger = logging.getLogger('luigi-interface')
            logger.info(samtools_cmd)
            samtools_cmd()

    def merge_prop_paired(self, bam_file, fbam, bbam):
        """Merge properly paired files."""
        tmp_dir = os.path.join(self.workdir, "tmp")
        if os.path.isdir(tmp_dir) is False:
            os.makedirs(tmp_dir)
        f99 = basename(bam_file).split(".bam")[0] + "_" + "f99.bam"
        f147 = basename(bam_file).split(".bam")[0] + "_" + "f147.bam"
        options = ["merge", fbam, os.path.join(tmp_dir, f99),
                   os.path.join(tmp_dir, f147), "-f"]
        merge_cmd = samtools[options]
        merge_cmd()
        logger = logging.getLogger('luigi-interface')
        logger.info(merge_cmd)
        f163 = basename(bam_file).split(".bam")[0] + "_" + "f163.bam"
        f83 = basename(bam_file).split(".bam")[0] + "_" + "f83.bam"
        options = ["merge", bbam, os.path.join(tmp_dir, f163),
                   os.path.join(tmp_dir, f83), "-f"]
        merge_cmd = samtools[options]
        logger.info(merge_cmd)
        merge_cmd()

    def sort_bam(self, bam_file):
        """Sort BAM file."""
        sorted_bam_file = bam_file.split(".bam")[0] + "_srt.bam"
        options = ["sort", bam_file,
                   "-o", sorted_bam_file]
        samtools_cmd = samtools[options]
        logger = logging.getLogger('luigi-interface')
        logger.info(samtools_cmd)
        samtools_cmd()


class ExtractPPW(luigi.WrapperTask):
    """A wrapper task for mapping."""

    fastq_dic = luigi.DictParameter()
    indexfile = luigi.Parameter()
    workdir = luigi.Parameter()
    num_cpus = luigi.IntParameter()
    kingdom = luigi.Parameter()

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
            trim_dir = os.path.join(self.workdir, "processes", "qc", samp)
            map_dir = os.path.join(self.workdir, "processes", "mapping", samp)
            if os.path.isdir(map_dir) is False:
                os.makedirs(map_dir)
            if self.kingdom in ['prokarya', 'eukarya']:
                yield ExtractPP(
                                num_cpus=self.num_cpus,
                                map_dir=map_dir,
                                sample=samp,
                                kingdom=self.kingdom,
                                workdir=self.workdir)
            elif self.kingdom == 'both':
                # prok_gff = os.path.basename(self.gff_file.split(";")[0]).split(".gff")[0]
                # euk_gff = os.path.basename(self.gff_file.split(";")[1]).split(".gff")[0]
                yield ExtractPP(
                                num_cpus=self.num_cpus,
                                map_dir=map_dir,
                                sample=samp,
                                kingdom=self.kingdom,
                                workdir=self.workdir)


# @requires(ExtractPP)
class FindNovelRegions(luigi.Task):
    """create a new gff file based on coverage of intergenic regions"""

    kingdom = luigi.Parameter()
    gff_file = luigi.Parameter()
    map_dir = luigi.Parameter()
    sample = luigi.Parameter()
    workdir = luigi.Parameter()

    def output(self):
        """Check for presence of bedfile."""
        if self.kingdom in ['prokarya', 'eukarya']:
            novels = self.map_dir + "/" + self.sample + "_bw_novel.bedfile"
        elif self.kingdom == "both":
            novels = self.map_dir + "/" + self.sample + "_prok_bw_novel.bedfile"
        return luigi.LocalTarget(novels)


    def run(self):
        """Based on eukarya or prokarya or both, run these commands"""        
        if self.kingdom in ['prokarya', 'eukarya']:
            srt_bam_file = self.map_dir + "/" + self.sample + "_srt.bam"
            chrom_sizes = self.map_dir + "/" + self.sample + "_size.txt"
            fw_srt_bam_file = self.map_dir + "/" + self.sample + "_srt.fw_srt.bam"
            bw_srt_bam_file = self.map_dir + "/" + self.sample + "_srt.bw_srt.bam"
            fw_gcov = self.map_dir + "/" + self.sample + "_fw_gcov.bedfile"
            bw_gcov = self.map_dir + "/" + self.sample + "_bw_gcov.bedfile"
            fw_novels = self.map_dir + "/" + self.sample + "_fw_novel.bedfile"
            bw_novels = self.map_dir + "/" + self.sample + "_bw_novel.bedfile"
            out_gff = os.path.join(self.workdir, "processes", "no_region.gff")
            self.get_genome_ref(srt_bam_file, chrom_sizes)
            self.remove_region(self.gff_file, out_gff)
            self.genome_coverage(fw_srt_bam_file, fw_gcov)
            self.genome_coverage(bw_srt_bam_file, bw_gcov)
            self.novel_regions(out_gff, fw_gcov, fw_novels)
            self.novel_regions(out_gff, bw_gcov, bw_novels)
        elif self.kingdom == "both":
            prok_gff = self.gff_file.split(",")[0]
            euk_gff = self.gff_file.split(",")[1]
            euk_srt_bam_file = self.map_dir + "/" + self.sample + "_srt_euk.bam"
            euk_chrom_sizes = self.map_dir + "/" + self.sample + "_euk_size.txt"
            prok_srt_bam_file = self.map_dir + "/" + self.sample + "_srt_prok.bam"
            prok_chrom_sizes = self.map_dir + "/" + self.sample + "_size_prok.txt"
            euk_fw_srt_bam_file = self.map_dir + "/" + self.sample + "_srt_euk.fw_srt.bam"
            prok_fw_srt_bam_file = self.map_dir + "/" + self.sample + "_srt_prok.fw_srt.bam"
            euk_bw_srt_bam_file = self.map_dir + "/" + self.sample + "_srt_euk.bw_srt.bam"
            prok_bw_srt_bam_file = self.map_dir + "/" + self.sample + "_srt_prok.bw_srt.bam"
            euk_fw_gcov = self.map_dir + "/" + self.sample + "_euk_fw_gcov.bedfile"
            prok_fw_gcov = self.map_dir + "/" + self.sample + "_prok_fw_gcov.bedfile"
            euk_bw_gcov = self.map_dir + "/" + self.sample + "_euk_bw_gcov.bedfile"
            prok_bw_gcov = self.map_dir + "/" + self.sample + "_prok_bw_gcov.bedfile"
            euk_fw_novels = self.map_dir + "/" + self.sample + "_euk_fw_novel.bedfile"
            prok_fw_novels = self.map_dir + "/" + self.sample + "_prok_fw_novel.bedfile"
            euk_bw_novels = self.map_dir + "/" + self.sample + "_euk_bw_novel.bedfile"
            prok_bw_novels = self.map_dir + "/" + self.sample + "_prok_bw_novel.bedfile"
            euk_out_gff = os.path.join(self.workdir, "euk_no_region.gff")
            prok_out_gff = os.path.join(self.workdir, "prok_no_region.gff")
            self.get_genome_ref(euk_srt_bam_file, euk_chrom_sizes)
            self.get_genome_ref(prok_srt_bam_file, prok_chrom_sizes)
            self.remove_region(euk_gff, euk_out_gff)
            self.remove_region(prok_gff, prok_out_gff)
            self.genome_coverage(euk_fw_srt_bam_file, euk_fw_gcov)
            self.genome_coverage(prok_fw_srt_bam_file, prok_fw_gcov)
            self.genome_coverage(euk_bw_srt_bam_file, euk_bw_gcov)
            self.genome_coverage(prok_bw_srt_bam_file, prok_bw_gcov)
            self.novel_regions(gff=euk_out_gff, gcov=euk_fw_gcov, novels=euk_fw_novels)
            self.novel_regions(gff=prok_out_gff, gcov=prok_fw_gcov, novels=prok_fw_novels)
            self.novel_regions(gff=euk_out_gff, gcov=euk_bw_gcov, novels=euk_bw_novels)
            self.novel_regions(gff=prok_out_gff, gcov=prok_bw_gcov, novels=prok_bw_novels)

    def get_genome_ref(self, sorted_bam_file, chrom_sizes):
        """Calculate the size of each chromosome/contig."""
        chain = samtools["view", "-H", sorted_bam_file] | grep["@SQ"] | sed["s/@SQ.*SN://g"] | sed["s/LN://g"] > chrom_sizes
        chain()

    def genome_coverage(self, srt_bam, gcov):
        """Get the coverage info."""
        (bedtools["genomecov", "-split", "-bg", "-ibam", srt_bam] > gcov)()

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


class FindNovelRegionsW(luigi.WrapperTask):
    fastq_dic = luigi.DictParameter()
    workdir = luigi.Parameter()
    gff_file = luigi.Parameter()
    kingdom = luigi.Parameter()

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
            trim_dir = os.path.join(self.workdir, "processes", "qc", samp)
            map_dir = os.path.join(self.workdir, "processes", "mapping", samp)
            yield FindNovelRegions(map_dir=map_dir,
                                   sample=samp,
                                   workdir=self.workdir,
                                   gff_file=self.gff_file,
                                   kingdom=self.kingdom)


# @inherits(FindNovelRegionsW)
class CompileGFF(luigi.Task):
    """Compile novel regions and create a gff."""
    fastq_dic = luigi.DictParameter()
    kingdom = luigi.Parameter()
    workdir = luigi.Parameter()
    gff_file = luigi.Parameter()

    def output(self):
        if self.kingdom in ["prokarya", "eukarya"]:
            bw_novel = os.path.join(self.workdir, "processes", "novel", "updated.gff")
        elif self.kingdom == "both":
            bw_novel = self.workdir + "/" + "prok_updated.gff"
        return luigi.LocalTarget(bw_novel)

    def run(self):
        if self.kingdom in ["prokarya", "eukarya"]:
            map_dir = os.path.join(self.workdir, "processes", "mapping")
            if os.path.exists(os.path.join(self.workdir, "processes", "novel")) is False:
                os.makedirs(os.path.join(self.workdir, "processes", "novel"))
            fw_novel = os.path.join(self.workdir, "processes", "novel", "fw_all_novel.txt")
            bw_novel = os.path.join(self.workdir, "processes", "novel", "bw_all_novel.txt")
            fw_beds = [os.path.join(map_dir, samp, samp + "_fw_novel.bedfile") for samp in self.fastq_dic.keys()]
            bw_beds = [os.path.join(map_dir, samp, samp + "_bw_novel.bedfile") for samp in self.fastq_dic.keys()]
            fw_gff = os.path.join(self.workdir, "processes", "novel", "fw_all_novel.gff")
            bw_gff = os.path.join(self.workdir, "processes", "novel", "bw_all_novel.gff")
            final_gff = os.path.join(self.workdir, "novel", "updated.gff")

            self.compile_novel_regions(fw_beds, fw_novel)
            self.compile_novel_regions(bw_beds, bw_novel)
            self.make_gff(fw_novel, "+", fw_gff)
            self.make_gff(bw_novel, "-", bw_gff)
            self.concat_gff(self.gff_file, fw_gff, bw_gff, final_gff)
        elif self.kingdom == "both":
            euk_fw_novel = os.path.join(self.workdir, "euk_fw_all_novel.txt")
            prok_fw_novel = os.path.join(self.workdir, "prok_fw_all_novel.txt")
            euk_bw_novel = os.path.join(self.workdir, "euk_bw_all_novel.txt")
            prok_bw_novel = os.path.join(self.workdir, "prok_bw_all_novel.txt")
            euk_fw_beds = [os.path.join(self.workdir, samp, "mapping_results", samp +
                                    "_euk_fw_novel.bedfile") for samp in self.fastq_dic.keys()]
            prok_fw_beds = [os.path.join(self.workdir, samp, "mapping_results", samp +
                                    "_prok_fw_novel.bedfile") for samp in self.fastq_dic.keys()]
            euk_bw_beds = [os.path.join(self.workdir, samp, "mapping_results", samp +
                                    "_euk_bw_novel.bedfile") for samp in self.fastq_dic.keys()]
            prok_bw_beds = [os.path.join(self.workdir, samp, "mapping_results", samp +
                                    "_prok_bw_novel.bedfile") for samp in self.fastq_dic.keys()]
            euk_fw_gff = os.path.join(self.workdir, "euk_fw_all_novel.gff")
            prok_fw_gff = os.path.join(self.workdir, "prok_fw_all_novel.gff")
            euk_bw_gff = os.path.join(self.workdir, "euk_bw_all_novel.gff")
            prok_bw_gff = os.path.join(self.workdir, "prok_bw_all_novel.gff")
            euk_final_gff = os.path.join(self.workdir, "euk_updated.gff")
            prok_final_gff = os.path.join(self.workdir, "prok_updated.gff")

            self.compile_novel_regions(euk_fw_beds, euk_fw_novel)
            self.compile_novel_regions(prok_fw_beds, prok_fw_novel)
            self.compile_novel_regions(euk_bw_beds, euk_bw_novel)
            self.compile_novel_regions(prok_bw_beds, prok_bw_novel)
            self.make_gff(euk_fw_novel, "+", euk_fw_gff)
            self.make_gff(prok_fw_novel, "+", prok_fw_gff)
            self.make_gff(euk_bw_novel, "-", euk_bw_gff)
            self.make_gff(prok_bw_novel, "-", prok_bw_gff)
            self.concat_gff(self.gff_file.split(",")[1], euk_fw_gff, euk_bw_gff, euk_final_gff)
            self.concat_gff(self.gff_file.split(",")[0], prok_fw_gff, prok_bw_gff, prok_final_gff)




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

