#! /usr/bin/env python

"""Luigi Tasks to perform various RNA seq functions, from Mapping to Counting.

Mapping is done using hisat2 and counting is done using featurecounts
and stringtie
"""
import os
import sys
import logging
import luigi
import pandas as pd
from piret.miscs import RefFile
from collections import defaultdict as dd, Counter
from sys import stderr, exit
from plumbum.cmd import STAR
from plumbum.cmd import samtools, stringtie, mv, awk
from plumbum.cmd import hisat2
from luigi.util import inherits, requires
from luigi import Parameter, IntParameter, DictParameter, ListParameter
from luigi import LocalTarget
from luigi import ExternalTask
from luigi.contrib.external_program import ExternalProgramTask
dir_path = os.path.dirname(os.path.realpath(__file__))
lib_path = os.path.abspath(os.path.join(dir_path, '..'))
sys.path.append(lib_path)


class STARindex(luigi.Task):
    """ Creates STAR index."""
    fasta = Parameter()
    num_cpus = IntParameter()
    gff_file = Parameter()
    stardb_dir = Parameter()
    kingdom = Parameter()

    def requires(self):
        if ',' in self.fasta:
            return[RefFile(self.fasta.split(",")[1])]
        else:
            return[RefFile(self.fasta)]

    def output(self):
        """Expected index output"""
        return LocalTarget(os.path.join(self.stardb_dir, "chrName.txt"))

    def run(self):
        if os.path.exists(self.stardb_dir) is False:
            os.makedirs(self.stardb_dir)
        if ',' in self.fasta:
            fnas = self.fasta.split(",")
            gffs = self.gff_file.split(",")
        else:
            fnas = [self.fasta]
            gffs = [self.gff_file]
        fnas.insert(0, "--genomeFastaFiles")
        gffs.insert(0, "--sjdbGTFfile")
        if self.kingdom == "eukarya":
            ind_opt = ["--runMode", "genomeGenerate",
                       "--runThreadN", self.num_cpus,
                       "--genomeDir", self.stardb_dir]
        elif self.kingdom == "prokarya":
            ind_opt = ["--runMode", "genomeGenerate",
                       "--runThreadN", self.num_cpus,
                       "--genomeDir", self.stardb_dir]
        elif self.kingdom == "both":
            ind_opt = ["--runMode", "genomeGenerate",
                       "--runThreadN", self.num_cpus,
                       "--genomeDir", self.stardb_dir]
        for f in fnas:
            ind_opt.append(f)
        star_cmd = STAR[ind_opt]
        logger = logging.getLogger('luigi-interface')
        logger.info(star_cmd)
        star_cmd()


class map_star(luigi.Task):
    """Mapping the QCed sequences to reference using star aligner."""
    fastqs = ListParameter()
    stardb_dir = Parameter()
    map_dir = Parameter()
    sample = Parameter()
    num_cpus = IntParameter()
    align_intron_min = luigi.IntParameter()
    align_intron_max = luigi.IntParameter()

    # def requires(self):
    #     """See if input file exist."""
    #     return [luigi.LocalTarget(self.fastqs[0]),
    #             luigi.LocalTarget(self.fastqs[1])]

    def output(self):
        """SAM file output of the mapping."""
        bam_file = os.path.join(self.map_dir, self.sample) + "_srt.bam"
        return luigi.LocalTarget(bam_file)

    def run(self):
        """Run star"""
        if os.path.exists(self.map_dir) is False:
            os.makedirs(self.map_dir)
        star_option = ["--genomeDir", self.stardb_dir,
                       "--runThreadN", self.num_cpus,
                       "--outFileNamePrefix", os.path.join(
                           self.map_dir, self.sample + "_"),
                       "--outSAMtype", "BAM", "SortedByCoordinate",
                       "--alignIntronMax", self.align_intron_max,
                       "alignIntronMin", self.align_intron_min,
                       "--readFilesIn", self.fastqs[0], self.fastqs[1]]
        star_cmd = STAR[star_option]
        star_cmd()
        bam_file = os.path.join(self.map_dir, self.sample) + \
            "_Aligned.sortedByCoord.out.bam"
        mv_list = [bam_file, os.path.join(
            self.map_dir, self.sample) + "_srt.bam"]
        mv_cmd = mv[mv_list]
        mv_cmd()


class map_starW(luigi.WrapperTask):
    """A wrapper task for mapping using star."""

    fastq_dic = luigi.DictParameter()
    stardb_dir = luigi.Parameter()
    workdir = luigi.Parameter()
    num_cpus = luigi.IntParameter()

    def requires(self):
        """A wrapper task for mapping using STAR."""
        for samp, fastq in self.fastq_dic.items():
            trim_dir = os.path.join(self.workdir, "processes", "qc", samp)
            map_dir = os.path.join(self.workdir, "processes", "mapping", samp)
            if os.path.isdir(map_dir) is False:
                os.makedirs(map_dir)
            yield map_star(fastqs=[trim_dir + "/" + samp + ".1.trimmed.fastq",
                                   trim_dir + "/" + samp + ".2.trimmed.fastq"],
                           stardb_dir=self.stardb_dir, map_dir=map_dir,
                           sample=samp, num_cpus=self.num_cpus)


# @ requires(map_starW)
class SummarizeStarMap(luigi.Task):
    """Summarizes mapping results of all samples into a table"""

    fastq_dic = luigi.DictParameter()
    workdir = luigi.Parameter()

    def requires(self):
        for samp, fastq in self.fastq_dic.items():
            out_file = os.path.join(
                self.workdir, "processes", "mapping", samp, samp+"_Log.final.out")
            return RefFile(out_file)

    def output(self):
        """Mapping Summary Output."""
        out_file = os.path.join(self.workdir, "processes",
                                "mapping", "MapSummary.csv")
        return luigi.LocalTarget(out_file)

    def run(self):
        """Parse the mapping stats."""
        summ_dic = {}
        for samp, fastq in self.fastq_dic.items():
            map_dir = os.path.join(self.workdir, "processes", "mapping", samp)
            filename = map_dir + "/" + samp + "_Log.final.out"
            with open(filename, 'r') as file:
                lines = file.readlines()
                total_reads = lines[5].split("|")[1].strip()
                unq_map_reads = lines[8].split("|")[1].strip()
                multi_map_reads = int(lines[23].split("|")[1].strip())
                multix_map_reads = int(lines[25].split("|")[1].strip())
                mm_reads = multi_map_reads + multix_map_reads
                unmap_rds1 = int(lines[28].split("|")[1].strip())
                unmap_rds2 = int(lines[30].split("|")[1].strip())
                unmap_rds3 = int(lines[32].split("|")[1].strip())
                unmap_rds = unmap_rds1 + unmap_rds2 + unmap_rds3
                summ_dic[samp] = [total_reads,
                                  unq_map_reads,
                                  mm_reads, unmap_rds]
        sm_tbl = pd.DataFrame.from_dict(summ_dic, orient='index')
        sm_tbl.columns = ["Paired reads", "Uniquely mapped",
                          "multi mapped",
                          "Unmapped"]
        sm_tbl = sm_tbl.astype('int32')
        sm_tbl["perc_aln"] = (sm_tbl["Uniquely mapped"] /
                              sm_tbl["Paired reads"])*100
        sm_tbl["perc_unaln"] = 100-sm_tbl["perc_aln"]
        out_file = os.path.join(self.workdir, "processes",
                                "mapping", "MapSummary.csv")
        sm_tbl.to_csv(out_file)
