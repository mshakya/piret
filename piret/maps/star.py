#! /usr/bin/env python

"""Luigi Tasks to perform various RNA seq functions, from Mapping to Counting.

Mapping is done using hisat2 and counting is done using featurecounts
and stringtie
"""

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
from plumbum.cmd import hisat2
from plumbum.cmd import samtools, stringtie, mv, awk
from plumbum.cmd import STAR
import pandas as pd
from sys import stderr, exit
import logging
from collections import defaultdict as dd, Counter
from piret.miscs  import RefFile




class STARindex(luigi.Task):
    """ Creates STAR index."""
    fasta = Parameter()
    num_cpus = IntParameter()
    gff_file = Parameter()
    stardb_dir= Parameter()
    kingdom = Parameter()

    def requires(self):
        if self.kingdom == "both":
            return[RefFile(self.fasta.split(",")[1])]
        else:
            return[RefFile(self.fasta)]

    def output(self):
        """Expected index output"""
        return LocalTarget(os.path.join(self.stardb_dir, "chrName.txt"))
    
    def run(self):
        if os.path.exists(self.stardb_dir) is False:
            os.makedirs(self.stardb_dir)
        if self.kingdom == "eukarya":
            ind_opt = ["--runMode", "genomeGenerate",
                       "--runThreadN", self.num_cpus,
                       "--genomeDir", self.stardb_dir,
                       "--genomeFastaFiles", self.fasta,
                       "--sjdbGTFfile", self.gff_file]
        elif self.kingdom == "prokarya":
            ind_opt = ["--runMode", "genomeGenerate",
                       "--runThreadN", self.num_cpus,
                       "--genomeDir", self.stardb_dir,
                       "--genomeFastaFiles", self.fasta]
        elif self.kingdom == "both":
            ind_opt = ["--runMode", "genomeGenerate",
                       "--runThreadN", self.num_cpus,
                       "--genomeDir", self.stardb_dir,
                       "--genomeFastaFiles", self.fasta.split(",")[0],
                       self.fasta.split(",")[1]]

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
                       "--outFileNamePrefix", os.path.join(self.map_dir, self.sample + "_"),
                       "--outSAMtype", "BAM", "SortedByCoordinate",
                       "--alignIntronMax", self.align_intron_max,
                       "alignIntronMin", self.align_intron_min,
                       "--readFilesIn", self.fastqs[0], self.fastqs[1]]
        star_cmd = STAR[star_option]
        star_cmd()
        bam_file = os.path.join(self.map_dir, self.sample) + "_Aligned.sortedByCoord.out.bam"
        mv_list = [bam_file, os.path.join(self.map_dir, self.sample) + "_srt.bam"]
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
            map_dir = os.path.join(self.workdir, "processes", "mapping", samp)
            if os.path.isdir(map_dir) is False:
                os.makedirs(map_dir)
            yield map_star(fastqs=[fastq.split(":")[0],
                                   fastq.split(":")[1]],
                                   stardb_dir=self.stardb_dir, map_dir=map_dir,
                                   sample=samp, num_cpus=self.num_cpus)

@requires(map_starW)
class SummarizeStarMap(luigi.Task):
    """Summarizes mapping results of all samples into a table"""

    def output(self):
        """Mapping Summary Output."""
        out_file = self.workdir + "/" + "MapSummary.csv"
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
                multi_map_reads = lines[23].split("|")[1].strip()
                multiX_map_reads = lines[25].split("|")[1].strip()
                summ_dic[samp] = [total_reads,
                                  unq_map_reads,
                                  multi_map_reads,
                                  multiX_map_reads]
        summ_table = pd.DataFrame.from_dict(summ_dic, orient='index')
        summ_table.columns = ["Paired reads", "Uniquely mapped reads",
                              "mapped to multiple loci", "mapped to too many loci"]
        out_file = self.workdir + "/" + "MapSummary.csv"
        summ_table.to_csv(out_file)
