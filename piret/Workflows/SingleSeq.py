from __future__ import print_function
import os
import sys
import argparse
import luigi
dir_path = os.path.dirname(os.path.realpath(__file__))
lib_path = os.path.abspath(os.path.join(dir_path, '..'))
bin_path = os.path.join(lib_path, 'bin')
sys.path.append(lib_path)
os.environ["PATH"] += os.pathsep + bin_path
from piret.checks.Design import CheckDesign
from piret.summary import summarize
from piret.qc import FaQC
from piret.maps import hisat2
from piret.maps import star
from piret.dge import edgeR
from piret.dge import DESeq2
from piret.dge import ballgown
from piret.counts import stringtie
from piret.counts import featurecounts
from piret.novel import srna
from piret.functions.function import RunEmapper, GetAAs
from piret.pathways.opaver import RunOpaver
from luigi.interface import build


class SingleSeq:
    """Class that pieces luigi task for prokarya."""

    def __init__(self, qc, fastq_dic, ref_fasta, num_cpus,
                 local_scheduler, hisat_index, stardb_dir, workdir, kingdom,
                 no_of_jobs, exp_desn_file,
                 p_value, aligner, gff_file, pathway):
        self.qc = qc
        self.ref_fasta = ref_fasta
        self.fastq_dic = fastq_dic
        self.num_cpus = num_cpus
        self.hisat_index = hisat_index
        self.stardb_dir = stardb_dir
        self.workdir = workdir
        self.kingdom = kingdom
        self.local_scheduler = local_scheduler
        self.no_of_jobs = no_of_jobs,
        self.exp_desn_file = exp_desn_file
        self.p_value = p_value
        self.aligner = aligner
        self.gff_file = gff_file
        self.pathway = pathway

    def run_faqc(self):
        """A function that calls the FaQC function.

        it returns QCed files in respective directory
        """
        if self.qc is True:
            build([FaQC.SummarizeQC(fastq_dic=self.fastq_dic,
                                num_cpus=self.num_cpus,
                                workdir=self.workdir)],
                local_scheduler=self.local_scheduler,
                workers=1)
            qc_dic = {}
            for samp, path in self.fastq_dic.items():
                trim_dir = os.path.join(self.workdir, "processes", "qc", samp)
                qc_dic[samp] = trim_dir + "/" + samp + ".1.trimmed.fastq" + ":" + \
                               trim_dir + "/" + samp + ".2.trimmed.fastq"          
            return qc_dic

        else:
            return self.fastq_dic

    def create_db(self):
        """Function to create hisat or STAR index."""
        if self.aligner in ["hisat2", "hisat"]:
            build([hisat2.HisatIndex(fasta=self.ref_fasta,
                                  hi_index=self.hisat_index,
                                  num_cpus=self.num_cpus)],
                  local_scheduler=self.local_scheduler)
        elif self.aligner in ["star", "STAR"]:
            build([star.STARindex(fasta=self.ref_fasta,
                                 num_cpus=self.num_cpus,
                                 gff_file=self.gff_file,
                                 stardb_dir=self.stardb_dir,
                                 kingdom=self.kingdom)],
                  local_scheduler=self.local_scheduler)

    def map_reads(self, qc_dic):
        """Function to map reads."""
        if self.aligner == "hisat2":
            build([hisat2.HisatMapW(fastq_dic=qc_dic, num_cpus=self.num_cpus,
                                 indexfile=self.hisat_index, workdir=self.workdir)],
              local_scheduler=self.local_scheduler)
        elif self.aligner in ["STAR", "star"]:
            build([star.map_starW(fastq_dic=qc_dic, num_cpus=self.num_cpus,
                                 stardb_dir=self.stardb_dir, workdir=self.workdir)],
              local_scheduler=self.local_scheduler)

    def map_summarize(self):
        """Summarize mapped reads into a table."""
        if self.aligner == "hisat2":
            build([hisat2.SummarizeHisatMap(fastq_dic=self.fastq_dic,
                                    workdir=self.workdir,
                                    indexfile=self.hisat_index,
                                    num_cpus=self.num_cpus)],
            local_scheduler=self.local_scheduler, workers=1)
        elif self.aligner in ["STAR", "star"]:
            build([star.SummarizeStarMap(fastq_dic=self.fastq_dic,
                                    workdir=self.workdir,
                                    stardb_dir=self.stardb_dir,
                                    num_cpus=self.num_cpus)],
            local_scheduler=self.local_scheduler, workers=1)

    def extract_pp(self):
        """Extract properly paired reads."""
        build([srna.ExtractPPW(fastq_dic=self.fastq_dic, num_cpus=self.num_cpus,
                                 indexfile=self.hisat_index, workdir=self.workdir,
                                 kingdom=self.kingdom)],
              local_scheduler=self.local_scheduler)

    def NovelRegions(self):
        """Find novel regions."""
        build([srna.FindNovelRegionsW(fastq_dic=self.fastq_dic,
                                  workdir=self.workdir,
                                 kingdom=self.kingdom,
                                 gff_file=self.gff_file)],
              local_scheduler=self.local_scheduler)

    def create_new_gff(self):
        build([srna.CompileGFF(fastq_dic=self.fastq_dic,
                               kingdom=self.kingdom,
                               workdir=self.workdir,
                               gff_file=self.gff_file)],
          local_scheduler=self.local_scheduler, workers=1)


    def map_hisat_summarize(self):
        no_of_jobs = 1
        build([srna.SummarizeHisatMap(fastq_dic=self.fastq_dic,
                                  num_cpus=self.num_cpus,
                                  indexfile=self.hisat_index,
                                  workdir=self.workdir,
                                  ref_file=self.ref_fasta,
                                  kingdom=self.kingdom)],
         local_scheduler=self.local_scheduler, workers=no_of_jobs)

    def find_NovelRegions(self):
        build([srna.FindNovelRegionsW(fastq_dic=self.fastq_dic,
                                      num_cpus=self.num_cpus,
                                      indexfile=self.hisat_index,
                                      kingdom=self.kingdom,
                                      workdir=self.workdir,
                                      ref_file=self.ref_fasta,
                                      gff_file=self.gff_file)],
        local_scheduler = self.local_scheduler,
        workers=1)

    def feature_count(self, new_gff):
        build([featurecounts.FeatureCounts(fastq_dic=self.fastq_dic,
                                  num_cpus=self.num_cpus,
                                  gff_file=new_gff,
                                  indexfile=self.hisat_index,
                                  kingdom=self.kingdom,
                                  workdir=self.workdir,
                                  ref_file=self.ref_fasta)],
        local_scheduler=self.local_scheduler, workers=1)

    def feature_count_updated(self, new_gff):
        build([featurecounts.FeatureCounts(fastq_dic=self.fastq_dic,
                                  num_cpus=self.num_cpus,
                                  gff_file=new_gff,
                                  indexfile=self.hisat_index,
                                  kingdom=self.kingdom,
                                  workdir=self.workdir,
                                  ref_file=self.ref_fasta)],
        local_scheduler=self.local_scheduler, workers=1)

    def run_stringtie(self):
        build([Map.StringTieScoresW(fastq_dic=self.fastq_dic,
                                    num_cpus=self.num_cpus,
                                    indexfile=self.hisat_index,
                                    workdir=self.workdir,
                                    ref_file=self.ref_fasta,
                                    gff_file=os.path.join(self.workdir, "updated.gff"),
                                    kingdom=self.kingdom)],
                          local_scheduler=self.local_scheduler,
                          workers=1)

    def run_edger(self):
        build([edgeR.edgeR(kingdom=self.kingdom,
                         workdir=self.workdir,
                        #  gff_file=self.gff_file,
                        #  pathway=pathway,
                        #  GAGE=GAGE,
                         exp_design=self.exp_desn_file,
                         p_value=self.p_value)],
        local_scheduler=self.local_scheduler, workers=1)

    def run_deseq2(self):
        build([DESeq2.DESeq2(workdir=self.workdir,
                          kingdom=self.kingdom,
                        #   gff_file=self.gff_file,
                        #   pathway=pathway,
                        #   GAGE=GAGE,
                          exp_design=self.exp_desn_file,
                          p_value=self.p_value)],
        local_scheduler=self.local_scheduler, workers=1)

    def merge_stringtie(self, new_gff):
        build([stringtie.MergeStringTies(fastq_dic=self.fastq_dic,
                                    num_cpus=self.num_cpus,
                                    workdir=self.workdir,
                                    gff_file=new_gff,
                                    kingdom=self.kingdom)],
                local_scheduler=self.local_scheduler, workers=1)

    def restringtie(self):
        build([stringtie.ReStringTieScoresW(fastq_dic=self.fastq_dic,
                                       num_cpus=self.num_cpus,
                                       workdir=self.workdir,
                                       kingdom=self.kingdom)],
              local_scheduler=self.local_scheduler,
              workers=1)

    def run_ballgown(self):
        build([ballgown.ballgown(
                            workdir=self.workdir,
                            kingdom=self.kingdom,
                            exp_design=self.exp_desn_file,
                            p_value=self.p_value)],
              local_scheduler=self.local_scheduler, workers=1)

    def run_emapper(self, new_gff):
        build([RunEmapper(workdir=self.workdir,
                          gff_file=new_gff,
                          fasta_file=self.ref_fasta,
                          kingdom=self.kingdom)],
              local_scheduler=self.local_scheduler, workers=1)

    def run_opaver(self, method):
        build([RunOpaver(workdir=self.workdir,
                         kingdom=self.kingdom,
                         method=method)],
              local_scheduler=self.local_scheduler, workers=1)

    def summ_json(self, new_gff, method, NovelRegions):
        build([summarize.conver2json(gff_file=new_gff,
                                     fasta_file=self.ref_fasta,
                                     pathway=self.pathway,
                                     workdir=self.workdir,
                                     kingdom=self.kingdom,
                                     method=method,
                                     NovelRegions=NovelRegions)],
              local_scheduler=self.local_scheduler, workers=1)
