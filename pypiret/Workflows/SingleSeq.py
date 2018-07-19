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
from pypiret import CheckDesign, FaQC, Map, Summ, DGE, srna
from luigi.interface import build



class SingleSeq:
    """Class that pieces luigi task for prokarya."""

    def __init__(self, fastq_dic, ref_fasta, ref_gff, num_cpus,
                 local_scheduler, hisat_index, workdir, kingdom,
                 no_of_jobs, bindir, exp_desn_file,
                 p_value ):
        self.ref_fasta = ref_fasta
        self.ref_gff = ref_gff
        self.fastq_dic = fastq_dic
        self.num_cpus = num_cpus
        self.hisat_index = hisat_index
        self.workdir = workdir
        self.kingdom = kingdom
        self.local_scheduler = local_scheduler
        self.no_of_jobs = no_of_jobs
        self.bindir = bindir
        self.exp_desn_file = exp_desn_file
        self.p_value = p_value


    def create_hisat_index(self):
        """Function to create hisat index."""
        build([Map.HisatIndex(fasta=self.ref_fasta, hi_index=self.hisat_index,
               num_cpus=self.num_cpus)],
            local_scheduler=self.local_scheduler)

    def map_summarize(self):
        no_of_jobs = 1
        build([srna.SummarizeMap(fastq_dic=self.fastq_dic,
                                  num_cpus=self.num_cpus,
                                  indexfile=self.hisat_index,
                                  workdir=self.workdir,
                                  ref_file=self.ref_fasta,
                                  kingdom=self.kingdom)],
         local_scheduler=self.local_scheduler, workers=no_of_jobs)

    def find_novel_regions(self):
        build([srna.FindNovelRegionsW(fastq_dic=self.fastq_dic,
                                      num_cpus=self.num_cpus,
                                      indexfile=self.hisat_index,
                                      kingdom=self.kingdom,
                                      workdir=self.workdir,
                                      ref_file=self.ref_fasta,
                                      gff_file=self.ref_gff)],
        local_scheduler = self.local_scheduler,
        workers=self.no_of_jobs)

    def create_new_gff(self):
        build([srna.CompileGFF(fastq_dic=self.fastq_dic,
                               num_cpus=self.num_cpus,
                               indexfile=self.hisat_index,
                               kingdom=self.kingdom,
                               workdir=self.workdir,
                               ref_file=self.ref_fasta,
                               gff_file=self.ref_gff)],
          local_scheduler=self.local_scheduler, workers=self.no_of_jobs)

    def feature_count(self):
        build([Summ.FeatureCounts(fastq_dic=self.fastq_dic,
                                  num_cpus=self.num_cpus,
                                  gff=os.path.join(self.workdir, "updated.gff"),
                                  indexfile=self.hisat_index,
                                  kingdom=self.kingdom,
                                  workdir=self.workdir,
                                  ref_file=self.ref_fasta)],
        local_scheduler=self.local_scheduler, workers=self.no_of_jobs)

    def run_stringtie(self):
        build([Map.StringTieScoresW(fastq_dic=self.fastq_dic,
                                    num_cpus=self.num_cpus,
                                    indexfile=self.hisat_index,
                                    workdir=self.workdir,
                                    ref_file=self.ref_fasta,
                                    gff_file=os.path.join(self.workdir, "updated.gff"),
                                    kingdom=self.kingdom)],
                          local_scheduler=self.local_scheduler,
                          workers=self.no_of_jobs)

    def run_edger(self):
        build([DGE.edgeR(fastq_dic=self.fastq_dic,
                         num_cpus=self.num_cpus,
                         indexfile=self.hisat_index,
                         bindir=self.bindir,
                         workdir=self.workdir,
                         ref_file=self.ref_fasta,
                         kingdom=self.kingdom,
                         gff=os.path.join(self.workdir, "updated.gff"),
                         exp_design=self.exp_desn_file,
                         p_value=self.p_value)],
        local_scheduler=self.local_scheduler, workers=self.no_of_jobs)

    def run_deseq2(self):
        build([DGE.DESeq2(fastq_dic=self.fastq_dic,
                          num_cpus=self.num_cpus,
                          indexfile=self.hisat_index,
                          bindir=self.bindir,
                          workdir=self.workdir,
                          ref_file=self.ref_fasta,
                          kingdom=self.kingdom,
                          gff=os.path.join(self.workdir, "updated.gff"),
                          exp_design=self.exp_desn_file,
                          p_value=self.p_value)],
        local_scheduler=self.local_scheduler, workers=self.no_of_jobs)

    def merge_stringtie(self):
        build([Summ.MergeStringTies(fastq_dic=self.fastq_dic,
                                    num_cpus=self.num_cpus,
                                    indexfile=self.hisat_index,
                                    workdir=self.workdir,
                                    ref_file=self.ref_fasta,
                                    gff_file=os.path.join(self.workdir, "updated.gff"),
                                    kingdom=self.kingdom)],
                local_scheduler=self.local_scheduler, workers=self.no_of_jobs)

    def restringtie(self):
        build([Summ.ReStringTieScoresW(fastq_dic=self.fastq_dic,
                                       num_cpus=self.num_cpus,
                                       indexfile=self.hisat_index,
                                       workdir=self.workdir,
                                       ref_file=self.ref_fasta,
                                       gff_file=os.path.join(self.workdir, "updated.gff"),
                                       kingdom=self.kingdom)],
              local_scheduler=self.local_scheduler,
              workers=self.no_of_jobs)


    def run_ballgown(self):
        build([DGE.ballgown(fastq_dic=self.fastq_dic,
                          num_cpus=self.num_cpus,
                          indexfile=self.hisat_index,
                          bindir=self.bindir,
                          workdir=self.workdir,
                          ref_file=self.ref_fasta,
                          kingdom=self.kingdom,
                          exp_design=self.exp_desn_file,
                          gff_file=os.path.join(self.workdir, "updated.gff"),
                          p_value=self.p_value)],
        local_scheduler=self.local_scheduler, workers=self.no_of_jobs)



