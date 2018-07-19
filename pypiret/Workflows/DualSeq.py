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



class DualSeq:
    """Class that pieces luigi task for dual seq."""

    def __init__(self, fastq_dic, ref_fastas, ref_gffs, num_cpus,
                 local_scheduler, hisat_index, workdir, kingdom,
                 no_of_jobs, bindir, exp_desn_file,
                 p_value):
        self.ref_fastas = ref_fastas
        self.ref_gffs = ref_gffs
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
        build([Map.HisatIndex(fasta=self.ref_fastas,
                              hi_index=self.hisat_index,
                              num_cpus=self.num_cpus),
               Map.SAMindex(fasta=self.ref_fastas, workdir=self.workdir),
               Map.CreateSplice(gff_file=self.ref_gffs.split(",")[1],
                            workdir=self.workdir),
               Map.GetChromName(prok_ref=self.ref_fastas.split(",")[0],
                            euk_ref=self.ref_fastas.split(",")[1],
                            workdir=self.workdir)],
                          local_scheduler=self.local_scheduler,
                          workers=self.no_of_jobs)

    def map_hisat(self):
        """Function that maps reads to hisat2 index."""
        build([Map.HisatMapW(fastq_dic=self.fastq_dic,
                                 num_cpus=self.num_cpus,
                                 indexfile=self.hisat_index,
                                 workdir=self.workdir,
                                 ref_file=self.ref_fastas)],
                   local_scheduler=self.local_scheduler,
                   workers=self.no_of_jobs)
            

    def split_prokeuk(self):
        build([Map.Split2ProkEukW(fastq_dic=self.fastq_dic,
                                   num_cpus=self.num_cpus,
                                   indexfile=self.hisat_index,
                                   workdir=self.workdir,
                                   ref_file=self.ref_fastas)],
                          local_scheduler=self.local_scheduler,
                          workers=self.no_of_jobs)

    def summarize_map(self):
        build([srna.SummarizeMap(fastq_dic=self.fastq_dic,
                                 num_cpus=self.num_cpus,
                                 indexfile=self.hisat_index,
                                 workdir=self.workdir,
                                 ref_file=self.ref_fastas,
                                 kingdom=self.kingdom)],
                          local_scheduler=self.local_scheduler,
                          workers=self.no_of_jobs)

    def novel_regions(self):
        build([srna.FindNovelRegionsW(fastq_dic=self.fastq_dic,
                                      num_cpus=self.num_cpus,
                                      indexfile=self.hisat_index,
                                      workdir=self.workdir,
                                      ref_file=self.ref_fastas,
                                      gff_file=self.ref_gffs,
                                      kingdom=self.kingdom)],
                          local_scheduler=self.local_scheduler,
                          workers=self.no_of_jobs)
    
    def compile_gff(self):
        build([srna.CompileGFF(fastq_dic=self.fastq_dic,
                               num_cpus=self.num_cpus,
                               indexfile=self.hisat_index,
                               kingdom=self.kingdom,
                               workdir=self.workdir,
                               ref_file=self.ref_fastas,
                               gff_file=self.ref_gffs)],
          local_scheduler=self.local_scheduler, workers=self.no_of_jobs)

    def merge_stringties(self):
        build([Summ.MergeStringTies(fastq_dic=self.fastq_dic, num_cpus=self.num_cpus,
                                indexfile=self.hisat_index, workdir=self.workdir,
                                ref_file=self.ref_fastas,
                                gff_file=self.ref_gffs,
                                kingdom="both")],
           local_scheduler=self.local_scheduler, workers=self.no_of_jobs)
    
    def restringtie(self):
        build([Summ.ReStringTieScoresW(fastq_dic=self.fastq_dic,
                                       num_cpus=self.num_cpus,
                                       indexfile=self.hisat_index,
                                       workdir=self.workdir,
                                       ref_file=self.ref_fastas,
                                       gff_file=self.ref_gffs,
                                       kingdom=self.kingdom)],
           local_scheduler=self.local_scheduler, workers=self.no_of_jobs
           )

    def run_ballgown(self):
        build([DGE.ballgown(fastq_dic=self.fastq_dic, num_cpus=self.num_cpus,
                            indexfile=self.hisat_index, workdir=self.workdir,
                            ref_file=self.ref_fastas, gff_file=self.ref_gffs,
                            kingdom="eukarya",
                            exp_design=self.exp_desn_file,
                            bindir=self.bindir, p_value=self.p_value)],
               local_scheduler=self.local_scheduler, workers=self.no_of_jobs)

    def feature_counts(self):
        build([Summ.FeatureCounts(fastq_dic=self.fastq_dic,
                                  num_cpus=self.num_cpus,
                                  gff=os.path.join(self.workdir, "prok_updated.gff"),
                                  indexfile=self.hisat_index,
                                  kingdom='prokarya',
                                  workdir=self.workdir,
                                  ref_file=self.ref_fastas.split(",")[0]),
          Summ.FeatureCounts(fastq_dic=self.fastq_dic,
                             num_cpus=self.num_cpus,
                             gff=os.path.join(self.workdir, "euk_updated.gff"),
                             indexfile=self.hisat_index,
                             kingdom='eukarya',
                             workdir=self.workdir,
                             ref_file=self.ref_fastas.split(",")[1])],
          local_scheduler=self.local_scheduler, workers=self.no_of_jobs)


    def run_edger(self):

        build([DGE.edgeR(fastq_dic=self.fastq_dic,
                   num_cpus=self.num_cpus,
                   indexfile=self.hisat_index,
                   workdir=self.workdir,
                   ref_file=self.ref_fastas.split(",")[0],
                   kingdom='prokarya',
                   gff=os.path.join(self.workdir, "prok_updated.gff"),
                   exp_design=self.exp_desn_file,
                   bindir=self.bindir,
                   p_value=self.p_value),
               DGE.edgeR(fastq_dic=self.fastq_dic, num_cpus=self.num_cpus,
                    indexfile=self.hisat_index,
                    bindir=self.bindir,
                    workdir=self.workdir,
                    ref_file=self.ref_fastas.split(",")[1],
                    kingdom='eukarya',
                    gff=os.path.join(self.workdir, "euk_updated.gff"),
                    exp_design=self.exp_desn_file,
                    p_value=self.p_value)], local_scheduler=self.local_scheduler,
                          workers=self.no_of_jobs)

    def run_deseq2(self):
        build([DGE.DESeq2(fastq_dic=self.fastq_dic,
                                      num_cpus=self.num_cpus,
                                      indexfile=self.hisat_index,
                                      bindir=self.bindir,
                                      workdir=self.workdir,
                                      ref_file=self.ref_fastas.split(",")[0],
                                      kingdom='prokarya',
                                      gff=os.path.join(self.workdir, "prok_updated.gff"),
                                      exp_design=self.exp_desn_file,
                                      p_value=self.p_value),
                DGE.DESeq2(fastq_dic=self.fastq_dic,
                                      num_cpus=self.num_cpus,
                                      indexfile=self.hisat_index,
                                      bindir=self.bindir,
                                      workdir=self.workdir,
                                      ref_file=self.ref_fastas.split(",")[1],
                                      kingdom='eukarya',
                                      gff=os.path.join(self.workdir, "euk_updated.gff"),
                                      exp_design=self.exp_desn_file,
                                      p_value=self.p_value)],
                          local_scheduler=self.local_scheduler,
                          workers=self.no_of_jobs)


