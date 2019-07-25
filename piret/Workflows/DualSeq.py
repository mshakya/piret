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
from piret.Checks.Design import CheckDesign
from piret.Runs import FaQC, Map, Summ, DGE, srna
from luigi.interface import build



class DualSeq:
    """Class that pieces luigi task for dual seq."""

    def __init__(self, fastq_dic, prok_org_code, euk_org_code, ref_fastas, ref_gffs, num_cpus,
                 local_scheduler, hisat_index, workdir, kingdom,
                 no_of_jobs, exp_desn_file,
                 p_value, stardb_dir, **kwargs):
        self.aligner = kwargs['aligner']
        self.ref_fastas = ref_fastas
        self.ref_gffs = ref_gffs
        self.fastq_dic = fastq_dic
        self.num_cpus = num_cpus
        self.hisat_index = hisat_index
        self.workdir = workdir
        self.kingdom = kingdom
        self.local_scheduler = local_scheduler
        self.no_of_jobs = no_of_jobs
        self.exp_desn_file = exp_desn_file
        self.p_value = p_value
        self.prok_org_code = prok_org_code
        self.euk_org_code = euk_org_code
        self.stardb_dir = stardb_dir

    def create_db(self, gff):
        """Function to create hisat index."""
        if self.aligner == "hisat2":
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
        elif self.aligner == "STAR":
            build([Map.STARindex(fasta=self.ref_fastas,
                             num_cpus=self.num_cpus,
                             gff_file=gff,
                             stardb_dir = self.stardb_dir,
                             kingdom=self.kingdom),
                   Map.SAMindex(fasta=self.ref_fastas, workdir=self.workdir),
                   Map.CreateSplice(gff_file=self.ref_gffs.split(",")[1],
                            workdir=self.workdir),
                   Map.GetChromName(prok_ref=self.ref_fastas.split(",")[0],
                            euk_ref=self.ref_fastas.split(",")[1],
                            workdir=self.workdir)],
                            local_scheduler=self.local_scheduler)
            
    def map_reads(self):
        """Function to map reads."""
        if self.aligner == "hisat2":
            build([Map.HisatMapW(fastq_dic=self.fastq_dic, num_cpus=self.num_cpus,
                                 indexfile=self.hisat_index, workdir=self.workdir)],
              local_scheduler=self.local_scheduler)
        elif self.aligner == "STAR":
            build([Map.map_starW(fastq_dic=self.fastq_dic, num_cpus=self.num_cpus,
                                 stardb_dir=self.stardb_dir, workdir=self.workdir)],
                  local_scheduler=self.local_scheduler)

    def map_summarize(self):
        """Summarize mapped reads into a table."""
        if self.aligner == "hisat2":
            build([Map.SummarizeHisatMap(fastq_dic=self.fastq_dic,
                                    workdir=self.workdir,
                                    indexfile=self.hisat_index,
                                    num_cpus=self.num_cpus)],
            local_scheduler=self.local_scheduler, workers=1)
        elif self.aligner == "STAR":
            build([Map.SummarizeStarMap(fastq_dic=self.fastq_dic,
                                    workdir=self.workdir,
                                    stardb_dir=self.stardb_dir,
                                    num_cpus=self.num_cpus)],
            local_scheduler=self.local_scheduler, workers=1)

    def split_prokeuk(self):
        build([Map.Split2ProkEukW(fastq_dic=self.fastq_dic,
                                workdir=self.workdir,
                                  ref_file=self.ref_fastas)],
                          local_scheduler=self.local_scheduler,
                          workers=self.no_of_jobs)

    def summarize_map(self):
        """Summarize mapped reads into a table."""
        if self.aligner == "hisat2":
            build([Map.SummarizeHisatMap(fastq_dic=self.fastq_dic,
                                    workdir=self.workdir,
                                    indexfile=self.hisat_index,
                                    num_cpus=self.num_cpus)],
            local_scheduler=self.local_scheduler, workers=1)
        elif self.aligner == "STAR":
            build([Map.SummarizeStarMap(fastq_dic=self.fastq_dic,
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

    def novel_regions(self, gff):
        """Find novel regions."""
        build([srna.FindNovelRegionsW(fastq_dic=self.fastq_dic,
                                  workdir=self.workdir,
                                 kingdom=self.kingdom, 
                                 gff_file=gff)],
              local_scheduler=self.local_scheduler)
    
    def create_new_gff(self, gff):
        build([srna.CompileGFF(fastq_dic=self.fastq_dic,
                               kingdom=self.kingdom,
                               workdir=self.workdir,
                               gff_file=gff)],
          local_scheduler=self.local_scheduler, workers=self.no_of_jobs)

    def merge_stringties(self):
        build([Summ.MergeStringTies(fastq_dic=self.fastq_dic, num_cpus=self.num_cpus,
                                indexfile=self.hisat_index, workdir=self.workdir,
                                gff_file=self.ref_gffs,
                                kingdom="both")],
           local_scheduler=self.local_scheduler, workers=self.no_of_jobs)
    
    def restringtie(self):
        build([Summ.ReStringTieScoresW(fastq_dic=self.fastq_dic,
                                       num_cpus=self.num_cpus,
                                       workdir=self.workdir,
                                       kingdom=self.kingdom)],
           local_scheduler=self.local_scheduler, workers=self.no_of_jobs
           )

    def run_ballgown(self):
        build([DGE.ballgown(fastq_dic=self.fastq_dic, num_cpus=self.num_cpus,
                            workdir=self.workdir,
                            kingdom="eukarya",
                            exp_design=self.exp_desn_file,
                            p_value=self.p_value)],
               local_scheduler=self.local_scheduler, workers=self.no_of_jobs)

    def feature_counts(self):
        build([Summ.FeatureCounts(fastq_dic=self.fastq_dic,
                                  num_cpus=self.num_cpus,
                                  gff_file=os.path.join(self.workdir, "prok_updated.gff"),
                                  indexfile=self.hisat_index,
                                  kingdom='prokarya',
                                  workdir=self.workdir,
                                  ref_file=self.ref_fastas.split(",")[0]),
          Summ.FeatureCounts(fastq_dic=self.fastq_dic,
                             num_cpus=self.num_cpus,
                             gff_file=os.path.join(self.workdir, "euk_updated.gff"),
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
                   gff_file=os.path.join(self.workdir, "prok_updated.gff"),
                   exp_design=self.exp_desn_file,
                   p_value=self.p_value,
                   prok_org_code=self.prok_org_code,
                   euk_org_code=None),
               DGE.edgeR(fastq_dic=self.fastq_dic, num_cpus=self.num_cpus,
                    indexfile=self.hisat_index,
                    workdir=self.workdir,
                    ref_file=self.ref_fastas.split(",")[1],
                    kingdom='eukarya',
                    gff_file=os.path.join(self.workdir, "euk_updated.gff"),
                    exp_design=self.exp_desn_file,
                    p_value=self.p_value,
                    euk_org_code=self.euk_org_code,
                    prok_org_code=None)], local_scheduler=self.local_scheduler,
                          workers=self.no_of_jobs)

    def run_deseq2(self):
        build([DGE.DESeq2(fastq_dic=self.fastq_dic,
                          num_cpus=self.num_cpus,
                          indexfile=self.hisat_index,
                          workdir=self.workdir,
                          ref_file=self.ref_fastas.split(",")[0],
                          kingdom='prokarya',
                          gff_file=os.path.join(self.workdir, "prok_updated.gff"),
                          exp_design=self.exp_desn_file,
                          p_value=self.p_value,
                          prok_org_code=self.prok_org_code,
                          euk_org_code=None),
                DGE.DESeq2(fastq_dic=self.fastq_dic,
                            num_cpus=self.num_cpus,
                            indexfile=self.hisat_index,
                            workdir=self.workdir,
                            ref_file=self.ref_fastas.split(",")[1],
                            kingdom='eukarya',
                            gff_file=os.path.join(self.workdir, "euk_updated.gff"),
                            exp_design=self.exp_desn_file,
                            p_value=self.p_value,
                            euk_org_code=self.euk_org_code,
                            prok_org_code=None)],
                local_scheduler=self.local_scheduler, workers=self.no_of_jobs)


